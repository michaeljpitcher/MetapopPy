import epyc
import math
from environment import *
import copy
import numpy
import itertools
import heapq

LAMBDA = lambda: 0


class Dynamics(epyc.Experiment, object):
    """
    MetapopPy Dynamics. Creates the dynamics that will occur over a networked metapopulation. Essentially a list of
    event and their respective rates of occurrence for each patch (node) on the network.

    For each parameter sample, an empty network is created and seeding calculated. For each repetition of the parameter
    sample, the network is set to an empty state and then seeded using the seed values. Patches can be determined to be
    "active" based on their values.

    All events calculate their rates for the active patches based on the patch contents. An event/patch combination is
    chosen probabilistically based on the individual rates, and a time-step for the event to occur if chosen based on
    the total rates, as per the Gillespie Algorithm. The event is performed, the network is updated and the updates are
    propagated to recalculate the rates of events. The process continues until a time limit is reached.
    """

    INITIAL_TIME = 'initial_time'
    MAX_TIME = 'max_time'

    EVENTS = 'events'

    # the default maximum simulation time
    DEFAULT_MAX_TIME = 100.0  #: Default maximum simulation time.
    DEFAULT_START_TIME = 0.0
    DEFAULT_RESULT_INTERVAL = 1.0

    def __init__(self, network=None):
        """
        Create MetapopPy dynamics to run over the given network.
        :param network: Network on which to run dynamics
        """
        epyc.Experiment.__init__(self)

        # Initialise variables
        self._network = self._rate_table = self._patch_seeding = self._edge_seeding =None

        self._row_for_patch = {}
        self._active_patches = []

        # Create the events
        self._events = self._create_events()
        assert self._events, "No events created"

        # Create dependency matrices
        self._comp_dependencies = {c: [] for c in network.compartments()}
        self._patch_att_dependencies = {a: [] for a in network.patch_attributes()}
        self._edge_att_dependencies = {a: [] for a in network.edge_attributes()}

        for col in range(len(self._events)):
            event = self._events[col]
            for c in event.get_dependent_compartments():
                self._comp_dependencies[c].append(col)
            for a in event.get_dependent_patch_attributes():
                self._patch_att_dependencies[a].append(col)
            for a in event.get_dependent_edge_attributes():
                self._edge_att_dependencies[a].append(col)

        # Set the network prototype if one has been provided - this will be the network used for all runs.
        # If not provided, a network must be created during configure stage.
        self._prototype_network = network

        # Default start and end times
        self._start_time = self.DEFAULT_START_TIME
        self._max_time = self.DEFAULT_MAX_TIME
        self._record_interval = self.DEFAULT_RESULT_INTERVAL

        # Posted events - will occur at set times
        self._posted_events = []

    def _create_events(self):
        """
        Create the events
        :return:
        """
        raise NotImplementedError

    def required_event_parameters(self):
        """
        All parameters which are required by the model
        :return:
        """
        params = []
        for e in self._events:
            params += e.parameter_keys()
        # Remove duplicates
        return list(set(params))

    def set_start_time(self, start_time):
        """
        Set the initial time of simulations
        :param start_time:
        :return:
        """
        self._start_time = float(start_time)

    def set_maximum_time(self, maximum_time):
        """
        Set the maximum simulated run-time
        :param maximum_time:
        :return:
        """
        assert maximum_time > self._start_time, "Maximum time must exceed start time"
        self._max_time = maximum_time

    def set_record_interval(self, record_interval):
        """
        Set the interval to record results
        :param record_interval:
        :return:
        """
        self._record_interval = record_interval

    def network(self):
        """
        The current state of the network this set of dynamics is running upon
        :return:
        """
        return self._network

    def configure(self, params):
        """
        Configure each param sample - runs once for every sample of parameter values. Creates a network (if needed),
        attaches the update propagation mechanism and obtains the initial conditions (seedings) required for each run.
        :param params:
        :return:
        """
        epyc.Experiment.configure(self, params)

        # No prototype provided so we must build a new network from the parameters
        if not self._prototype_network:
            self._network = self._build_network(params)
        # If a prototype has been provided and hasn't been set yet
        elif not self._network:
            self._network = self._prototype_network

        assert isinstance(self._network, Environment), "Graph must be instance of MetapopPy Network class"
        assert self._network.nodes(), "Empty network is invalid"

        # Attach the update handler to the network
        self._network.set_handlers(lambda p, c, a: self._propagate_patch_update(p, c, a),
                                   lambda u, v, a: self._propagate_edge_update(u, v, a))

        # Get the initial conditions
        self._patch_seeding = self._get_initial_patch_seeding(params)
        self._edge_seeding = self._get_initial_edge_seeding(params)

        # Configure events
        for e in self._events:
            e.set_parameters(params)

        # Configure time
        if Dynamics.INITIAL_TIME in params:
            self._start_time = params[Dynamics.INITIAL_TIME]

    def _build_network(self, params):
        raise NotImplementedError

    def _get_initial_patch_seeding(self, params):
        raise NotImplementedError

    def _get_initial_edge_seeding(self, params):
        raise NotImplementedError

    def setUp(self, params):
        """
        Set up the run for each repetition. Runs once for every repetition within a parameter sample. Resets the network
        to an empty state and then seeds it with the already-calculated values.
        :param params:
        :return:
        """

        # Default setup
        epyc.Experiment.setUp(self, params)

        # TODO - resetting the network only works on the assumption that the network structure (edges) has not changed
        # Reset the network
        self._network.reset()

        # Seed the network using the pre-calculated seeding
        for n in self._network.nodes:
            # Patch has a seeding
            if self._patch_seeding and n in self._patch_seeding:
                seed = self._patch_seeding[n]
                if Environment.COMPARTMENTS in seed:
                    comp_seed = seed[Environment.COMPARTMENTS]
                else:
                    comp_seed = {}
                if Environment.ATTRIBUTES in seed:
                    att_seed = seed[Environment.ATTRIBUTES]
                else:
                    att_seed = {}
                self._network.update_patch(n, comp_seed, att_seed)
            # Patch does not have a seeding, need to check if it is active
            elif not n in self._active_patches and self._patch_is_active(n):
                self._activate_patch(n)

        if self._edge_seeding:
            for (u, v), seed in self._edge_seeding.iteritems():
                self._network.update_edge(u, v, seed)

        # Check that at least one patch is active
        assert any(self._active_patches), "No patches are active"


    def _propagate_patch_update(self, patch_id, compartment_changes, patch_attribute_changes):
        """
        When a patch is changed, update the relevant entries in the rate table. This function is passed as a lambda
        function to the network, and is called whenever a change is made.
        :param patch_id:
        :param compartment_changes:
        :param patch_attribute_changes:
        :return:
        """
        # Check if patch is active
        # TODO - more efficient way of checking if patch is active
        active = patch_id in self._active_patches
        # If patch is already active
        if active:
            row = self._row_for_patch[patch_id]
            # Determine columns (events) to update by finding events which have dependencies on the items changed
            cols_to_update = set(itertools.chain(*[self._comp_dependencies[c] for c in compartment_changes] +
                                                  [self._patch_att_dependencies[c] for c in patch_attribute_changes]))
            for col in cols_to_update:
                event = self._events[col]
                self._rate_table[row][col] = event.calculate_rate_at_patch(self._network, patch_id)
        # Patch is not previously active but should become active from this update
        elif self._patch_is_active(patch_id):
            self._activate_patch(patch_id)

    def _propagate_edge_update(self, patch_u, patch_v, edge_attribute_changes):
        """
        If an edge has its attributes changed, propagate the update to events occurring at either end of the edge
        that are dependent on the edge attributes changed
        :param patch_u: 
        :param patch_v: 
        :param edge_attribute_changes:
        :return: 
        """
        for patch_id in [patch_u, patch_v]:
            # Check if patch is active
            # TODO - more efficient way of checking if patch is active
            active = patch_id in self._active_patches
            # If patch is already active
            if active:
                row = self._row_for_patch[patch_id]
                # Determine columns (events) to update by finding events which have dependencies on the items changed
                cols_to_update = set(itertools.chain(*[self._edge_att_dependencies[a] for a in edge_attribute_changes]))
                for col in cols_to_update:
                    event = self._events[col]
                    self._rate_table[row][col] = event.calculate_rate_at_patch(self._network, patch_id)
            # Patch is not previously active but should become active from this update
            elif self._patch_is_active(patch_id):
                self._activate_patch(patch_id)

    def update_parameter(self, parameter, change):
        """
        Change the value of a parameter (e.g. if time-dependent). Will update the relevant columns of the rate table
        for all events which depend on this parameter.
        :param parameter:
        :param change:
        :return:
        """
        params = self.parameters()
        params[parameter] += change

        # TODO: could be better with a parameter dependency table
        # Loop through all events
        for col in range(self._rate_table.shape[1]):
            event = self._events[col]
            # Check if changed parameter is needed by the event
            if parameter in event.parameter_keys():
                # Update the parameter value on the event
                event.set_parameters(self._parameters)
                for row in range(self._rate_table.shape[0]):
                    # Recalculate the event rate at every patch
                    patch_id = self._active_patches[row]
                    self._rate_table[row][col] = event.calculate_rate_at_patch(self._network, patch_id)

    def _patch_is_active(self, patch_id):
        """
        Determine if the given patch is active (from the network). Default is that patches are always active, can be
        overridden to only process patches based on a given condition.
        :param patch_id:
        :return:
        """
        return True

    def _activate_patch(self, patch_id):
        """
        A patch has become active, so create a new row in the rate table for it and determine rates of events there.
        :param patch_id:
        :return:
        """
        # Calculate the row number
        self._row_for_patch[patch_id] = len(self._active_patches)
        # Add to active patch list
        self._active_patches.append(patch_id)

        # Create a row of rates - value in each column is rate of an event at this patch
        rates = numpy.array([e.calculate_rate_at_patch(self._network, patch_id) for e in self._events]).reshape(1, len(
            self._events))
        if self._rate_table is None:
            # This row becomes the rate table if it's currently empty
            self._rate_table = rates
        else:
            # Add the rates for this patch as a new row (build a new table by concatenation)
            self._rate_table = numpy.concatenate((self._rate_table, rates), 0)

        # Patch is activated, so seed it
        # Get seeding
        seeding = self._seed_activated_patch(patch_id, self.parameters())
        if Environment.COMPARTMENTS in seeding:
            comp_seeding = seeding[Environment.COMPARTMENTS]
        else:
            comp_seeding = None
        if Environment.ATTRIBUTES in seeding:
            att_seeding = seeding[Environment.ATTRIBUTES]
        else:
            att_seeding = None
        # Update the patch with the seeding values
        self._network.update_patch(patch_id, comp_seeding, att_seeding)

    def _seed_activated_patch(self, patch_id, params):
        raise NotImplementedError

    def post_event(self, t, event):
        """
        Post a event to occur at a set time at a given patch
        :param t:
        :param event:
        :param patch_id:
        :return:
        """
        assert isinstance(event, type(LAMBDA)) and event.__name__ == LAMBDA.__name__, \
            "Posted event must be a lambda function"
        heapq.heappush(self._posted_events, (t, event))

    def do(self, params):
        """
        Run a MetapopPy simulation. Uses Gillespie simulation - all combinations of events and patches are given a rate
        based on the state of the network. An event and patch combination are chosen and performed, the event is
        performed, updating the patch (and others). Time is incremented (based on total rates) and new rates calculated.
        :param params:
        :return:
        """
        results = {}

        time = self._start_time

        def record_results(results, record_time):
            # print "t=", record_time
            # TODO - we don't record edges / non-active patches
            results[record_time] = {}
            for p in self._active_patches:
                results[record_time][p] = copy.deepcopy(self._network.node[p])
            return results

        results = record_results(results, time)
        # Avoid rounding issues with time interval by rounding to 7 decimal places
        next_record_interval = round(time + self._record_interval, 7)

        total_network_rate = numpy.sum(self._rate_table)
        assert total_network_rate, "No events possible at start of simulation"

        num_events = len(self._events)

        while not self._at_equilibrium(time):
            # Calculate the timestep delta
            dt = (1.0 / total_network_rate) * math.log(1.0 / numpy.random.random())

            # If there's posted events scheduled to occur before the next event occurs - pick one and process it
            if self._posted_events:
                next_time = self._posted_events[0][0]
                if time + dt > next_time:
                    next_time, next_event = heapq.heappop(self._posted_events)
                    # Perform the event
                    next_event()
                    # Time progresses to the time of posted event
                    time = next_time
                    # Event has been executed, go to next loop
                    # NOTE: cannot continue processing events as this event may have changed rates of dynamic events
                    continue

            # Choose an event and patch based on the values in the rate table
            # TODO - numpy multinomial is faster than numpy choice (in python 2, maybe not in 3?)
            index_choice = numpy.random.multinomial(1, self._rate_table.flatten() / total_network_rate).argmax()
            # Given the index, find the patch (row) and event (column) this refers to
            patch_id = self._active_patches[index_choice / num_events]
            event = self._events[index_choice % num_events]

            # Perform the event. Handler will propagate the effects of any network updates
            event.perform(self._network, patch_id)

            # Move simulated time forward
            time += dt

            # Record results if interval(s) exceeded
            while time >= next_record_interval and next_record_interval <= self._max_time:
                record_results(results, next_record_interval)
                # Avoid rounding issues
                next_record_interval = round(next_record_interval + self._record_interval, 7)

            # Get the total rate by summing rates of all events at all patches
            total_network_rate = numpy.sum(self._rate_table)

            # If no events can occur, then end
            if total_network_rate == 0:
                break

        return results

    def _at_equilibrium(self, t):
        """
        Function to end simulation. Default is when max time is exceeded, can be overriden to end on a certain
        condition.
        :param t: Current simulated time
        :return: True to finish simulation
        """
        return t >= self._max_time

    def tearDown( self ):
        """
        Finish a run for a repetition. Runs once for every repetition within a parameter sample. Resets all values ready
        for the next run.
        :return:
        """
        # Perform the default tear-down
        epyc.Experiment.tearDown(self)

        # Reset rate table and lookups
        self._rate_table = None
        self._active_patches = []
        self._row_for_patch = {}

        # Reset posted events
        self._posted_events = []
