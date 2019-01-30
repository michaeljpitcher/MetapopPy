import epyc
import math
from network import *
import copy
import numpy
import itertools


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

    def _create_events(self):
        """
        Create the events
        :return:
        """
        raise NotImplementedError

    def required_event_parameters(self):
        params = []
        for e in self._events:
            params += e.parameter_keys()
        return list(set(params))

    def set_start_time(self, start_time):
        """
        Set the initial time of simulations
        :param start_time:
        :return:
        """
        self._start_time = start_time

    def set_maximum_time(self, maximum_time):
        """
        Set the maximum simulated run-time
        :param maximum_time:
        :return:
        """
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
        Configure each param sample. Run once for every sample of parameter values. Creates a network (if needed),
        attaches the update propagation mechanism and obtains the seedings required for each run.
        :param params:
        :return:
        """
        epyc.Experiment.configure(self, params)

        #  If a prototype has been provided (and hasn't already been set as the network i.e. first configure)
        if self._prototype_network and not self._network:
            self._network = self._prototype_network
        # No prototype provided so we must build a new network from the parameters
        elif not self._prototype_network:
            self._network = self._build_network(params)

        assert isinstance(self._network, Network), "Graph must be instance of MetapopPy Network class"
        assert self._network.nodes(), "Empty network is invalid"

        # Attach the update handler to the network
        self._network.set_handlers(lambda p, c, a: self._propagate_patch_update(p, c, a),
                                   lambda u, v, a: self._propagate_edge_update(u, v, a))

        self._patch_seeding = self._get_patch_seeding(params)
        self._edge_seeding = self._get_edge_seeding(params)

        # Configure events
        for e in self._events:
            e.set_parameters(params)

        # Configure time
        if Dynamics.INITIAL_TIME in params:
            self._start_time = params[Dynamics.INITIAL_TIME]

    def _build_network(self, params):
        raise NotImplementedError

    def _get_patch_seeding(self, params):
        raise NotImplementedError

    def _get_edge_seeding(self, params):
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

        # TODO - resetting the network only works on the assumption that the network structure has not changed
        # Reset the network
        self._network.reset()

        # Seed the network using the pre-calculated seeding
        for n, seed in self._patch_seeding.iteritems():
            if Network.COMPARTMENTS in seed:
                comp_seed = seed[Network.COMPARTMENTS]
            else:
                comp_seed = {}
            if Network.ATTRIBUTES in seed:
                att_seed = seed[Network.ATTRIBUTES]
            else:
                att_seed = {}
            self._network.update_patch(n, comp_seed, att_seed)
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

    def _propagate_edge_update(self, patch, neighbour, edge_attribute_changes):
        for patch_id in [patch, neighbour]:
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
