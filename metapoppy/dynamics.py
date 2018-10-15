import epyc
import math
from network import *
import copy


class Dynamics(epyc.Experiment, object):
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
        self._network_prototype = self._rate_table = self._patch_for_column = self._column_for_patch = \
            self._network = None

        # Create the events
        self._events = self._create_events()
        assert self._events, "No events created"

        # Set the network prototype if one has been provided (if not provided, will be required during the configure
        # stage)
        if network:
            self.set_network_prototype(network)

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

    def set_network_prototype(self, prototype):
        """
        Set the network prototype and configure rate tables
        :param prototype:
        :return:
        """
        # Prepare the network
        self._network_prototype = prototype
        assert isinstance(prototype, Network), "Graph must be instance of MetapopPy Network class"
        assert self._network_prototype.nodes(), "Empty network is invalid"

        # Obtain the patches of the network (for lookup purposes)
        self._patch_for_column = self._network_prototype.nodes()[:]
        self._column_for_patch = dict([(self._patch_for_column[k], k) for k in range(len(self._patch_for_column))])

        # Create a rate table. Rows are events, columns are patches
        self._rate_table = numpy.zeros([len(self._events), len(self._patch_for_column)], dtype=numpy.float)

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
        epyc.Experiment.configure(self, params)

        # Set the event reaction parameters using the initial conditions
        for e in self._events:
            e.set_parameters(params)

        # Allow for designated start time (used for time/age specific events)
        if Dynamics.INITIAL_TIME in params:
            self._start_time = params[Dynamics.INITIAL_TIME]

        # Prepare network - reset all values to zero
        assert self._network_prototype, "Network has not been configured"
        self._network_prototype.prepare()

        # Seed the prototype network
        self._seed_prototype_network(params)

    def _seed_prototype_network(self, params):
        """
        Seed the network in some manner based on the parameters.
        :param params:
        :return:
        """
        raise NotImplementedError

    def setUp(self, params):
        """
        Configure the dynamics and the network ready for a simulation.
        :param params: Initial parameters - initial conditions of network and event reaction parameters.
        :return:
        """
        # Default setup
        epyc.Experiment.setUp(self, params)

        # Use a copy of the network prototype (must be done on every run in case network has changed)
        self._network = self._network_prototype.copy()

        # Attach the update handler
        self._network.set_handler(lambda a, b, c, d: self._propagate_updates(a, b, c, d))

        # Set initial rate values for all event/patch combinations
        for col in range(len(self._patch_for_column)):
            for row in range(len(self._events)):
                patch_id = self._patch_for_column[col]
                event = self._events[row]
                self._rate_table[row][col] = event.calculate_rate_at_patch(self._network, patch_id)

    def _propagate_updates(self, patch_id, compartment_changes, attribute_changes, edge_changes):
        """
        When a patch ID is changed, update the relevant entries in the rate table. This function is passed as a lambda
        function to the network, and is called whenever a change is made.
        :param patch_id:
        :param compartment_changes:
        :param attribute_changes:
        :param edge_changes:
        :return:
        """
        col = self._column_for_patch[patch_id]
        for row in range(len(self._events)):
            event = self._events[row]
            self._rate_table[row][col] = event.calculate_rate_at_patch(self._network, patch_id)
        # TODO - only update events dependent on atts/comps changed

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
            results["t=" + str(record_time)] = {}
            for p, data in self.network().nodes(data=True):
                results["t=" + str(record_time)][p] = copy.deepcopy(data)
            return results

        number_patches = len(self._patch_for_column)
        # indices = range(0, self._rate_table.size)

        results = record_results(results, time)
        next_record_interval = time + self._record_interval

        while not self._at_equilibrium(time):

            # Get the total rate by summing rates of all events at all patches
            total_network_rate = numpy.sum(self._rate_table)

            # If no events can occur, then end
            if total_network_rate == 0:
                break

            # Calculate the timestep delta
            dt = (1.0 / total_network_rate) * math.log(1.0 / numpy.random.random())

            # Choose an event and patch based on the values in the rate table
            # TODO - numpy multinomial is faster than numpy choice (in python 2, maybe not in 3?)
            # index_choice = numpy.random.choice(indices, p=self._rate_table.flatten() / total_network_rate)
            index_choice = numpy.random.multinomial(1, self._rate_table.flatten() / total_network_rate).argmax()
            # Given the index, find the event and patch this refers to
            event = self._events[index_choice / number_patches]
            patch_id = self._patch_for_column[index_choice % number_patches]

            # Perform the event. Handler will propagate the effects of any network updates
            event.perform(self._network, patch_id)
            # Move simulated time forward
            time += dt

            # Record results if interval(s) exceeded
            while time >= next_record_interval and next_record_interval <= self._max_time:
                record_results(results, next_record_interval)
                next_record_interval += self._record_interval

        return results

    def _at_equilibrium(self, t):
        """
        Function to end simulation. Default is when max time is exceeded, can be overriden to end on a certain
        condition.
        :param t: Current simulated time
        :return: True to finish simulation
        """
        return t >= self._max_time

    def tearDown(self):
        """
        After a simulation ends, remove the used graph.
        :return:
        """
        # Perform the default tear-down
        epyc.Experiment.tearDown(self)

        # Remove the worked-on model
        self._network = None

        # Reset rate table
        self._rate_table.fill(0.0)
