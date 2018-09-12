import epyc
import numpy
import math
from network import *


class Dynamics(epyc.Experiment, object):
    INITIAL_PATCHES = 'initial_patches'
    INITIAL_EDGES = 'initial_edges'
    INITIAL_TIME = 'initial_time'
    MAX_TIME = 'max_time'

    EVENTS = 'events'

    # the default maximum simulation time
    DEFAULT_MAX_TIME = 20000  #: Default maximum simulation time.
    DEFAULT_START_TIME = 0

    def __init__(self, g):
        """
        Create MetapopPy dynamics to run over the given network.
        :param g: Network on which to run dynamics
        """
        assert isinstance(g, Network), "Graph must be instance of MetapopPy Network class"
        epyc.Experiment.__init__(self)
        self._compartments = g.compartments()
        self._patch_attributes = g.patch_attributes()
        self._edge_attributes = g.edge_attributes()

        # Prepare the network
        self.network_prototype = g
        self.network_prototype.prepare(lambda a, b, c, e: self._propagate_updates(a, b, c, e))

        self._events = []
        self._create_events(self.network_prototype)

        self._patch_ids = self.network_prototype.nodes()

        self._rate_table = numpy.array([[d, ] * len(self._patch_ids) for d in [0.0, ] * len(self._events)],
                                       dtype=numpy.float)

        self._network = None
        self._start_time = self.DEFAULT_START_TIME
        self._max_time = self.DEFAULT_MAX_TIME

    def _propagate_updates(self, patch_id, compartment_changes, attribute_changes, edge_changes):
        for e in range(len(self._events)):
            event = self._events[e]
            self._rate_table[e][patch_id] = event.calculate_rate_at_patch(self._network, patch_id)
        # TODO - only update events dependent on atts/comps changed

    def _create_events(self, network):
        raise NotImplementedError

    def network(self):
        """
        The network this set of dynamics runs upon
        :return:
        """
        return self._network

    def _at_equilibrium(self, t):
        """
        Function to end simulation. Default is when max time is exceeded, can be overriden to end on a certain
        condition.
        :param t: Current simulated time
        :return: True to finish simulation
        """
        return t >= self._max_time

    def _seed_events(self, event_parameters):
        """
        Set the event reaction parameters
        :param event_parameters:
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

        # Make a copy of the network prototype
        self._network = self.network_prototype.copy()

        # Seed the network patches
        for patch_id, data in params[Dynamics.INITIAL_PATCHES].iteritems():
            self._network.update_patch(patch_id, data[Network.COMPARTMENTS], data[Network.ATTRIBUTES])

        # Seed the network edges
        for (u,v), data in params[Dynamics.INITIAL_EDGES].iteritems():
            self._network.update_edge(u, v, data[Dynamics.INITIAL_EDGES])

        # Set the event reactions parameters using the initial conditions
        self._seed_events(params[Dynamics.EVENTS])

        # Allow for designated start time (used for time/age specific events)
        if Dynamics.INITIAL_TIME in params:
            self._start_time = params[Dynamics.INITIAL_TIME]

        # Set the maximum run time
        if Dynamics.MAX_TIME in params:
            self._max_time = params[Dynamics.MAX_TIME]

    def do(self, params):
        """
        Run a metapoppy simulation. Uses Gillespie simulation - all events are given a rate based on the state of the
        network. An event is chosen and performed, the network is updated, time is incremented and new rates calculated.
        :param params:
        :return:
        """

        time = self._start_time

        while not self._at_equilibrium(time):
            # Get the total rate by summing rates of all events at all patches
            total_network_rate = numpy.sum(self._rate_table)

            # If no events can occur, then end
            if total_network_rate == 0:
                break

            # Calculate the timestep delta
            dt = (1.0 / total_network_rate) * math.log(1.0 / numpy.random.random())

            # Choose an event and patch based on the values in the rate table
            index_choice = numpy.random.choice(range(0, self._rate_table.size),
                                               p=self._rate_table.flatten() / total_network_rate)
            event = self._events[index_choice / len(self._events)]
            patch = self._patch_ids[index_choice % len(self._patch_ids)]

            # Perform the event. Handler will propagate the effects of any network updates
            event.perform(self._network.node[patch], self._network.edges([patch]))

            # Move simulated time forward
            time += dt

        # TODO - experimental results. Maybe need to track populations over time
        return self.get_results()

    def get_results(self):
        raise NotImplementedError

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
        self._rate_table = numpy.array([[d, ] * len(self._patch_ids) for d in [0.0, ] * len(self._events)],
                                       dtype=numpy.float)
