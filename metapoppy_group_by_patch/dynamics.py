import epyc
import numpy
import math
from network import *


class Dynamics(epyc.Experiment, object):
    INITIAL_PATCHES = 'initial_patches'
    INITIAL_EDGES = 'initial_edges'
    INITIAL_TIME = 'initial_time'

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
        self._graphPrototype = g
        self._graphPrototype.prepare()

        self._create_events(self._graphPrototype)

        self._graph = None
        self._maxTime = self.DEFAULT_MAX_TIME

    def _create_events(self, network):
        raise NotImplementedError

    def network(self):
        """
        The network this set of dynamics runs upon
        :return:
        """
        return self._graph

    def set_maximum_time(self, t):
        """
        Set the maximum simulated time
        :param t: The maximum time
        """
        self._maxTime = t

    def _at_equilibrium(self, t):
        """
        Function to end simulation. Default is when max time is exceeded, can be overriden to end on a certain
        condition.
        :param t: Current simulated time
        :return: True to finish simulation
        """
        return t >= self._maxTime

    def setUp(self, params):
        """
        Configure the dynamics and the network ready for a simulation.
        :param params: Initial parameters - initial conditions of network and event reaction parameters.
        :return:
        """
        # Default setup
        epyc.Experiment.setUp(self, params)

        # Make a copy of the network prototype
        self._graph = self._graphPrototype.copy()

        # Seed the network patches
        for patch_id, data in params[Dynamics.INITIAL_PATCHES].iteritems():
            self._graph.update_patch(patch_id, data[Network.COMPARTMENTS], data[Network.ATTRIBUTES])

        # Seed the network edges
        for (u,v), data in params[Dynamics.INITIAL_EDGES].iteritems():
            self._graph.update_edge(u,v,data)

        # Set the event reactions parameters using the initial conditions
        self._seed_events(params[Network.EVENTS])

    def _seed_events(self, event_parameters):
        """
        Set the event reaction parameters
        :param params:
        :return:
        """
        raise NotImplementedError

    def do(self, params):
        """
        Run a metapoppy simulation. Uses Gillespie simulation - all events are given a rate based on the state of the
        network. An event is chosen and performed, the network is updated, time is incremented and new rates calculated.
        :param params:
        :return:
        """

        # Allow for designated start time (used for time/age specific events)
        if Dynamics.INITIAL_TIME in params:
            time = params[Dynamics.INITIAL_TIME]
        else:
            time = Dynamics.DEFAULT_START_TIME

        while not self._at_equilibrium(time):
            # TODO - this is inefficient, should post changes to dynamics
            # Get the total rate by summing rates at all patches
            total_network_rate = sum([d[Network.TOTAL_RATE] for _, d in self._graph.nodes(data=True)])

            # No events can occur, so end
            if total_network_rate == 0:
                break

            # Calculate the timestep delta
            dt = (1.0 / total_network_rate) * math.log(1.0 / numpy.random.random())

            # Choose an event
            r = numpy.random.random() * total_network_rate
            count = 0

            # Loop through all patches adding the total rate at each patch
            for n, d in self._graph.nodes(data=True):
                # If the total rate at this patch would exceed the r number, this is where an event will occur
                if count + d[Network.TOTAL_RATE] > r:
                    # Find the event which tips the count over r
                    for e in d[Network.EVENTS]:
                        count += e.rate()
                        if count >= r:
                            # Perform the event and stop iterations
                            e.perform(self._graph)
                            # Stop looping through events
                            break
                    # Stop looping through patches
                    break
                # Have not exceeded threshold, so continue
                else:
                    count += d[Network.TOTAL_RATE]
            # Move simulated time forward
            time += dt

        # TODO - experimental results. Maybe need to track populations over time
        return self.get_results()

    def get_results(self):
        return None

    def tearDown(self):
        """
        After a simulation ends, remove the used graph.
        :return:
        """
        # Perform the default tear-down
        epyc.Experiment.tearDown(self)

        # Remove the worked-on model
        self._graph = None
