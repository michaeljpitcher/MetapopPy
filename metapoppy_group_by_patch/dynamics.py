import epyc
import numpy
import math
from network import *


class Dynamics(epyc.Experiment, object):
    INITIAL_PATCHES = 'initial_patches'
    INITIAL_EDGES = 'initial_edges'

    # the default maximum simulation time
    DEFAULT_MAX_TIME = 20000  #: Default maximum simulation time.

    def __init__(self, g):
        epyc.Experiment.__init__(self)
        self._compartments = g.compartments()
        self._patch_attributes = g.patch_attributes()
        self._edge_attributes = g.edge_attributes()

        self._graphPrototype = g
        self._graphPrototype.prepare()
        self._create_events(self._graphPrototype)

        self._graph = None
        self._maxTime = self.DEFAULT_MAX_TIME

    def _create_events(self, network):
        raise NotImplementedError

    def network(self):
        """
        The network this dynamics runs upon
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

    def _seed_events(self, params):
        raise NotImplementedError

    def do(self, params):
        """
        Run a metapoppy simulation. Uses Gillespie simulation - all events are given a rate based on the state of the
        network. An event is chosen and performed, the network is updated, time is incremented and new rates calculated.
        :param params:
        :return:
        """

        # run the dynamics
        time = 0

        while not self._at_equilibrium(time):
            total_network_rate = sum([d[Network.TOTAL_RATE] for _, d in self._graph.nodes(data=True)])

            if total_network_rate == 0:
                break

            # calculate the timestep delta
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
                            return
                else:
                    count += d[Network.TOTAL_RATE]
            time += dt

    def tearDown(self):
        """
        After a simulation ends, destroy the used graph.
        :return:
        """
        # perform the default tear-down
        epyc.Experiment.tearDown(self)

        # throw away the worked-on model
        self._graph = None
