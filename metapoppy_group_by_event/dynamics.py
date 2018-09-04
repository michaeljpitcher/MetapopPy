import epyc
import numpy


class Dynamics(epyc.Experiment, object):
    # Additional metadata elements
    TIME = 'simulation_time'  #: Metadata element holding the logical simulation end-time.
    EVENTS = 'simulation_events'  #: Metadata element holding the number of events that happened.

    # the default maximum simulation time
    DEFAULT_MAX_TIME = 20000  #: Default maximum simulation time.

    def __init__(self, g):
        epyc.Experiment.__init__(self)
        self._graphPrototype = g
        self._graph = None
        self._maxTime = self.DEFAULT_MAX_TIME
        self._posted = []

    def network(self):
        '''Return the network this dynamics is running over.
        :returns: the network'''
        return self._graph

    def setMaximumTime(self, t):
        '''Set the maximum default simulation time. The default is given
        by :attr:`DEFAULT_MAX_TIME`.

        param: t: the maximum time'''
        self._maxTime = t

    def at_equilibrium(self, t):
        '''Test whether the model is an equilibrium. Override this method to provide
        alternative and/or faster simulations.

        :param t: the current simulation timestep
        :returns: True if we're done'''
        return (t >= self._maxTime)

    def setUp(self, params):
         # perform the default setup
        epyc.Experiment.setUp(self, params)

        # make a copy of the network prototype
        self._graph = self._graphPrototype.copy()

    def do(self, params):
        '''Run the simulation using Gillespie dynamics. The process terminates
        when either there are no events with zero rates or when :meth:`at_equilibrium`
        return True.

        :param params: the experimental parameters
        :returns: the experimental results dict'''

        # run the dynamics
        g = self.network()
        t = 0
        events = 0
        while not self.at_equilibrium(t):
            # pull the transition dynamics at this timestep
            transitions = self.eventRateDistribution(t)

            # compute the total rate of transitions for the entire network
            a = 0.0
            for (_, r, _) in transitions:
                a = a + r
            if a == 0:
                break  # no events with non-zero rates

            # calculate the timestep delta
            r1 = numpy.random.random()
            dt = (1.0 / a) * math.log(1.0 / r1)

            # calculate which event happens
            if len(transitions) == 1:
                # if there's only one, that's the one that happens
                (l, _, ef) = transitions[0]
            else:
                # otherwise, choose one at random based on the rates
                r2 = numpy.random.random()
                xc = r2 * a
                k = 0
                (l, xs, ef) = transitions[k]
                while xs < xc:
                    k = k + 1
                    (l, xsp, ef) = transitions[k]
                    xs = xs + xsp

            # increment the time
            t = t + dt

            # fire any events posted for at or before this time
            events = events + self.runPendingEvents(t)

            # it's possible that posted events have removed all elements
            # from the chosen locus, in which case we simply continue
            # with the next event selection
            # sd: is this correct? or does it mess up the statistics too much?
            if len(l) > 0:
                # draw a random element from the chosen locus
                e = l.draw()

                # perform the event by calling the event function,
                # passing the dynamics, event time, network, and element
                ef(self, t, g, e)

                # increment the event counter
                events = events + 1

        # run any events posted for before the maximum simulation time
        self.runPendingEvents(self._maxTime)

        # add some more metadata
        (self.metadata())[self.TIME] = t
        (self.metadata())[self.EVENTS] = events

        # report results
        rc = self.experimentalResults()
        return rc


    def tearDown(self):
        '''At the end of each experiment, throw away the copy.'''

        # perform the default tear-down
        epyc.Experiment.tearDown(self)

        # throw away the worked-on model
        self._graph = None

