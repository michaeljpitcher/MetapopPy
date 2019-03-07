import json
import epyc

import matplotlib.pyplot as plt
import matplotlib.animation

from ..network import *


class MetapoppyResults(object):
    def __init__(self):
        self._repetitions = 0
        self._parameters = []
        self.timesteps = []

        self.results = []
        self._parameter_samples = []

        self._pos = None

    def load_epyc_results_from_json(self, filename):
        """
        Given filename for an epyc JSON output file, loads the data from that file. Data can be aggregated for cases
        of multiple repetitions
        :param filename:
        :return:
        """
        with open(filename) as data_file:
            results = json.load(data_file)[epyc.Experiment.RESULTS]

        self._parameters = results.values()[0][0][epyc.Experiment.PARAMETERS].keys()
        # TODO - assumes same timesteps in every run
        self.timesteps = [float(n) for n in results.values()[0][0][epyc.Experiment.RESULTS].keys()]
        self.timesteps.sort()

        for param_sample_key, sample_results in results.iteritems():
            self._repetitions = len(sample_results)
            # Parameter values are held at repetition level so have to access the first repetition
            self._parameter_samples.append(sample_results[0][epyc.Experiment.PARAMETERS])
            self.results.append([{float(k): v for k, v in repetition[epyc.Experiment.RESULTS].iteritems()}
                                 for repetition in sample_results])

    def set_node_positions(self, pos):
        self._pos = pos

    def create_network(self, vis_folder):
        param_sample = 0
        repetition = 0

        assert self._pos, "Node positions not set"

        # Create Graph
        G = networkx.Graph()
        for k,v in self._pos.iteritems():
            G.add_node(k, Position=v, Active=False)

        # Build plot
        fig, ax = plt.subplots(figsize=(6, 4))

        def update(timestep_index):
            ax.clear()

            timestep = self.timesteps[timestep_index]

            # Scale plot ax
            ax.set_title("t={0}".format(timestep), fontweight="bold")
            ax.set_xticks([])
            ax.set_yticks([])

            self.setup_frame(ax, timestep)

            data = self.results[param_sample][repetition][timestep]
            for n in data.keys():
                G.nodes[n]['Active'] = True
                self.draw_active_patch(G, n, data[n])

            inactive_patches = [n for n in G.nodes if not G.nodes[n]['Active']]
            self.draw_inactive_patches(G, inactive_patches)

        ani = matplotlib.animation.FuncAnimation(fig, update, frames=len(self.timesteps), interval=1, repeat=True)
        # TODO - save (issue with ffmpeg)
        plt.show()

    def setup_frame(self, ax, timestep):
        pass

    def draw_active_patch(self, G, patch, data):
        networkx.draw_networkx_nodes(G, self._pos, nodelist=[patch], node_color='red', node_size=10)

    def draw_inactive_patches(self, G, patches):
        networkx.draw_networkx_nodes(G, self._pos, nodelist=patches, node_color='grey', node_size=5, alpha=0.1)
