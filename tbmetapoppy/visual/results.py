import networkx
import matplotlib.pyplot as plt

from metapoppy.visual.results import MetapoppyResults


class TBMetapoppyResults(MetapoppyResults):
    def __init__(self):
        MetapoppyResults.__init__(self)
        self._boundary = None

    def set_boundary(self, boundary):
        self._boundary = boundary

    def setup_frame(self, ax, timestep):
        boundary = plt.Polygon(self._boundary, fc='grey', alpha=0.3)
        plt.gca().add_patch(boundary)

    def draw_active_patch(self, G, patch, data):
        n = data['compartments']['c']
        networkx.draw_networkx_nodes(G, self._pos, nodelist=[patch], node_color=(1.0, 0.0, 0.0), node_size=10)

    def draw_inactive_patches(self, G, patches):
        # networkx.draw_networkx_nodes(G, self._pos, nodelist=[p], node_color=(a, 0.0, 0.0), node_size=5)
        pass