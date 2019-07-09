import networkx
import matplotlib.animation
import matplotlib.pyplot as plt
from matplotlib import gridspec

from metapoppy.visual.results import MetapoppyResults


class TBMetapoppyResults(MetapoppyResults):
    def __init__(self):
        MetapoppyResults.__init__(self)
        self._boundary = None
        self.plot_vals = {}
        self._cav = {}

    def set_boundary(self, boundary):
        self._boundary = boundary

    def draw_active_patch(self, G, patch, data, ax):
        val = sum([data['compartments'][k] for k in ['b_ed', 'b_er', 'b_im', 'b_id']]) / 7000.0

        if val < 0.5:
            col = 'orange'
        else:
            col = 'red'
            if patch not in self._cav:
                self._cav[patch] = 1
            else:
                self._cav[patch] += 1

        if patch in self._cav:
            networkx.draw_networkx_nodes(G, self._pos, ax=ax, nodelist=[patch], node_color='black',
                                         node_size=(200 * val)+(5*self._cav[patch]))

        networkx.draw_networkx_nodes(G, self._pos, ax=ax, nodelist=[patch], node_color=col, node_size=200 * val)


    def draw_inactive_patches(self, G, patches, ax):
        # networkx.draw_networkx_nodes(G, self._pos, nodelist=[p], node_color=(a, 0.0, 0.0), node_size=5)
        pass

    def create_tbmetapoppy_network(self):
        self.plot_vals = {}
        param_sample = 0
        repetition = 0

        # Create Graph
        G = networkx.Graph()
        for k, v in self._pos.iteritems():
            G.add_node(k, Position=v, Active=False)

        # Build plot
        fig = plt.figure(figsize=(12, 10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])

        lung_rec = plt.Rectangle((0, 0), 50, 100, fc='0.9')
        w = 2.5
        h = 5
        lymph_rec = plt.Rectangle((60-w/2.0, 50-h/2.0), w, h, fc='0.9')

        def update(timestep_index):
            timestep = self.timesteps[timestep_index]

            ax0.clear()
            ax1.clear()

            ax0.add_patch(lung_rec)
            ax0.add_patch(lymph_rec)

            ax0.set_title('t={0}'.format(timestep))
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.text(25, 101, 'Lung')
            ax0.text(58.5, 53.5, 'Lymph')

            ax1.set_xlim(0, max(self.timesteps))
            ax1.set_ylim(0, 7000)
            ax1.set_ylabel('Bacteria')
            ax1.set_xlabel('Time (days)')

            data = self.results[param_sample][repetition][timestep]
            ax1.plot([0], [0], color='blue', label='Lung')
            for n in data.keys():
                G.nodes[n]['Active'] = True
                self.draw_active_patch(G, n, data[n], ax0)

                val = sum([data[n]['compartments'][k] for k in ['b_ed', 'b_er', 'b_im', 'b_id']])
                if n in self.plot_vals:
                    self.plot_vals[n].append((timestep, val))
                else:
                    self.plot_vals[n] = [(timestep, val)]
                if n == 'lymph_patch':
                    ax1.plot([t for t,_ in self.plot_vals[n]], [b for _,b in self.plot_vals[n]], color='green', label='Lymph')
                else:
                    ax1.plot([t for t, _ in self.plot_vals[n]], [b for _, b in self.plot_vals[n]], color='blue')
            ax1.plot([365,365],[0,7000], '--', linewidth=2, color='black')
            ax1.text(368, 6400, 'T-cell deficiency')

            ax1.legend()

        ani = matplotlib.animation.FuncAnimation(fig, update, frames=len(self.timesteps), interval=400, repeat=True)
        # TODO - save (issue with ffmpeg)
        plt.show()
