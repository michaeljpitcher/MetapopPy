import networkx
import numpy


class Network(networkx.Graph):

    INITIAL = 'initial'
    COMPARTMENTS = 'compartments'
    ATTRIBUTES = 'attributes'
    EVENTS = 'events'
    TOTAL_RATE = 'total_rate'

    def __init__(self, compartments, patch_attributes, edge_attributes):
        self._compartments = compartments
        self._patch_attributes = patch_attributes
        self._edge_attributes = edge_attributes
        networkx.Graph.__init__(self)

    def compartments(self):
        return self._compartments

    def patch_attributes(self):
        return self._patch_attributes

    def edge_attributes(self):
        return self._edge_attributes

    def prepare(self):
        for n in self.nodes():
            self.node[n][Network.COMPARTMENTS] = dict([(c, 0) for c in self._compartments])
            self.node[n][Network.ATTRIBUTES] = dict([(c, 0) for c in self._patch_attributes])
            self.node[n][Network.EVENTS] = []
            self.node[n][Network.TOTAL_RATE] = 0.0
        for u, v in self.edges():
            for a in self._edge_attributes:
                self.edge[u][v][a] = 0.0

    # TODO - update the add node/edge functions with prepare (as above) to allow for dynamic graphs

    def get_compartment_value(self, patch_id, compartment):
        return self.node[patch_id][Network.COMPARTMENTS][compartment]

    def get_attribute_value(self, patch_id, attribute):
        return self.node[patch_id][Network.ATTRIBUTES][attribute]

    def attach_event(self, event, patch_id):
        event.set_patch_id(patch_id)
        self.node[patch_id][Network.EVENTS].append(event)

    def update_rates_at_patch(self, patch_id):
        patch_data = self.node[patch_id]
        for event in self.node[patch_id][Network.EVENTS]:
            rate_before = event.rate()
            event.update_rate(patch_data, [d for _,_,d in self.edges([patch_id], data=True)])
            rate_difference = event.rate() - rate_before
            patch_data[Network.TOTAL_RATE] += rate_difference

    def perform_an_event(self):
        total_network_rate = sum([d[Network.TOTAL_RATE] for _,d in self.nodes(data=True)])
        r = numpy.random.random() * total_network_rate
        count = 0
        # Loop through all patches adding the total rate at each patch
        for n,d in self.nodes(data=True):
            # If the total rate at this patch would exceed the r number, this is where an event will occur
            if count + d[Network.TOTAL_RATE] > r:
                # Find the event which tips the count over r
                for e in d[Network.EVENTS]:
                    count += e.rate()
                    if count >= r:
                        # Perform the event and stop iterations
                        e.perform(self)
                        return
            else:
                count += d[Network.TOTAL_RATE]

    def update_patch(self, patch_id, compartment_changes=None, attribute_changes=None):
        patch_data = self.node[patch_id]
        if compartment_changes:
            for (c,v) in compartment_changes.iteritems():
                patch_data[Network.COMPARTMENTS][c] += v
        if attribute_changes:
            for (a,v) in attribute_changes.iteritems():
                patch_data[Network.ATTRIBUTES][a] += v
        self.update_rates_at_patch(patch_id)

    def update_edge(self, u, v, attribute_changes):
        for att,value in attribute_changes.iteritems():
            self.edge[u][v][att] += value
        self.update_rates_at_patch(u)
        self.update_rates_at_patch(v)
