import networkx
import numpy


class Network(networkx.Graph):
    """
    A MetapopPy network. Extends networkX graph, adding data for patch subpopulations and environmental attributes.
    """

    INITIAL = 'initial'
    COMPARTMENTS = 'compartments'
    ATTRIBUTES = 'attributes'
    EVENTS = 'events'
    TOTAL_RATE = 'total_rate'

    def __init__(self, compartments, patch_attributes, edge_attributes):
        """
        Create a network
        :param compartments: List of population compartments
        :param patch_attributes: List of patch attributes
        :param edge_attributes: List of edge attributes
        """
        self._compartments = compartments
        self._patch_attributes = patch_attributes
        self._edge_attributes = edge_attributes
        networkx.Graph.__init__(self)

    def compartments(self):
        """
        Get function for compartments
        :return:
        """
        return self._compartments

    def patch_attributes(self):
        """
        Get function for patch attributes
        :return:
        """
        return self._patch_attributes

    def edge_attributes(self):
        """
        Get function for edge attributes
        :return:
        """
        return self._edge_attributes

    def prepare(self):
        """
        Prepare the network. Adds data to the node dictionary, for subpopulation, environmental attributes and events
        and rates. All set to zero, will be seeded with a value once simulation runs.
        :return:
        """
        # Prepare the network
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
        """
        Get function for finding a compartment value at a patch
        :param patch_id:
        :param compartment:
        :return:
        """
        return self.node[patch_id][Network.COMPARTMENTS][compartment]

    def get_attribute_value(self, patch_id, attribute):
        """
        Get function for finding an environmental attribute value at a patch
        :param patch_id:
        :param attribute:
        :return:
        """
        return self.node[patch_id][Network.ATTRIBUTES][attribute]

    def attach_event(self, event, patch_id):
        """
        Attach the given event to the patch
        :param event:
        :param patch_id:
        :return:
        """
        event.set_patch_id(patch_id)
        self.node[patch_id][Network.EVENTS].append(event)

    def update_rates_at_patch(self, patch_id):
        """
        Loop through all events at a patch and trigger a rate update
        :param patch_id:
        :return:
        """
        # TODO - dependencies to improve efficiency here
        patch_data = self.node[patch_id]
        for event in patch_data[Network.EVENTS]:
            rate_before = event.rate()
            event.update_rate(patch_data, [d for _,_,d in self.edges([patch_id], data=True)])
            rate_difference = event.rate() - rate_before
            patch_data[Network.TOTAL_RATE] += rate_difference

    def update_patch(self, patch_id, compartment_changes=None, attribute_changes=None):
        """
        Update the values at a patch. Also triggers updates of the rates of events at the patch.
        :param patch_id: Patch ID in the network
        :param compartment_changes: Changes to compartments
        :param attribute_changes: Changes to environmental attributes.
        :return:
        """
        patch_data = self.node[patch_id]
        if compartment_changes:
            for (c,v) in compartment_changes.iteritems():
                patch_data[Network.COMPARTMENTS][c] += v
        if attribute_changes:
            for (a,v) in attribute_changes.iteritems():
                patch_data[Network.ATTRIBUTES][a] += v
        self.update_rates_at_patch(patch_id)

    def update_edge(self, u, v, attribute_changes):
        """
        Update the environmental attributes of the edge. Also triggers rate update of the events at the connected
        patches.
        :param u: Patch in network
        :param v: Patch in network
        :param attribute_changes: Changes to attributes
        :return:
        """
        for att,value in attribute_changes.iteritems():
            self.edge[u][v][att] += value
        self.update_rates_at_patch(u)
        self.update_rates_at_patch(v)
