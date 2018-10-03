import networkx
import numpy


class Network(networkx.Graph):
    """
    A MetapopPy network. Extends networkX graph, adding data for patch subpopulations and environmental attributes.
    """

    COMPARTMENTS = 'compartments'
    ATTRIBUTES = 'attributes'

    def __init__(self, compartments, patch_attributes, edge_attributes, template=None):
        """
        Create a network
        :param compartments: List of population compartments
        :param patch_attributes: List of patch attributes
        :param edge_attributes: List of edge attributes
        """
        self._compartments = compartments
        self._patch_attributes = patch_attributes
        self._edge_attributes = edge_attributes
        self._handler = None
        networkx.Graph.__init__(self)

        if template:
            self.add_nodes_from(template.nodes())
            self.add_edges_from(template.edges())

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

    def prepare(self, handler=None):
        """
        Prepare the network. Adds data to the node dictionary for subpopulation, environmental attributes and events
        and rates. All set to zero, will be seeded with a value once simulation runs.
        :return:
        """
        self._prepare_compartments()
        self._prepare_attributes()
        self._prepare_edges()
        self._handler = handler

    def _prepare_compartments(self):
        """
        Prepare the network compartments (assumes all patches have same compartments)
        :return:
        """
        for _, patch_data in self.nodes(data=True):
            patch_data[Network.COMPARTMENTS] = dict([(c, 0) for c in self._compartments])

    def _prepare_attributes(self):
        """
        Prepare the network attributes (assumes all patches have same attributes)
        :return:
        """
        for _, patch_data in self.nodes(data=True):
            patch_data[Network.ATTRIBUTES] = dict([(a, 0) for a in self._patch_attributes])

    def _prepare_edges(self):
        """
        Prepare edge attributes (assumes all edges have the same attributes)
        :return:
        """
        for u, v in self.edges():
            for a in self._edge_attributes:
                self.edge[u][v][a] = 0.0

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

    def update_patch(self, patch_id, compartment_changes=None, attribute_changes=None):
        """
        Update the given patch with the given changes. If a handler is attached, this will be called with the changes
        in order to propagate the updates.
        :param patch_id: ID of patch changed
        :param compartment_changes: dict of Key:compartment, Value: amount changed
        :param attribute_changes: dict of Key:attribute, Value: amount changed
        :return:
        """
        patch_data = self.node[patch_id]
        if compartment_changes:
            for comp, change in compartment_changes.iteritems():
                patch_data[Network.COMPARTMENTS][comp] += change
                assert patch_data[Network.COMPARTMENTS][comp] >= 0, "Compartment {0} cannot drop below zero".format(comp)
        if attribute_changes:
            for attr, change in attribute_changes.iteritems():
                patch_data[Network.ATTRIBUTES][attr] += change
        # Propagate the changes
        if self._handler:
            if not compartment_changes:
                compartment_changes = {}
            if not attribute_changes:
                attribute_changes = {}
            self._handler(patch_id, compartment_changes.keys(), attribute_changes.keys(), {})

    def update_edge(self, u,v, attribute_changes):
        edge = self.edge[u][v]
        for attr, change in attribute_changes.iteritems():
            edge[attr] += change
        if self._handler:
            self._handler(u, [], [], attribute_changes.keys())
            self._handler(v, [], [], attribute_changes.keys())


class TypedNetwork(Network):
    """
    A MetapopPy network where patches are assigned a "type", which can restrict which dynamics occurs there.
    """

    PATCH_TYPE = 'patch_type'

    def __init__(self, compartments, patch_attributes_by_type, edge_attributes):
        """
        Create a network
        :param compartments: List of population compartments
        :param patch_attributes_by_type: List of patch attributes, grouped by patch type
        :param edge_attributes: List of edge attributes
        """

        self._attribute_by_type = patch_attributes_by_type

        all_patch_attributes = []
        for a in patch_attributes_by_type.values():
            all_patch_attributes += a
        all_patch_attributes.append(TypedNetwork.PATCH_TYPE)
        Network.__init__(self, compartments, all_patch_attributes, edge_attributes)
        self._patch_types = {}

    def set_patch_type(self, patch_id, patch_type):
        self.node[patch_id][TypedNetwork.PATCH_TYPE] = patch_type

    def get_patches_by_type(self, patch_type, data=False):
        if patch_type not in self._patch_types:
            return []
        elif data:
            return [(n, self.node[n]) for n in self._patch_types[patch_type]]
        else:
            return self._patch_types[patch_type]

    def prepare(self, handler=None):
        # Ensure every patch has been given a type
        for patch_id, patch_data in self.nodes(data=True):
            assert TypedNetwork.PATCH_TYPE in patch_data, "Node {0} must be assigned a patch type".format(patch_id)
            # Create the shortcut list for this patch type if we haven't seen it yet
            if patch_data[TypedNetwork.PATCH_TYPE] not in self._patch_types:
                self._patch_types[patch_data[TypedNetwork.PATCH_TYPE]] = []
            # Add to the list
            self._patch_types[patch_data[TypedNetwork.PATCH_TYPE]].append(patch_id)
        Network.prepare(self, handler)

    def _prepare_attributes(self):
        # Prepare the network attributes
        for _, patch_data in self.nodes(data=True):
            p_type = patch_data[TypedNetwork.PATCH_TYPE]
            if p_type in self._attribute_by_type:
                attributes = self._attribute_by_type[p_type]
                patch_data[Network.ATTRIBUTES] = dict([(a, 0) for a in attributes])
            else:
                patch_data[Network.ATTRIBUTES] = {}