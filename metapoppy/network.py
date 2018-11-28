import networkx


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

    def set_handler(self, handler):
        """
        Given a lambda function as an update handler, assign it to the network. Function will be called when an update
        is performed
        :param handler:
        :return:
        """
        self._handler = handler

    def reset(self):
        self._reset_patches()
        self._reset_edges()

    def _reset_patches(self):
        networkx.set_node_attributes(self, {n: {TypedNetwork.COMPARTMENTS: {c: 0 for c in self._compartments},
                                                TypedNetwork.ATTRIBUTES: {a: 0.0 for a in self._patch_attributes}}
                                            for n in self.nodes})

    def _reset_edges(self):
        networkx.set_edge_attributes(self, {(u,v): {a: 0.0 for a in self._edge_attributes} for (u,v) in self.edges})

    def get_compartment_value(self, patch_id, compartment):
        """
        Get function for finding a compartment value (or values) at a patch
        :param patch_id:
        :param compartment:
        :return:
        """
        if isinstance(compartment, list):
            data = self._node[patch_id][Network.COMPARTMENTS]
            return sum([data[c] for c in compartment])
        else:
            return self._node[patch_id][Network.COMPARTMENTS][compartment]

    def get_attribute_value(self, patch_id, attribute):
        """
        Get function for finding an environmental attribute value at a patch
        :param patch_id:
        :param attribute:
        :return:
        """
        if isinstance(attribute, list):
            data = self._node[patch_id][Network.ATTRIBUTES]
            return sum([data[c] for c in attribute])
        else:
            return self._node[patch_id][Network.ATTRIBUTES][attribute]

    def update_patch(self, patch_id, compartment_changes=None, attribute_changes=None):
        """
        Update the given patch with the given changes. If a handler is attached, this will be called with the changes
        in order to propagate the updates.
        :param patch_id: ID of patch changed
        :param compartment_changes: dict of Key:compartment, Value: amount changed
        :param attribute_changes: dict of Key:attribute, Value: amount changed
        :return:
        """
        patch_data = self.nodes[patch_id]
        if compartment_changes:
            for comp, change in compartment_changes.iteritems():
                patch_data[Network.COMPARTMENTS][comp] += change
                assert patch_data[Network.COMPARTMENTS][comp] >= 0, \
                    "Compartment {0} cannot drop below zero".format(comp)
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

    def update_edge(self, u, v, attribute_changes):
        """
        Update the attributes of an edge
        :param u: Patch 1
        :param v: Patch 2
        :param attribute_changes: dict of Key:attribute, Value: amount changed
        :return:
        """
        edge = self.get_edge_data(u, v)
        for attr, change in attribute_changes.iteritems():
            edge[attr] += change
        # If a handler exists, propagate the updates
        if self._handler:
            self._handler(u, [], [], attribute_changes.keys())
            self._handler(v, [], [], attribute_changes.keys())


class TypedNetwork(Network):
    """
    A MetapopPy network where patches are assigned a "type", which can be used to restrict which dynamics occurs there.
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
        """
        Assign the patch type to the patch
        :param patch_id:
        :param patch_type:
        :return:
        """
        self.nodes[patch_id][TypedNetwork.PATCH_TYPE] = patch_type
        if patch_type not in self._patch_types:
            self._patch_types[patch_type] = []
        self._patch_types[patch_type].append(patch_id)
        if patch_type not in self._attribute_by_type:
            self._attribute_by_type[patch_type] = []

    def get_patches_by_type(self, patch_type, data=False):
        """
        Return a list of patch IDs of the given type
        :param patch_type:
        :param data: If true, also returns the patch data
        :return:
        """
        if patch_type not in self._patch_types:
            return []
        elif data:
            return [(n, self.nodes[n]) for n in self._patch_types[patch_type]]
        else:
            return self._patch_types[patch_type]

    def _reset_patches(self):
        networkx.set_node_attributes(self, {n: {TypedNetwork.PATCH_TYPE: self.node[n][TypedNetwork.PATCH_TYPE],
                                                TypedNetwork.COMPARTMENTS: {c: 0 for c in self._compartments},
                                                TypedNetwork.ATTRIBUTES: {a: 0 for a in
                                                    self._attribute_by_type[self.node[n][TypedNetwork.PATCH_TYPE]]}}
                                            for n in self.nodes})
