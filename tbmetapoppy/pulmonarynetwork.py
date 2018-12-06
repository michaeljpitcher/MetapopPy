from metapoppy.network import TypedNetwork
import numpy


class PulmonaryNetwork(TypedNetwork):
    # Patch types
    ALVEOLAR_PATCH = 'alveolar_patch'
    LYMPH_PATCH = 'lymph_patch'
    # Patch attributes
    VENTILATION = 'ventilation'
    PERFUSION = 'perfusion'
    OXYGEN_TENSION = 'oxygen_tension'
    DRAINAGE = 'drainage'
    PATCH_ATTRIBUTES = {ALVEOLAR_PATCH: [VENTILATION, PERFUSION, OXYGEN_TENSION, DRAINAGE]}

    # Edge attributes
    EDGE_ATTRIBUTES = [PERFUSION]

    # Configuration
    TOPOLOGY = 'topology'
    SINGLE_PATCH = 'single_patch'
    SPACE_FILLING_TREE_2D = 'space_filling_tree_2d'
    BOUNDARY = 'boundary'
    LENGTH_DIVISOR = 'length_divisor'
    MINIMUM_AREA = 'minimum_area'

    def __init__(self, network_config, compartments):
        TypedNetwork.__init__(self, compartments, PulmonaryNetwork.PATCH_ATTRIBUTES, PulmonaryNetwork.EDGE_ATTRIBUTES)

        self._alveolar_positions = {}
        if network_config[PulmonaryNetwork.TOPOLOGY] == PulmonaryNetwork.SINGLE_PATCH:
            self._build_single_patch_network()
        elif network_config[PulmonaryNetwork.TOPOLOGY] == PulmonaryNetwork.SPACE_FILLING_TREE_2D:
            self._build_2d_space_filling_tree(network_config)
            # Get horizontal positions and max/min values
            ys = [y for _, y in self._alveolar_positions.values()]
            self._y_max = max(ys)
            y_min = min(ys)
            self._y_range = self._y_max - y_min

    def _build_single_patch_network(self):
        # Create lymph patch
        self.add_node(PulmonaryNetwork.LYMPH_PATCH)
        self.set_patch_type(PulmonaryNetwork.LYMPH_PATCH, PulmonaryNetwork.LYMPH_PATCH)

        # Create the alveolar patches
        self.add_node(PulmonaryNetwork.ALVEOLAR_PATCH)
        self.set_patch_type(PulmonaryNetwork.ALVEOLAR_PATCH, PulmonaryNetwork.ALVEOLAR_PATCH)
        self._alveolar_positions[PulmonaryNetwork.ALVEOLAR_PATCH] = (1, 1)
        # Add an edge from this alveolar patch to lymph patch
        self.add_edge(PulmonaryNetwork.LYMPH_PATCH, PulmonaryNetwork.ALVEOLAR_PATCH)

    def _build_2d_space_filling_tree(self, network_config):
        boundary = network_config[PulmonaryNetwork.BOUNDARY]
        if isinstance(boundary, str):
            boundary = [[float(a) for a in n[1:-1].split(",")] for n in boundary.split(":")]
        length_divisor = float(network_config[PulmonaryNetwork.LENGTH_DIVISOR])
        minimum_area = float(network_config[PulmonaryNetwork.MINIMUM_AREA])

        # apex = max([b[1] for b in boundary])
        # base = min([b[1] for b in boundary])

        def polygon_area(points):
            """
            Calculate the area of a polygon using the shoelace formula
            :param points:
            :return:
            """
            x = [a for (a, b) in points]  # x coordinates
            y = [b for (a, b) in points]  # y coordinates
            return 0.5 * numpy.abs(numpy.dot(x, numpy.roll(y, 1)) - numpy.dot(y, numpy.roll(x, 1)))

        def _find_halfway(parent_polygon, parent_size, tipping_point):
            # TODO - method is not perfect and requires rounding of values
            mid_point = ((parent_polygon[tipping_point][0] + parent_polygon[tipping_point - 1][0]) / 2.0,
                         (parent_polygon[tipping_point][1] + parent_polygon[tipping_point - 1][1]) / 2.0)

            area_using_midpoint = round(polygon_area(parent_polygon[0: tipping_point] + [mid_point]), 10)

            if area_using_midpoint == round(parent_size / 2.0, 10):
                return mid_point
            else:
                new_polygon = parent_polygon[0: tipping_point] + [mid_point] + parent_polygon[tipping_point:]
                if area_using_midpoint < parent_size / 2.0:
                    tipping_point += 1
                return _find_halfway(new_polygon, parent_size, tipping_point)

        max_id = 0
        bronchial_positions = {str(max_id): boundary[0]}
        alveolar_positions = []
        bronchial_edges = []

        # Initialise the data set - start point, the parent boundary, the parent size
        data = [(max_id, boundary, polygon_area(boundary))]

        while data:
            split_patch_id, parent_boundary, parent_area = data.pop(0)
            # Loop through all points on the boundary (except split point and points immediately closest to the split
            # point (i.e. first two and last on boundary)
            new_split_point = None
            child_boundaries = None
            i = 2
            while not new_split_point:
                # Need to find a point on the boundary whereby the polygon from split point up to this point, following
                # all points on the parent boundary is of half the size of the parent boundary.
                # Create a boundary - from split point up to this point
                child_boundaries = [parent_boundary[0:i + 1], parent_boundary[i:] + [parent_boundary[0]]]
                # Calculate area
                new_area = polygon_area(child_boundaries[0])
                # New area is exactly half the parent area
                if new_area == parent_area / 2.0:
                    new_split_point = (parent_boundary[0][0] + (parent_boundary[i][0] - parent_boundary[0][0]) /
                                       length_divisor,
                                       parent_boundary[0][1] + (parent_boundary[i][1] - parent_boundary[0][1]) /
                                       length_divisor)
                    child_boundaries[0] = child_boundaries[0] + [new_split_point]
                    child_boundaries[1] = [new_split_point] + child_boundaries[1]
                # New area exceeds half the parent area - mid point lies between this point and the previous point
                elif new_area > parent_area / 2.0:
                    mid_point = _find_halfway(parent_boundary, parent_area, i)
                    new_split_point = (parent_boundary[0][0] + (mid_point[0] - parent_boundary[0][0]) / length_divisor,
                                       parent_boundary[0][1] + (mid_point[1] - parent_boundary[0][1]) / length_divisor)
                    child_boundaries[0] = [new_split_point] + parent_boundary[0:i] + [mid_point]
                    child_boundaries[1] = [new_split_point] + [mid_point] + parent_boundary[i:] + [parent_boundary[0]]
                # If new area is less than parent area / 2.0, then continue to the next point
                i += 1
            new_area = polygon_area(child_boundaries[0])
            max_id += 1
            bronchial_edges.append((split_patch_id, max_id))
            if new_area >= minimum_area:
                bronchial_positions[str(max_id)] = new_split_point
                data.append((max_id, child_boundaries[0], new_area))
                data.append((max_id, child_boundaries[1], new_area))
            else:
                alveolar_positions.append(new_split_point)

        # Create lymph patch
        self.add_node(PulmonaryNetwork.LYMPH_PATCH)
        self.set_patch_type(PulmonaryNetwork.LYMPH_PATCH, PulmonaryNetwork.LYMPH_PATCH)

        # Create the alveolar patches
        index = 1
        for pos in alveolar_positions:
            self.add_node(index)
            self.set_patch_type(index, PulmonaryNetwork.ALVEOLAR_PATCH)
            self._alveolar_positions[index] = pos
            # Add an edge from this alveolar patch to lymph patch
            self.add_edge(PulmonaryNetwork.LYMPH_PATCH, index)
            index += 1

    def get_pulmonary_att_seeding(self, ventilation_skew, perfusion_skew, drainage_skew):
        """
        Set the initial values for pulmonary attributes
        :param ventilation_skew:
        :param perfusion_skew:
        :param drainage_skew:
        :return:
        """
        # Only 1 alveolar patch so just set everything to 1
        if len(self._alveolar_positions) == 1:
            seeding = {PulmonaryNetwork.ALVEOLAR_PATCH: {TypedNetwork.ATTRIBUTES: {PulmonaryNetwork.VENTILATION: 1,
                                                                                   PulmonaryNetwork.PERFUSION: 1,
                                                                                   PulmonaryNetwork.OXYGEN_TENSION: 1,
                                                                                   PulmonaryNetwork.DRAINAGE: 1}}}
            return seeding

        temp_seeding = {}

        total_v = 0
        total_q = 0

        # For ventilation and perfusion, value at apex will be 1, at base will be equal to skew. These are normalised
        # later to bring sum of all values for whole lung == 1
        # For drainage, value in middle of lung will equal 1, values for base and apex are derived from this and the
        # skew, and used to calculate values in between

        mins = {PulmonaryNetwork.VENTILATION: 1.0, PulmonaryNetwork.PERFUSION: 1.0,
                PulmonaryNetwork.DRAINAGE: 1.0 / (1 + (drainage_skew-1)/2.0)}
        maxs = {PulmonaryNetwork.VENTILATION: ventilation_skew, PulmonaryNetwork.PERFUSION: perfusion_skew,
                PulmonaryNetwork.DRAINAGE: mins[PulmonaryNetwork.DRAINAGE] * drainage_skew}
        diffs = {k: maxs[k] - mins[k] for k in mins.keys()}

        # Calculate values - 1 + (y_max - y) * (att_max - att_min)/(y_max - y_min)
        for index, (_,y_pos) in self._alveolar_positions.iteritems():
            temp_seeding[index] = {}
            vent_value = mins[PulmonaryNetwork.VENTILATION] + (self._y_max - y_pos) * (
                        diffs[PulmonaryNetwork.VENTILATION] / self._y_range)
            perf_value = mins[PulmonaryNetwork.PERFUSION] + (self._y_max - y_pos) * (
                        diffs[PulmonaryNetwork.PERFUSION] / self._y_range)
            drain_value = mins[PulmonaryNetwork.DRAINAGE] + (self._y_max - y_pos) * (
                        diffs[PulmonaryNetwork.DRAINAGE] / self._y_range)
            total_v += vent_value
            total_q += perf_value
            # Temp values - will be normalised below
            temp_seeding[index] = {PulmonaryNetwork.VENTILATION: vent_value,
                                   PulmonaryNetwork.PERFUSION: perf_value,
                                   PulmonaryNetwork.DRAINAGE: drain_value}

        seeding = {}
        # Normalise and assign
        for patch_id, values in temp_seeding.iteritems():
            seeding[patch_id] = {TypedNetwork.ATTRIBUTES: {}}
            v = values[PulmonaryNetwork.VENTILATION] / total_v
            q = values[PulmonaryNetwork.PERFUSION] / total_q
            o2 = v / q
            seeding[patch_id][TypedNetwork.ATTRIBUTES] = {PulmonaryNetwork.VENTILATION: v,
                                                          PulmonaryNetwork.PERFUSION: q,
                                                          PulmonaryNetwork.OXYGEN_TENSION: o2,
                                                          PulmonaryNetwork.DRAINAGE: values[PulmonaryNetwork.DRAINAGE]}
        return seeding


# TODO change in patch perfusion needs to feed through to edge
