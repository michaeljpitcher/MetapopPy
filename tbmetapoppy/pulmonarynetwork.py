from metapoppy.network import TypedNetwork
import numpy


class PulmonaryNetwork(TypedNetwork):

    TOPOLOGY = 'topology'
    SPACE_FILLING_TREE_2D = 'space_filling_tree_2d'

    BOUNDARY = 'boundary'
    LENGTH_DIVISOR = 'length_divisor'
    MINIMUM_AREA = 'minimum_area'

    ALVEOLAR_PATCH = 'alveolar_patch'
    LYMPH_PATCH = 'lymph_patch'

    VENTILATION = 'ventilation'
    PERFUSION = 'perfusion'
    OXYGEN_TENSION = 'oxygen_tension'
    DRAINAGE = 'drainage'

    PATCH_ATTRIBUTES = {ALVEOLAR_PATCH: [VENTILATION, PERFUSION, OXYGEN_TENSION, DRAINAGE]}

    # TODO - more edge attributes
    EDGE_ATTRIBUTES = [PERFUSION]

    def __init__(self, network_config, compartments):
        TypedNetwork.__init__(self, compartments, PulmonaryNetwork.PATCH_ATTRIBUTES, PulmonaryNetwork.EDGE_ATTRIBUTES)

        self._alveolar_positions = []
        if network_config[PulmonaryNetwork.TOPOLOGY] == PulmonaryNetwork.SPACE_FILLING_TREE_2D:
            self._build_2d_space_filling_tree(network_config)

    def build_bps_network(self):
        # TODO - build bps_network
        pass

    def _build_2d_space_filling_tree(self, network_config):

        boundary = network_config[PulmonaryNetwork.BOUNDARY]
        length_divisor = network_config[PulmonaryNetwork.LENGTH_DIVISOR]
        minimum_area = network_config[PulmonaryNetwork.MINIMUM_AREA]

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
            # Loop through all points on the boundary (except split point and points immediately closest to the split point
            # (i.e. first two and last on boundary)
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

        self.add_node(PulmonaryNetwork.LYMPH_PATCH)
        self.set_patch_type(PulmonaryNetwork.LYMPH_PATCH, PulmonaryNetwork.LYMPH_PATCH)

        index = 1
        for pos in alveolar_positions:
            self.add_node(index)
            self.set_patch_type(index, PulmonaryNetwork.ALVEOLAR_PATCH)
            self._alveolar_positions.append(pos)

            self.add_edge(PulmonaryNetwork.LYMPH_PATCH, index)
            index += 1
