from metapoppy.environment import TypedEnvironment
import numpy
import ConfigParser


class TBPulmonaryEnvironment(TypedEnvironment):
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
    CYTOKINE = 'cytokine'
    EDGE_ATTRIBUTES = [PERFUSION, CYTOKINE]

    # Configuration
    TOPOLOGY = 'topology'
    SINGLE_PATCH = 'single_patch'
    SPACE_FILLING_TREE_2D = 'space_filling_tree_2d'
    BOUNDARY = 'boundary'
    LENGTH_DIVISOR = 'length_divisor'
    MINIMUM_AREA = 'minimum_area'

    VENTILATION_SKEW = 'ventilation_skew'
    PERFUSION_SKEW = 'perfusion_skew'
    DRAINAGE_SKEW = 'drainage_skew'

    # Compartments
    BACTERIUM_EXTRACELLULAR_REPLICATING = 'b_er'
    BACTERIUM_EXTRACELLULAR_DORMANT = 'b_ed'
    EXTRACELLULAR_BACTERIA = [BACTERIUM_EXTRACELLULAR_REPLICATING, BACTERIUM_EXTRACELLULAR_DORMANT]
    BACTERIUM_INTRACELLULAR_DENDRITIC = 'b_id'
    BACTERIUM_INTRACELLULAR_MACROPHAGE = 'b_im'
    INTRACELLULAR_BACTERIA = [BACTERIUM_INTRACELLULAR_DENDRITIC, BACTERIUM_INTRACELLULAR_MACROPHAGE]
    BACTERIA = EXTRACELLULAR_BACTERIA + INTRACELLULAR_BACTERIA

    MACROPHAGE_RESTING = 'm_r'
    MACROPHAGE_INFECTED = 'm_i'
    MACROPHAGE_ACTIVATED = 'm_a'
    MACROPHAGES = [MACROPHAGE_RESTING, MACROPHAGE_INFECTED, MACROPHAGE_ACTIVATED]

    DENDRITIC_CELL_IMMATURE = 'd_i'
    DENDRITIC_CELL_MATURE = 'd_m'
    DENDRITIC_CELLS = [DENDRITIC_CELL_IMMATURE, DENDRITIC_CELL_MATURE]

    T_CELL_NAIVE = 't_n'
    T_CELL_ACTIVATED = 't_a'
    T_CELLS = [T_CELL_NAIVE, T_CELL_ACTIVATED]

    CASEUM = 'c'

    TB_COMPARTMENTS = BACTERIA + MACROPHAGES + DENDRITIC_CELLS + T_CELLS + [CASEUM]

    ACTIVATED_CELL = {MACROPHAGE_RESTING: MACROPHAGE_ACTIVATED, T_CELL_NAIVE: T_CELL_ACTIVATED}
    INFECTED_CELL = {MACROPHAGE_RESTING: MACROPHAGE_INFECTED, DENDRITIC_CELL_IMMATURE: DENDRITIC_CELL_MATURE}
    INTERNAL_BACTERIA_FOR_CELL = {MACROPHAGE_INFECTED: BACTERIUM_INTRACELLULAR_MACROPHAGE,
                                  DENDRITIC_CELL_MATURE: BACTERIUM_INTRACELLULAR_DENDRITIC}

    def __init__(self, network_config):
        """
        Create a pulmonary environment for studying TB.
        :param network_config:
        """
        TypedEnvironment.__init__(self, TBPulmonaryEnvironment.TB_COMPARTMENTS, TBPulmonaryEnvironment.PATCH_ATTRIBUTES,
                                  TBPulmonaryEnvironment.EDGE_ATTRIBUTES)

        self._alveolar_positions = {}
        self._pulmonary_att_seeding = {}
        self._topology = network_config[TBPulmonaryEnvironment.TOPOLOGY]
        if self._topology == TBPulmonaryEnvironment.SINGLE_PATCH:
            self._build_single_patch_network()
        elif self._topology == TBPulmonaryEnvironment.SPACE_FILLING_TREE_2D:
            self._build_2d_space_filling_tree(network_config)
            # Get horizontal positions and max/min values
            ys = [y for _, y in self._alveolar_positions.values()]
            self._y_max = max(ys)
            y_min = min(ys)
            self._y_range = self._y_max - y_min

        self._infected_patches = []

    def output_positions(self, filename):
        """
        Write the positions of all nodes to a config file (to avoid building a topology over and over)
        :param filename:
        :return:
        """
        cp = ConfigParser.ConfigParser()
        cp.add_section(TypedEnvironment.POSITION)

        minx = maxx = miny = maxy = None
        for n, v in self._alveolar_positions.iteritems():
            if not minx:
                minx = v[0]
                maxx = v[0]
                miny = v[1]
                maxy = v[1]
            else:
                minx = min(minx, v[0])
                maxx = max(maxx, v[0])
                miny = min(miny, v[1])
                maxy = max(maxy, v[1])

            cp.set(TypedEnvironment.POSITION, str(n), str(v)[1:-1])

        cp.set(TypedEnvironment.POSITION, TBPulmonaryEnvironment.LYMPH_PATCH,
               str(maxx + (maxx-minx)/2.0) + ', ' + str(miny + (maxy-miny)/2.0))

        with open(filename, 'w') as file_pos:
            cp.write(file_pos)

    def infected_patches(self):
        """
        Get all infected patches
        :return:
        """
        return self._infected_patches

    def _build_single_patch_network(self):
        """
        Build a topology where the lung is a single patch
        :return:
        """
        # Create lymph patch
        self.add_node(TBPulmonaryEnvironment.LYMPH_PATCH)
        self.set_patch_type(TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.LYMPH_PATCH)

        # Create the alveolar patches
        self.add_node(TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self.set_patch_type(TBPulmonaryEnvironment.ALVEOLAR_PATCH, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self._alveolar_positions[TBPulmonaryEnvironment.ALVEOLAR_PATCH] = (1, 1)
        # Add an edge from this alveolar patch to lymph patch
        self.add_edge(TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.ALVEOLAR_PATCH)

    def _build_2d_space_filling_tree(self, network_config):
        """
        Create a 2D space filling tree to form the topology. Process:
        - Takes a boundary value, a start point on the boundary, a length divisor and a minimum area.
        - From start point, find another point on the boundary such that line from start point evenly bisects area.
        - Place a new point on the line, of distance length divisor of whole length.
        - Now we have (evenly sized) areas, and a start point for both
        - Repeat, until the size of the areas drops below the minimum area
        :param network_config:
        :return:
        """
        boundary = network_config[TBPulmonaryEnvironment.BOUNDARY]
        if isinstance(boundary, str):
            # TODO - unnecessary removal of brackets - just don't include in the first place
            boundary = [[float(a) for a in n[1:-1].split(",")] for n in boundary.split(":")]
        length_divisor = float(network_config[TBPulmonaryEnvironment.LENGTH_DIVISOR])
        minimum_area = float(network_config[TBPulmonaryEnvironment.MINIMUM_AREA])

        # apex = max([b[1] for b in boundary])
        # base = min([b[1] for b in boundary])

        def polygon_area(points):
            """
            Calculate the area of a polygon using the shoelace formula
            :param points:
            :return:
            """
            x = [q[0] for q in points]  # x coordinates
            y = [q[1] for q in points]  # y coordinates
            return 0.5 * numpy.abs(numpy.dot(x, numpy.roll(y, 1)) - numpy.dot(y, numpy.roll(x, 1)))

        def _find_halfway(parent_polygon, parent_size, tipping_point):
            # TODO - method is not perfect and requires rounding of values
            middle_point = ((parent_polygon[tipping_point][0] + parent_polygon[tipping_point - 1][0]) / 2.0,
                            (parent_polygon[tipping_point][1] + parent_polygon[tipping_point - 1][1]) / 2.0)

            area_using_midpoint = round(polygon_area(parent_polygon[0: tipping_point] + [middle_point]), 10)

            if area_using_midpoint == round(parent_size / 2.0, 10):
                return middle_point
            else:
                new_polygon = parent_polygon[0: tipping_point] + [middle_point] + parent_polygon[tipping_point:]
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
        self.add_node(TBPulmonaryEnvironment.LYMPH_PATCH)
        self.set_patch_type(TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.LYMPH_PATCH)

        # Create the alveolar patches
        index = 1
        for pos in alveolar_positions:
            self.add_node(index)
            self.set_patch_type(index, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
            self._alveolar_positions[index] = pos
            # Add an edge from this alveolar patch to lymph patch
            self.add_edge(TBPulmonaryEnvironment.LYMPH_PATCH, index)
            index += 1

    def calculate_pulmonary_attribute_values(self, params):
        """
        Calculate the initial values for pulmonary attributes - to be set when a patch becomes active
        :param params:
        :return:
        """
        # Only 1 alveolar patch so just set everything to 1
        if len(self._alveolar_positions) == 1:
            vent = params[TBPulmonaryEnvironment.VENTILATION]
            perf = params[TBPulmonaryEnvironment.PERFUSION]
            o2 = float(vent)/perf
            drain = params[TBPulmonaryEnvironment.DRAINAGE]
            seeding = {TBPulmonaryEnvironment.ALVEOLAR_PATCH: {TypedEnvironment.ATTRIBUTES:
                                                               {TBPulmonaryEnvironment.VENTILATION: vent,
                                                                TBPulmonaryEnvironment.PERFUSION: perf,
                                                                TBPulmonaryEnvironment.OXYGEN_TENSION: o2,
                                                                TBPulmonaryEnvironment.DRAINAGE: drain}}}
            self._pulmonary_att_seeding = seeding
            return seeding

        ventilation_skew = params[TBPulmonaryEnvironment.VENTILATION_SKEW]
        perfusion_skew = params[TBPulmonaryEnvironment.PERFUSION_SKEW]
        drainage_skew = params[TBPulmonaryEnvironment.DRAINAGE_SKEW]

        temp_seeding = {}

        total_v = 0
        total_q = 0

        # For ventilation and perfusion, value at apex will be 1, at base will be equal to skew. These are normalised
        # later to bring sum of all values for whole lung == 1
        # For drainage, value in middle of lung will equal 1, values for base and apex are derived from this and the
        # skew, and used to calculate values in between

        mins = {TBPulmonaryEnvironment.VENTILATION: 1.0, TBPulmonaryEnvironment.PERFUSION: 1.0,
                TBPulmonaryEnvironment.DRAINAGE: 1.0 / (1 + (drainage_skew - 1) / 2.0)}
        maxs = {TBPulmonaryEnvironment.VENTILATION: ventilation_skew, TBPulmonaryEnvironment.PERFUSION: perfusion_skew,
                TBPulmonaryEnvironment.DRAINAGE: mins[TBPulmonaryEnvironment.DRAINAGE] * drainage_skew}
        diffs = {k: maxs[k] - mins[k] for k in mins.keys()}

        # Calculate values - 1 + (y_max - y) * (att_max - att_min)/(y_max - y_min)
        for index, (_, y_pos) in self._alveolar_positions.iteritems():
            temp_seeding[index] = {}
            vent_value = mins[TBPulmonaryEnvironment.VENTILATION] + (self._y_max - y_pos) * (
                    diffs[TBPulmonaryEnvironment.VENTILATION] / self._y_range)
            perf_value = mins[TBPulmonaryEnvironment.PERFUSION] + (self._y_max - y_pos) * (
                    diffs[TBPulmonaryEnvironment.PERFUSION] / self._y_range)
            drain_value = mins[TBPulmonaryEnvironment.DRAINAGE] + (self._y_max - y_pos) * (
                    diffs[TBPulmonaryEnvironment.DRAINAGE] / self._y_range)
            total_v += vent_value
            total_q += perf_value
            # Temp values - will be normalised below
            temp_seeding[index] = {TBPulmonaryEnvironment.VENTILATION: vent_value,
                                   TBPulmonaryEnvironment.PERFUSION: perf_value,
                                   TBPulmonaryEnvironment.DRAINAGE: drain_value}

        self._pulmonary_att_seeding = {}
        # Normalise and assign
        for patch_id, values in temp_seeding.iteritems():
            self._pulmonary_att_seeding[patch_id] = {TypedEnvironment.ATTRIBUTES: {}}
            v = values[TBPulmonaryEnvironment.VENTILATION] / total_v
            q = values[TBPulmonaryEnvironment.PERFUSION] / total_q
            o2 = v / q
            self._pulmonary_att_seeding[patch_id][TypedEnvironment.ATTRIBUTES] = \
                {TBPulmonaryEnvironment.VENTILATION: v,
                 TBPulmonaryEnvironment.PERFUSION: q,
                 TBPulmonaryEnvironment.OXYGEN_TENSION: o2,
                 TBPulmonaryEnvironment.DRAINAGE: values[TBPulmonaryEnvironment.DRAINAGE]}
        return self._pulmonary_att_seeding

    def update_patch(self, patch_id, compartment_changes=None, attribute_changes=None):
        """
        Update a patch on the network. Also adds patch to list of infected patches, and updates the edge if the
        perfusion value of a lung patch has been amended
        :param patch_id: ID of patch
        :param compartment_changes: compartments changed
        :param attribute_changes: attributes changed
        :return:
        """
        if (patch_id not in self._infected_patches and
           self._node[patch_id][TypedEnvironment.PATCH_TYPE] == TBPulmonaryEnvironment.ALVEOLAR_PATCH and
           compartment_changes and any(x in TBPulmonaryEnvironment.BACTERIA for x in compartment_changes)):
            # TODO - awkward way of doing this
            self._infected_patches.append(patch_id)
        TypedEnvironment.update_patch(self, patch_id, compartment_changes, attribute_changes)
        if self._node[patch_id][TypedEnvironment.PATCH_TYPE] == TBPulmonaryEnvironment.ALVEOLAR_PATCH:
            if compartment_changes and TBPulmonaryEnvironment.MACROPHAGE_INFECTED in compartment_changes:
                self.update_edge(patch_id, TBPulmonaryEnvironment.LYMPH_PATCH,
                        {TBPulmonaryEnvironment.CYTOKINE: compartment_changes[TBPulmonaryEnvironment.MACROPHAGE_INFECTED]})
            if attribute_changes and TBPulmonaryEnvironment.PERFUSION in attribute_changes:
                val = attribute_changes[TBPulmonaryEnvironment.PERFUSION]
                self.update_edge(patch_id, TBPulmonaryEnvironment.LYMPH_PATCH, {TBPulmonaryEnvironment.PERFUSION: val})

    def reset(self):
        """
        Reset the environment. Also clears the infected patches list.
        :return:
        """
        self._infected_patches = []
        TypedEnvironment.reset(self)
