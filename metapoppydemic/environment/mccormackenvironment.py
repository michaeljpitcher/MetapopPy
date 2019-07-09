from metapoppy import Environment
import networkx


class McCormackEnvironment(Environment):

    BIRTH_RATE = 'birth_rate'
    BASE_DEATH_RATE = 'base_death_rate'
    POPULATION_DEATH_RATE = 'population_death_rate'
    INFECTION_LAMBDA = 'infection_lambda' # lambda
    CARRYING_CAPACITY = 'carrying_capacity' # K
    INFECTION_DEATH_RATE = 'infection_death_rate' # alpha
    RECOVERY_RATE = 'recovery_rate' # gamma

    def __init__(self, compartments):
        template = networkx.Graph()
        template.add_edges_from([(1, 2), (2, 3), (1, 3)])

        patch_attributes = [McCormackEnvironment.BIRTH_RATE,
                            McCormackEnvironment.BASE_DEATH_RATE, McCormackEnvironment.POPULATION_DEATH_RATE,
                            McCormackEnvironment.INFECTION_LAMBDA, McCormackEnvironment.CARRYING_CAPACITY,
                            McCormackEnvironment.INFECTION_DEATH_RATE, McCormackEnvironment.RECOVERY_RATE]
        edge_attributes = []

        Environment.__init__(self, compartments, patch_attributes, edge_attributes, template)