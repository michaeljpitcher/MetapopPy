import json

# TODO - repetitions
class MetapoppyOutput(object):
    def __init__(self, data):
        self.parameters = data['parameters']
        self.timesteps = [float(k) for k in data['results'].keys()]
        self.timesteps.sort()
        self.data = [data['results'][str(t)] for t in self.timesteps]

    def all_data_by_compartments(self):
        comps = self.data[0].values()[0]['compartments'].keys()
        comps.sort()
        return self.timesteps, {c: self.data_by_compartment(c) for c in comps}

    def data_by_compartment(self, compartment):
        nodes = set([item for sublist in [d.keys() for d in self.data] for item in sublist])
        c_data = {k: [] for k in nodes}

        for t in self.data:
            for n in nodes:
                if n in t:
                    c_data[n].append(t[n]['compartments'][compartment])
                else:
                    c_data[n].append(0)
        return c_data


class MetapoppyResultSet(object):
    def __init__(self, json_filename):
        with open(json_filename) as data_file:
            json_file = json.load(data_file)
        data = json_file['results'].values()
        self.param_variations = len(data)
        self.repetitions = len(data[0])
        print '{0}: {1} parameter variations with {2} repetitions'.format(self.__class__.__name__,
                                                                          self.param_variations, self.repetitions)
        self.results = []

        if self.repetitions > 1:
            pass
        else:
            for d in [n[0] for n in data]:
                self.results.append(MetapoppyOutput(d))

    def all_data_by_compartments(self):
        return [(o.parameters, o.all_data_by_compartments()) for o in self.results]