

class Parameters:
    '''Create a container for a set of parameters'''
    #Todo Can I sub class a dictionary? Just use NamedTuple? Dataclass?
    # TODO Possibly use configParser? and incorporate as specialized methods in CalculationBase
    def __init__(self, par: dict):
        for key, value in par.items():
            setattr(self, key, value)

    @classmethod
    def from_file(cls, path):
        '''
        Load a parameter set from a .txt file with each line as 'param = value'
        '''
        file_dict = dict()
        with open(path, "r") as file:
            lines = file.readlines()
            lines = [x.strip("\n") for x in lines]
        for line in lines:
            line = line.split(" = ")
            file_dict[line[0]] = line[1]
        return cls(file_dict)


class OrcaParam(Parameters):
    """Parameter validation specific to ORCA."""

    def __init__(self,
                 specialpar: dict=None,
                 param_set: Parameters=None,
                 ):
        """Takes unique input params (specialpar) that override the param_set defaults.

        specialpar: dict of kwargs to override or set certain parameters
        param_set: default paramset
        """
        par = self.with_validation(specialpar, param_set)
        super().__init__(self, par)

    @staticmethod
    def with_validation(overrides: dict, defaults: Parameters):
        for p in overrides:
            pass
        return NotImplementedError

    def par_from_filename(filename):
        """This could be a user setting area? where they pass a schema?
        Could pass a basic f-string.

        """
        par = dict()
        parts = filename.split('_')
        par['ccdc_id'] = parts[0]
        parts.remove(par["ccdc_id"])
        poss_par = dict(structure='',
                        lig_charge='lc',
                        ox_state='os',
                        charge='c',
                        multiplicity='m',
                        )
        for key, value in poss_par.items():
            par[key] = next((part.strip(value) for part in parts
                            if part.startswith(value)), None)
            if par[key]:
                par[key] = n_to_neg(par[key])
        return par