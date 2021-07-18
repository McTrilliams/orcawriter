import pathlib
from typing import Dict, List
from jinja2 import Template

INP_TEMPLATE = r'orcawriter/orca_input_template.txt'
SH_TEMPLATE = r'my_class/orca_sh_template.txt'


class CalculationBase:
    """Base Class for an ORCA DFT calculation."""

    def __init__(self, args, **kwargs):
        """Create an empty object, unless kwargs are passed."""
        if kwargs:
            for key, value in kwargs.items():
                setattr(self, key, value)

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        """Return str with the class name, and wanted attrs."""
        repr_str = f'{self.__class__.__name__}: '
        wanted_attrs = ['id', 'name']
        for attr in wanted_attrs:
            repr_str += f'{attr}={getattr(self, attr, __default=None)}'
        return repr_str


class Calculation(CalculationBase):

    @classmethod
    def from_webapp(cls, param_dct: Dict, files: List):
        """
        Plan return a CalculationBase from the flask app.
        param_dct will come from the form data
        files should be a list
        :param param_dct:
        :param files: List[str or file objects]
        :return:
        """

        par = dict()
        return cls(**par)

    def write_inp(self, dest_dir):
        """
            Writes a .inp file for an ORCA calculation for LongLeaf cluster @UNC-CH
            Uses the parameters from the two dictionaries passed
            std_param is generated from the input standard parameter file
            inp_param is generated from user input
            :param param: Dictionary passed
            :return:
            """
        dest_dir = pathlib.Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)

        with open(
                f"{dest_dir}\\{self.name}.inp", "w+", newline="\n"
        ) as inp:

            inp.write(
                    f"! {self.stdpar.spin_restriction} {self.stdpar.functional} {self.stdpar.basis_set} "
            )
            if "SinglePoint" != self.stdpar.calc_type:
                inp.write(f"{self.stdpar.calc_type} ")
            # if "NumFreq" == inp["calc_type"]:
            #     inp.write("MOREAD ")
            inp.write(
                    f"{self.stdpar.relativistic} {self.stdpar.grid1} {self.stdpar.final_grid} \n"
                    f"! {self.stdpar.scf_level} {self.stdpar.dispersion} "
                    f"{self.stdpar.solvent_model}("
                    f"{self.stdpar.solvent}) \n\n"
                    f"%pal nprocs {self.stdpar.num_procs} end\n\n"
                    f"%rel OneCenter {self.stdpar.one_center_value} \nend\n\n"
                    f"%method\n\tRI {self.stdpar.resolution_id} \nend\n\n"
            )
            # if 'NumFreq' == inp_param['calc_type']:
            #     inp.write(f"%moinp \"{inp_param['gbw_filename']}\"\n\n")
            inp.write(f"%basis\n\tAux \"{self.stdpar.aux_basis_set}\"\nend\n\n")
            # if param['constrain']:
            #     inpfile.write(
            #         f"%geom Constraints\n\t{{{param['constraint_type']} {param[
            #         'atom_metal_index']} "
            #         f"{param['atom_hydride_index']} "
            #         f"{param['constraint_value']} C}}\n\tend"
            #     )
            #     if True == param['invert_constraints']:
            #         inpfile.write(f"invertConstraints {param['invert_constraints']}")
            #     inpfile.write(f"\nend\n")
            if True == self.stdpar.write_MOs:
                inp.write(f"%output\nPrint [ P_Basis ] 2\nPrint [ P_MOs ] 1\nend\n")

            inp.write(
                    f"%scf MaxIter {self.stdpar.max_iter} \n\tshift shift 0.3 erroff 0.0 "
                    f"end\nend\n\n"
                    f"* xyzfile {self.charge} {self.multiplicity} "
                    f"{self.xyzname}.xyz\n\n"
            )

    def write_shell(self, dest_dir):
        """
        Writes a .sh file for an ORCA calculation for LongLeaf cluster
        Uses the parameters from the dictionaries passed to it
        :param param: Dictionary passed
        :return:
        """
        dest_dir = pathlib.Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        with open(
                f"{dest_dir}\\{self.name}.sh", "w+", newline="\n"
        ) as shFile:
            shFile.write(
                    f"#!/bin/bash\n#SBATCH -N {self.stdpar.num_nodes}"
                    f"\n#SBATCH -n {self.stdpar.num_procs}"
                    f"\n#SBATCH -t {self.stdpar.wall_time}"
                    f"\n#SBATCH -p {self.stdpar.queue_name}"
                    f"\n\nexport OMPI_MCA_btl=vader,self,tcp\nexport"
                    f"PATH={self.stdpar.mpi_path}"
                    f"\nexport {self.stdpar.mpi_path_lib}"
                    f":$LD_LIBRARY_PATH"
                    f"\n\nmkdir {self.stdpar.scratch_path}/$SLURM_JOBID"
                    f"\nworkdir={self.stdpar.scratch_path}/$SLURM_JOBID"
                    f"\nexport ORC_SCRDIR=$workdir"
                    f"\ntdir=$ORCA_SCRDIR"
                    f"\n\nexport ORC_SCRDIR=$workdir"
                    f"\n\ncp {self.xyzname}.xyz $workdir"
            )
            if "NumFreq" == self.stdpar.calc_type:
                shFile.write(f"\ncp {self.gbw_name}.gbw")
            shFile.write(
                    f"\ncp {self.name}.inp $workdir"
                    f"\n\ncd $workdir\npwd"
                    f"\n\ncat $[SLURM_JOB_NODELIST] > $workdir/$job.nodes"
                    f"\n\nls -ltr"
                    f"\n\n{self.stdpar.orca_path} {self.name}.inp >& "
                    f"{self.name}.out\n\n\n"
            )

    def to_inp_dict(self) -> Dict:
        """Create a dict of params for the inp file writer."""
        inp_param_dct = dict()
        return inp_param_dct

    def to_sh_dict(self) -> Dict:
        """Create a dict of params for the inp file writer."""
        sh_param_dct = dict()
        return sh_param_dct

    def to_inp(self, file: str, template=INP_TEMPLATE):
        """Replaces 'write_inp' and creates the .inp text file for an ORCA calculation."""
        #creat dict
        inp_params = self.to_inp_dict()
        with open(template, 'r') as f:
            template_str = f.read()
        tm = Template(template_str)
        output_str = tm.render(**inp_params)
        with open(file, 'w+') as f:
            f.write(output_str)

    def to_sh(self, file: str, template=SH_TEMPLATE):
        """Replaces 'write_inp' and creates the .sh text file for an ORCA calculation."""
        #creat dict
        sh_params = self.to_sh_dict()
        with open(template, 'r') as f:
            template_str = f.read()
        tm = Template(template_str)
        output_str = tm.render(**sh_params)
        with open(file, 'w+') as f:
            f.write(output_str)

