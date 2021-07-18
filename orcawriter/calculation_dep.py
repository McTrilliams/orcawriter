import pathlib
import csv
import copy
from typing import List

import pandas as pd
import numpy as np
from myclass.util import n_to_neg, neg_to_n
from myclass.util import normalize_vector, find_ax_lig
from scripts.xyzreader import read_xyz
from temp_dev.find_ax_lig_temp_dev import find_axial_lig_idx

def ox_state_tuples(low, high):
    """Create tuples of desired oxidation states"""
    return tuple(range(low, high+1))

d_counts = dict(Sc=3,
                Ti=4,
                V=5,
                Cr=6,
                Mn=7,
                Fe=8,
                Co=9,
                Ni=10,
                Cu=11,
                Zn=12)

transition_metals = dict(V=ox_state_tuples(2, 5),
                         Cr=ox_state_tuples(0, 6),
                         Mn=ox_state_tuples(0, 5),
                         Fe=ox_state_tuples(0, 3),
                         Co=ox_state_tuples(0, 3),
                         Ni=ox_state_tuples(0, 2),
                         Cu=ox_state_tuples(1, 2),
                         )


class Xyz(object):
    """
    Parent class for loading a .xyz file into a dataframe, centers the main
    atom, currently looks for first-row transition metals as default.
    Use .from_file to load from an xyz file, default initialization is by
    passing a dataframe.
    Other functionality:
    - Atom Count
    - Writing a new .xyz file
    """

    #
    def __init__(self,
                 df: pd.DataFrame,
                 main_atom: str=None
        ) -> None:
        """From a pd.DataFrame create an XYZ object that holds atoms and their x, y, z coordinates."""
        self.df = df # TODO does it make sense to write this with pandas?
        self._main_atom = main_atom

    @classmethod
    def fromfile(
            cls,
            path: str or pathlib.Path,
            main_atom: str=None
    ):
        """
        #Read an .xyz file and create an Xyz object.

        :param path:
        :return:
        """
        return cls(df=read_xyz(path), main_atom=main_atom)


    def atom_count(self) -> int:
        """Return the number of atoms in the XYZ object."""
        return self.df.shape[0]

    def find_atoms(self, element_symbols: List[str]) -> pd.DataFrame:
        """Return the atoms that match the given elemental symbol as a Series or DataFrame."""
        # TODO check if the symbol is in the periodic module.
        atom = self.df.loc[self.df['element'].isin(element_symbols)]
        return atom

    def center_atom(self, element_symbol: str) -> None or ValueError:
        """Center the xyz coord (self.df) on the given element_symbol (0,0,0), else return a ValueError.

        This will update the coordinate df

        #TODO could I move this into an outer function 'new_df = center_atom(Xyz)'"""
        atom_to_center = self.find_atoms([element_symbol])
        if atom_to_center.shape[0] > 1:
            return ValueError(f'More than one of {element_symbol} is present.')
        elif atom_to_center.shape[0] == 0:
            return ValueError(f'No {element_symbol} was found in the Xyz object.')
        for col in ['x', 'y', 'z']:
            self.df[col] = self.df[col].apply(lambda x: x - atom_to_center[col])

    def calc_dist_to_coord(self, x:int=0, y:int=0, z:int=0) -> None:
        """Calculate the distance to the origin or x y z coordinates input)"""
        self.df['dist2main'] = (((self.df["x"]-x)**2) + ((self.df["y"]-y)**2) + ((self.df["z"]-z)**2))**0.5

    #


    def write2csv(self, filepath):
        filepath = pathlib.Path(filepath)
        filepath = filepath.with_suffix('.xyz')
        dfcopy = self.df.copy()
        dfcopy = dfcopy.drop(['atom_id', 'dist2main'], axis=1)

        try:
            dfcopy.to_csv(
                            filepath,
                            sep="\t",
                            header=False,
                            encoding="ascii",
                            float_format="%10.6f",
                            quoting=csv.QUOTE_NONE,
                            index=False,
                            line_terminator="\n",
            )

            setattr(self, 'path', filepath)
            print(f"XYZ file: {filepath.stem} successfully written to file!")

            print(f"XYZ file: {filepath} successfully written to file!")

        except OSError:
            print(f"Unable to write the .xyz file.")

        shape = next(iter(list(dfcopy.shape)))

        with open(filepath, "r+") as file:
            content = file.read()
            if not content.startswith(str(shape)):
                file.seek(0, 0)
                file.write(f"{shape}\n{filepath.stem}\n" + content) # \n{filepath.stem}\n" +
                # content)

    def main_atom_symbol(self):
        return next(iter(list(self.main_atom['element'])))

    @property
    def main_atom(self):
        return self._main_atom

    @main_atom.setter
    def main_atom(self, new_atom, new_coord=None):
        idx = self.main_atom.index.to_list()
        idx = next(iter(idx))
        self.df.at[idx, 'element'] = new_atom
        self._main_atom = self.find_main_atom()

    def add_hydride(self, coord_num, hyd_dist=1.45):
        df = self.df
        df.sort_values(by=['dist2main'], inplace=True)
        df.reset_index(inplace=True, drop=True)
        # Find the donor atoms
        donor_vec = df.loc[1:(coord_num), "x":"z"]
        coord_vec_idx = find_axial_lig_idx(donor_vec, df)
        coord_vec = df.loc[coord_vec_idx, 'x':'z']
        coord_vec_np = coord_vec * -1
        coord_vector_np = normalize_vector(coord_vec_np, hyd_dist)

        # make the series to append to the xyz df
        h_vec = pd.Series(["H"], index=["element"])
        atom_id = pd.Series([df['atom_id'].max() + 1], index=['atom_id'])
        h_vec = pd.concat([atom_id, h_vec, coord_vector_np, pd.Series(np.nan, index=["dist2main"])])

        # Append the Hydride and sort based on distance.
        df = df.append(h_vec, ignore_index=True)
        return df

    def get_atom_coordinates(self, index: int) -> np.array:
        """Return the coordinates of the atom at the given index."""
        return np.array(self.df.loc[index, 'x':'z'])


class Parameters:
    '''Create a Parameter object from a dictionary'''
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


class Output:

    def __init__(self, out_path):
        self.out_path = out_path
        self.read_output_file()

    def read_output_file(self):
        pass


class Calculation:

    def __init__(self, path, **kwargs):
        self.path = path
        if path.suffix == '.xyz':
            self.xyz = Xyz.from_file(path)
        for key, value in kwargs.items():
            setattr(self, key, value)

    @staticmethod
    def par_from_filename(filename):
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

    @staticmethod
    def par_from_filename_old(filename):
        par = dict()
        parts = filename.split('_')
        par['ccdc_id'] = parts[0]
        parts.remove(par["ccdc_id"])
        poss_par = dict(
                        structure='',
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

    @classmethod
    def from_file(cls, path):
        par = cls.par_from_filename(path.stem)
        return cls(path, **par)

    @classmethod
    def from_file_old(cls, path):
        par = cls.par_from_filename_old(path.stem)
        return cls(path, **par)

    def load_stdpar(self, stdpar_path):
        self.stdpar = Parameters.from_file(stdpar_path)

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

    def gen_calc_name(self):
        ox_state = neg_to_n(self.ox_state)
        charge = neg_to_n(self.charge)
        multiplicity = neg_to_n(self.multiplicity)
        lig_charge = neg_to_n(self.lig_charge)
        metal_id = self.xyz.main_atom_symbol()

        name = f"{self.ccdc_id}_{self.structure}_{metal_id}_os{ox_state}_c" \
               f"{charge}_m{multiplicity}_lc{lig_charge}"
        # if 'SinglePoint' == self.calc_type:
        #     name = f'{name}_SP'

        setattr(self, 'name', name)

    def gen_xyz_name(self):
        setattr(self, 'xyzname', f'{self.name}_initxyz')

    def __str__(self):
        if not self.name:
            self.name = self.gen_calc_name()
        return str(self.name)

    def __repr__(self):
        return self.__str__()
    #
    # def __eq__(self, other):
    #     self_dict = copy.copy(self.__dict__)
    #     other_dict = copy.copy(other.__dict__)
    #     drop_ele = ['path', 'xyz', 'stdpar']
    #     for ele in drop_ele:
    #         del self_dict[ele]
    #         del other_dict[ele]
    #     return self_dict == other_dict

'''
    def par_from_filename(self):
        if 'SP' in part_list:
            par['calc_type'] = 'single_point'
        elif 'NF' in part_list:
            par['calc_type'] = 'NumFreq'
        else:
            par['calc_type'] = 'TightOpt'

        # Initialize other defaults from phase 1
        par['hydride_constrain'] = False
        return par
'''

class ReorgEnergyEstimate:
    def __init__(self, met_calc, hyd_calc, eMet=None, eHyd=None):
        self.met = met_calc
        self.hyd = hyd_calc
        self.eMet = eMet
        self.eHyd = eHyd



