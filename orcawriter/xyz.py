# xyz.py
# Definitions of the Xyz object and it's children.

# If the implementation is easy to explain, it may be a good idea.
import pandas as pd
import pathlib
import csv
from scripts.xyzreader import read_xyz
from typing import List


class XyzBase(object):
    """
    Parent class for loading a .xyz file into a dataframe, centers the main
    atom, currently looks for first-row transition metals as default.
    Use .from_file to load from an xyz file, default initialization is by
    passing a dataframe.
    Other functionality:
    - Atom Count
    - Writing a new .xyz file
    """

    def __init__(self,
                 df: pd.DataFrame,
                 *args,
                 **kwargs,
        ) -> None:
        """From a pd.DataFrame create an XYZ object that holds atoms and their x, y, z coordinates."""
        # TODO add validation methods to the xyz class below
        self.df = df
        for kwarg in kwargs:
            setattr(self, kwarg, kwargs[kwarg])

    def __repr__(self):
        if hasattr(self, 'path'):
            return f'path = {self.path}\n{self.df.head(5)}'
        else:
            return f'{self.df.head(5)}'

    @classmethod
    def fromfile(
            cls,
            path: str or pathlib.Path=None,
            *args,
            **kwargs
    ):
        """
        #Read an .xyz file and create an Xyz object.

        :param path:
        :return:
        """
        df = read_xyz(path)
        return cls(df=df, path=path, *args, **kwargs)

    def atom_count(self) -> int:
        """Return the number of atoms in the XYZ object."""
        return self.df.shape[0]

    def find_atoms(self, element_symbols: List[str] or str) -> pd.DataFrame:
        """Return the atoms that match the given elemental symbol as a Series or DataFrame."""
        # TODO check if the symbol is in the periodic module.
        if isinstance(element_symbols, str):
            element_symbols = [element_symbols]
        atom = self.df.loc[self.df['element'].isin(element_symbols)]
        return atom

    def center_atom(self, element_symbol: str) -> None or ValueError:
        """Center the xyz coord (self.df) on the given element_symbol (0,0,0), else return a ValueError.

        This will update the coordinate df

        #TODO could I move this into an outer function 'new_df = center_atom(Xyz)'"""
        atom_to_center = self.find_atoms([element_symbol])
        if atom_to_center.shape[0] > 1:
            return ValueError(f'More than one of {element_symbol} atom is present.')
        elif atom_to_center.shape[0] == 0:
            return ValueError(f'No {element_symbol} atoms were found in the Xyz object.')
        for col in ['x', 'y', 'z']:
            self.df[col] = self.df[col].apply(lambda x: x - atom_to_center[col])

    def calc_dist_to_coord(self, x:int=0, y:int=0, z:int=0, column:str='dist2origin') -> None:
        """Calculate the distance to the origin or x y z coordinates input in a new column.

        x, y, z coordinates to calc distance to.
        :param column str name for new column. default='dist2origin'
        """
        self.df[column] = (((self.df["x"]-x)**2) + ((self.df["y"]-y)**2) + ((self.df["z"]-z)**2))**0.5

## Area in Work - Proceed with caution and compassion. ##

    def write2csv(self,
                  filepath: str or pathlib.Path,
                  setattr_path: bool=False
                  ) -> None:
        """Write the object to a .xyz file for transfer to a chemical program, namely ORCA."""
        if not isinstance(filepath, pathlib.Path):
            filepath = pathlib.Path(filepath)
        if not filepath.suffix == '.xyz':
            filepath = filepath.with_suffix('.xyz')
        dfcopy = self.df.copy()
        dropcols = [col for col in self.df.columns if col not in {'element', 'x', 'y', 'z'}]
        dfcopy = dfcopy.drop(dropcols, axis=1)
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
            if setattr_path:
                setattr(self, 'path', filepath)
            print(f"XYZ file: {filepath} successfully written to file!")
        except OSError:
            print(f"Unable to write {filepath})")
        shape = dfcopy.shape[0]

        with open(filepath, "w") as file:
            content = file.read()
            # Todo check if the description is there, track the description through the stack.
            if not content.startswith(str(shape)):
                file.seek(0, 0)
                file.write(f"{shape}\n{filepath.stem}\n" + content)


class Xyz(XyzBase):
    """Subclass for name change. Ain't no longer a basic B$%..."""
    pass


