
"""
Tests for the Calculation object that is the basis for ORCA/DFT script generation.
"""
import pytest
from typing import Dict
from orcawriter.calculation import Calculation
from orcawriter.xyz import Xyz


@pytest.fixture
def xyz_data_file() -> str:
    r"""Path for the xyz test data file."""
    return r'src/tests/data/tst_xyz_data.xyz'


@pytest.fixture
def test_load_xyz(xyz_data_file: str) -> Xyz:
    r"""Load the xyz file data and return a populated Xyz object.
    :param xyz_data_file: Path for an ``.xyz file`` to use for testing.
    :rtype: orcawriter.Xyz
    """
    return Xyz.fromfile(xyz_data_file)


@pytest.fixture
def good_form_data() -> Dict:
    r"""An incomplete set of parameters for testing calculation methods.

    These simulate a full set of params from a web form.
    :rtype: Dict
    """
    params = {
        'calc_type': 'TightOpt',
        'functional' : 'B3LYP',
        'basis_set' : 'def2-TZVP',
        'charge': 1,
        'multiplicity': 6,
        'relativistic' : 'ZORA',
        'dispersion_correction' : 'D3BJ',
        'solvent_model' : 'CPCMC',
        'solvent' : 'toluene',
        'resolution_id' : True,
        'aux_basis_set' : 'def2-SVP',
        'xyz_name' : 'tst_xyz_data.xyz'
    }
    return params

@pytest.fixture
def test_build_calculation(
           good_form_data: Dict,
           test_load_xyz: Xyz,
    ):
    """Return a Calculation created with good form data and an Xyz object.

    :param good_form_data: Dict of complete form data for testing.
    :param test_load_xyz: Xyz object with data for testing.
    :return:
    """

    calc = Calculation(**good_form_data)     # Put in some parameters
    calc.xyz = test_load_xyz # Add xyz data
    calc.set_defaults(r'src/tests/data/tst_default_params.ini')
    return calc


def test_to_inp(test_build_calculation):
    with open('src/tests/data/inp_tst_template.txt', 'r') as file:
        inp_tst_template = file.read()
    assert test_build_calculation.to_inp() == inp_tst_template



@pytest.fixture
def bad_form_data() -> Dict:
    """An incomplete set of parameters for testing calculation methods.

    These simulate a full set of params from a web form."""
    params = {
        'calc_type': 'Opt',
        'functional' : 'B3LYP',
        'basis_set' : 'def2-TVP', #typo
        'charge': 1,
        'multiplicity': 6,
        'relativistic' : 'ZORA',
        'dispersion_correction' : 'D3BJ',
        'solvent_model' : 'CPCMC',
        'solvent' : '1,2-dichlorobenzene', # bad solvent choice
        'resolution_id' : True,
        'aux_basis_set' : 'def2-SVP',
        'xyz_name' : 'test_xyz_mol',
    }
    return params
