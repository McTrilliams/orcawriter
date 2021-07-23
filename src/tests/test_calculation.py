
"""
Tests for the Calculation object that is the basis for ORCA/DFT script generation.
"""

import pytest
from typing import Dict
from orcawriter.calculation import Calculation
from orcawriter.xyz import Xyz

@pytest.fixture
def good_form_data() -> Dict:
    """An incomplete set of parameters for testing calculation methods.

    These simulate a full set of params from a web form."""
    params = {
        'calc_type': 'Opt',
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
        'xyz_name' : 'test_xyz_mol',
    }
    return params

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

@pytest.fixture
def xyz_data_file() -> Xyz:
    return r'orcawriter/templates/tst_xyz_data.xyz'

# @pytest.fixture
# def additional_params():
#     return NotImplementedError

def test_build_calculation(
           form_data: Dict,
           xyz_data: Xyz,
           additional_params: Dict
    ):
    calc = Calculation()
    return calc

# def test_to_inp(calc=test_create_calculation):
#     with open('inp_tst_template.txt', 'r') as file:
#         inp_tst_template = file.read()
#     assert calc.to_inp() == inp_tst_template


