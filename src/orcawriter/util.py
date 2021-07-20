import numpy as np
import math, cmath

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

d_count_mults = dict(
    d0=(1),
    d1=(2),
    d2=(1, 3),
    d3=(2, 4),
    d4=(1, 3, 5),
    d5=(2, 4, 6),
    d6=(1, 3, 5),
    d7=(2, 4),
    d8=(1, 3),
    d9=(2),
    d10=(1),
)

def neg_to_n(value: int) -> str or None:
    """Returns a str from int, and converts negative '-' -> 'n'. i.e. '-6' becomes 'n6'.

    Examples:
        +2 -> LC2
        0  -> LC0
        -7 -> LCn7

    If an integer is not entered a string of the input will be returned.
    None returns None.
    #TODO check on how to properly document this.
    :param value:
    :return: str
    """
    if value is None:
        return None
    elif not isinstance(value, int) or value >= 0:
        return str(value)
    else:
        return f"n{abs(int(value))}"


def n_to_neg(string: str) -> int or ValueError:
    """Returns an int from a str, and converts 'n' to negative 'n' -> '-'. i.e. 'n6' becomes '-6'.

    Examples:
        LC2 -> 2
        LC0  -> 0
        LCn7 -> -7

    A ValueError will be returned if the given string cannot be converted to int
    after stripping a leading 'n'.
    :param string: str
    :return: int
    # TODO Improve the quality of the logic here.
    """
    if string is None:
        return None
    elif not isinstance(string, str):
        return string
    else:
        try:
            strip = False
            if string.startswith('n'):
                string = string.strip('n')
                strip = True
                return int(string) * -1
            else:
                return int(string)
        except ValueError:
            if strip:
                return f'n{string}'
            else:
                return string




def normalize_vector(
        vector: np.array,
        length: float
):
    # TODO TEST
    """Normalize a vector to given length"""
    return vector * (length / np.linalg.norm(vector))


def find_ax_lig(donor_vec, df):
    max = donor_vec.index.max()
    angles = np.zeros((max, max))
    i = 0
    while i < max:
        j = 0
        vec_i = donor_vec.iloc[i]
        norm_i = np.linalg.norm(vec_i)
        while j < max:
            vec_j = donor_vec.iloc[j]
            norm_j = np.linalg.norm(vec_j)
            dot = np.dot(vec_i, vec_j)
            angle_rad = cmath.acos(dot / (norm_i * norm_j))
            angle_rad = angle_rad.real
            angle_deg = round(math.degrees(angle_rad), 2)
            angles[i][j] = angle_deg
            j += 1
        i += 1
    # TODO fix this section to be by column name and not by index.
    angle_limit = 150
    # print(angles)
    check_array = np.all(angles < angle_limit, axis=0)

    while len(np.where(check_array == True)[0]) > 1:
        check_array = np.all(angles < angle_limit, axis=0)
        angle_limit -= 5
    axial_lig_idx = np.where(check_array == True)
    axial_lig_idx = int(axial_lig_idx[0]) + 1  # Get the int of the vector that
    coord_vector = df.iloc[axial_lig_idx, 2:5] * -1 # TODO update this to not be inverse
    donor_id = df.iloc[axial_lig_idx, 1]
    # print(coord_vector, donor_id)
    return (coord_vector, donor_id)


def get_bond_length(xyzdf, atom_idx, other_idx):
    pass


def total_charge(
        ox_state: int,
        lig_charge: int,
        structure: str or int = None
):
    """Return the total charge from sum of oxidation state, ligand charge, and presence of a ligand in structure.

    Currently only 'hyd' - hydride enabled as -1 from structure. Alternatively, an int can be passed
    and added to oxidation state and ligand charge.
    """

    if 'hyd' == structure:
        hydride_charge = -1
    elif isinstance(structure, int):
        hydride_charge = int(structure)
    else:
        hydride_charge = 0
    return ox_state + lig_charge + hydride_charge

def high_spin_state(d_elec_count: float):
    """"""
    if d_elec_count > 5:
        high_spin_state = 2.5 - ((d_elec_count - 5) * 0.5)
    else:
        high_spin_state = d_elec_count * 0.5
    return high_spin_state

def calc_multiplicities(metal_id, ox_state, ):
    d_elec_count = d_counts[metal_id] - ox_state
    high_spin = high_spin_state(d_elec_count)
    poss_spin_states = []
    while high_spin >= 0:
        poss_spin_states.append(high_spin)
        high_spin -= 1
    poss_multiplicities = [int((2 * spin_state) + 1)
                           for spin_state in poss_spin_states]
    return poss_multiplicities


def submission_script_writer(file_list, path):
    with open(f"{path}\\submission_script.txt", "w+", newline="\n") as sub_file:
        sub_file.write("#!/bin/sh\n")
        for filename in file_list:
            sub_file.write(f"sbatch {filename}.sh\n")

    ### NEEDS UPDATING ###
def get_spin_state(multiplicity):
    """
    Determine previous spin state (i.e. ls, is, hs)
    low spin: i = 0
    intermediate spin i = 1
    high spin i = 2
    These correspond to multiplicities in tuples for each d-count.
    """
    multiplicity = int(multiplicity)
    i = 0
    while multiplicity > 0:
        i += 1
        multiplicity -= 2
    return i - 1

def get_d_count(metal_id, ox_state):
    d_zero_count = d_counts[metal_id]
    d_count = d_zero_count - ox_state
    return d_count

def get_new_multiplicity(metal_id, old_mult, new_ox_state):
    old_spin_state = get_spin_state(old_mult)
    new_d_count = get_d_count(metal_id, new_ox_state)
    if new_d_count > 10:
        return -10
    new_mult_tuple = d_count_mults[f"d{new_d_count}"]
    if isinstance(new_mult_tuple, tuple):
        try:
            new_mult = new_mult_tuple[old_spin_state]
        except IndexError:
            return new_mult_tuple[-1]
    else:
        new_mult = new_mult_tuple
    return new_mult