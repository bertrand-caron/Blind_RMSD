from typing import Dict, Callable, Sequence, Any, List
from functools import reduce

from chemistry_helpers.pdb import substitute_coordinates_in, is_pdb_atom_line

def should_keep_atom(atom, united=False):
    return (united and 'uindex' in atom) or (not united)

def atom_types_and_indices_for_indexes(indexes_list, atoms, united=False):
    return [
        [
            atoms[index]['type'] + str(len(atoms[index]['conn']))
            for index in indexes
            if should_keep_atom(atoms[index], united)
        ]
        for indexes in indexes_list
    ]

def connected_atom_indexes_for_indexes(indexes_list, atoms, united=True):
    return reduce(
        lambda x, y: x + y,
        [
            atoms[index]['conn']
            for index in indexes_list
            if should_keep_atom(atoms[index], united)
        ],
        [],
    )

def sorted_atom_types_for_indexes(indexes_list, atoms, united=False):
    return [
        sorted(indexes)
        for indexes in atom_types_and_indices_for_indexes(indexes_list, atoms, united)
    ]

def joined_sorted_atom_types_for_indexes(indexes_list, atoms, united=False):
    return [
        ','.join(indexes)
        for indexes in sorted_atom_types_for_indexes(indexes_list, atoms, united)
    ]

def nth_order_neighbour_elements(data, n, united=False):
    atoms = data['atoms']
    nth_order_indexes = []
    nth_order_indexes.append(
        [
            [index]
            for (index, atom) in list(atoms.items())
            if should_keep_atom(atom, united) ]
    )

    for i in range(1, n + 1):
        nth_order_indexes.append([
            connected_atom_indexes_for_indexes(indexes_list, atoms, united)
            for indexes_list in nth_order_indexes[i-1]
        ])

    return [
        (i, joined_sorted_atom_types_for_indexes(v, atoms, united))
        for (i, v) in enumerate(nth_order_indexes)
    ]

def equivalence_list(data, united=False):
    return split_equivalence_group([ atom['equivalenceGroup'] for index, atom in list(data['atoms'].items()) if should_keep_atom(atom, united) ])

FLAVOUR_LIST_SHELL_NUMBER = 4

def flavour_list(data, united=False):
    eq_list = equivalence_list(data, united)
    nth_neighbour_list = [a[1] for a in nth_order_neighbour_elements(data, FLAVOUR_LIST_SHELL_NUMBER, united)]
    grouped_eq_list = group_by(equivalence_list(data), lambda x:x)

    flavours =  [
        (
            '|'.join(
                nth_neighbours[0:FLAVOUR_LIST_SHELL_NUMBER] + ('EQ{0}'.format(len(grouped_eq_list[eq])),)
            ),
        )
        for (eq, nth_neighbours) in zip(
            eq_list,
            list(zip(*nth_neighbour_list)),
        )
    ]

    return flavours

def element_list(data, united=False):
    return [ atom['type'] for index, atom in list(data['atoms'].items()) if should_keep_atom(atom, united) ]

def pdb_lines(data, united=False):
    return [atom['pdb'] for index, atom in list(data['atoms'].items()) if should_keep_atom(atom, united) ]

def connect_lines(data):
    return ['CONECT{0:5d}{1:5d}'.format(*bond['atoms']) for bond in data['bonds']]

def pdb_str(data, united=False):
    return '\n'.join(pdb_lines(data, united))

def nm_to_A(x):
    return 10*x

def point_list(data, united=False):
    return [ list(map(nm_to_A, atom['ocoord'] if 'ocoord' in atom else atom['coord'])) for index, atom in list(data['atoms'].items()) if should_keep_atom(atom, united) ]

def united_hydrogens_point_list(data, united=False):
    return [ list(map(nm_to_A, atom['ocoord'] if 'ocoord' in atom else atom['coord'])) for index, atom in list(data['atoms'].items()) if not should_keep_atom(atom, united) ]

def get_united_hydrogens_pdb_lines(data, united=False):
    return [ atom['pdb'] for _, atom in list(data['atoms'].items()) if not should_keep_atom(atom, united) ]

def permutated_list(a_list, permutation):
    on_j = lambda x:x[1]

    new_list = []
    assert len(a_list) == len(permutation), '{0} != {1}'.format(a_list, permutation)
    for (i,j) in sorted(permutation, key=on_j):
        new_list.append(a_list[i])
    return new_list

def map_to_str(a_list):
    return [str(x) for x in a_list]

def aligned_pdb_str(data, alignment, united=False):
    from io import StringIO
    pdb_str = StringIO()

    alignment_coordinates, _, united_H_coordinates, final_permutation = alignment
    heavy_atoms_pdb_lines = pdb_lines(data, united)

    if final_permutation:
        heavy_atoms_pdb_lines, alignment_coordinates = list(map(
            lambda a_list: permutated_list(a_list, final_permutation),
            (heavy_atoms_pdb_lines, alignment_coordinates),
        ))

    assert len(heavy_atoms_pdb_lines) == len(alignment_coordinates)
    for line, coordinates in zip(heavy_atoms_pdb_lines, alignment_coordinates):
        print(substitute_coordinates_in(line, tuple(coordinates)), file=pdb_str)

    atom_count = 0
    for line in get_united_hydrogens_pdb_lines(data, united):
        if is_pdb_atom_line(line):
            print(substitute_coordinates_in(
                line,
                tuple(united_H_coordinates[atom_count]),
            ), file=pdb_str)
            atom_count += 1

    for connect_line in connect_lines(data):
        print(connect_line, file=pdb_str)

    return pdb_str.getvalue()

# Differentiate -1's
def split_equivalence_group(eq_list):
    accu = 0
    split_eq_list = []
    for eq in eq_list:
        if eq != -1: split_eq_list.append(eq)
        else:
            split_eq_list.append(eq-accu)
            accu += 1
    return split_eq_list

def group_by(iterable: Sequence[Any], key: Callable[[Any], Any]) -> Dict[Any, List[Any]]:
    group_dict = {}
    for obj in iterable:
        group_dict.setdefault(key(obj), [])
        group_dict[key(obj)].append(obj)
    return group_dict

if __name__ == '__main__':
    from Blind_RMSD.pdb import pdb_data_for
    data = pdb_data_for(open('data/1.pdb').read())
    print(data.flavour_lists)
