from align import group_by

def atom_types_for_indexes(indexes_list, atoms):
    return [ [atoms[index]['type'] for index in indexes] for indexes in indexes_list ]

def connected_atom_indexes_for_indexes(indexes_list, atoms):
    return reduce(lambda x,y: x+y, [atoms[index]['conn'] for index in indexes_list], [])

def sorted_atom_types_for_indexes(indexes_list, atoms):
    return [sorted(indexes) for indexes in atom_types_for_indexes(indexes_list, atoms) ]

def joined_sorted_atom_types_for_indexes(indexes_list, atoms):
    return [ ''.join(indexes) for indexes in  sorted_atom_types_for_indexes(indexes_list, atoms) ]

def nth_order_neighbour_elements(data, n):
    atoms = data['atoms']
    nth_order_indexes = {}
    nth_order_indexes[0] = [ [index] for index in atoms.keys()]
    for i in range(1, n+1):
        nth_order_indexes[i] = [ connected_atom_indexes_for_indexes(indexes_list, atoms) for indexes_list in nth_order_indexes[i-1] ]
    return [(k, joined_sorted_atom_types_for_indexes(v, atoms)) for (k, v) in nth_order_indexes.items() if k != 0]

def equivalence_list(data):
    return split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data['atoms'].items() ])

def flavour_list(data):
    eq_list = equivalence_list(data)
    nth_neighbour_list = [a[1] for a in nth_order_neighbour_elements(data, 3)]
    grouped_eq_list = group_by(equivalence_list(data), lambda x:x)
    return [ str(len(grouped_eq_list[eq])) + '|' + a + '|' + b + '|' + c for eq, a, b, c in zip(eq_list, *nth_neighbour_list ) ]

def element_list(data):
    return [ atom['type'] for index, atom in data['atoms'].items() ]

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
