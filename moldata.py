from align import group_by

def first_order_neighbour_elements(data):
    return [ reduce(lambda x,y:x+y, sorted(map(lambda con: data['atoms'][con]['type'], conn_list)), '') for conn_list in [ atom['conn'] for index, atom in data['atoms'].items() ] ]

def second_order_neighbour_elements(data):
    return [ reduce(lambda x,y:x+y, sorted([ data['atoms'][second_neighbour_index]['type'] for second_neighbour_index in reduce(lambda x,y: x+y, [ data['atoms'][conn]['conn'] for conn in conn_list ])])) for conn_list in [ atom['conn'] for index, atom in data['atoms'].items() ] ]

def equivalence_list(data):
    return split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data['atoms'].items() ])

def flavour_list(data):
    eq_list = equivalence_list(data)
    first_neighbour_list, second_neighbour_list = first_order_neighbour_elements(data), second_order_neighbour_elements(data)
    grouped_eq_list = group_by(equivalence_list(data), lambda x:x)
    return [ str(len(grouped_eq_list[eq])) + '|' + first_neighbours + '|' + second_neighbours for eq, first_neighbours, second_neighbours in zip(eq_list, first_neighbour_list, second_neighbour_list) ]

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
