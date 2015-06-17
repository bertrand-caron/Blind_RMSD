import numpy as np
from Vector import Vector, rotmat, m2rotaxis
from itertools import product, groupby, permutations
from functools import partial
from charnley_rmsd import kabsch
from copy import deepcopy
from scoring import rmsd_array, ad_array, rmsd, ad, rmsd_array_for_loop
from permutations import N_amongst_array
import pprint
from ChemicalPoint import ChemicalPoint, on_elements, on_coords, on_canonical_rep, ELEMENT_NUMBERS

pp = pprint.PrettyPrinter(indent=2)

on_self, on_first_element, on_second_element = lambda x:x, lambda x:x[0], lambda x:x[1]
on_third_element, on_fourth_element = lambda x: x[2], lambda x: x[3]
on_second_element_and_flavour = lambda grouped_flavours, x: str(x[1]) + str(len(grouped_flavours[ x[2] ]))

# Kabsch Algorithm options
DEFAULT_MIN_N_UNIQUE_POINTS = 4
MAX_N_COMPLEXITY = 10 # Maximum number of permutations is MAX_N_COMPLEXITY^(N_UNIQUE_POINTS - MIN_N_UNIQUE_POINTS)

ALLOW_SHORTCUTS = False
DEFAULT_SCORE_TOLERANCE = 0.01

ON_BOTH_LISTS = [0,1]

DISABLE_BRUTEFORCE_METHOD = True

# Align points on points
def pointsOnPoints(point_lists, silent=True, use_AD=False, element_lists=None, flavour_lists=None, show_graph=False, bonds=None, score_tolerance=DEFAULT_SCORE_TOLERANCE, soft_fail=False):

    # Initializers
    has_elements = True if element_lists and all(element_lists) else False
    has_flavours = True if flavour_lists and all(flavour_lists) else False
    if not has_flavours:
        flavour_lists = map(lambda a_list: range(len(a_list)), element_lists)
        has_flavours = True
    has_bonds = True if bonds else False
    if has_bonds: bonds = map(np.array, bonds)

    # Assert that the fitting make sense
    assert len(point_lists[0]) == len(point_lists[1]), "Error: Size of point lists doesn't match: {0} and {1}".format(*map(len, point_lists))
    if has_flavours:
        assert len(flavour_lists[0]) == len(flavour_lists[1]), "Error: Size of flavour lists doesn't match: {0} and {1}".format(*map(len, flavour_lists))
        assert len(flavour_lists[0]) == len(point_lists[1]), "Error: Size of flavour lists doesn't match size of point lists: {0} and {1}".format(*map(len, [flavour_lists[0], point_lists[1]]))
        get_sorted_eqgroup_lengths = lambda flavour_list: sorted(map(len, group_by(flavour_list, on_self).values()))
        assert get_sorted_eqgroup_lengths(flavour_lists[0]) == get_sorted_eqgroup_lengths(flavour_lists[1]), "Error: There is not a one to one mapping between the lengths of the flavour sets: {0} and {1}".format(*map(get_sorted_eqgroup_lengths, flavour_lists))
    if has_elements:
        assert len(element_lists[0]) == len(element_lists[1]), "Error: Size of element lists doesn't match: {0} and {1}".format(*map(len, element_lists))
        assert len(element_lists[0]) == len(point_lists[1]), "Error: Size of element lists doesn't match size of point lists: {0} and {1}".format(*map(len, [element_lists[0], point_lists[1]]))
        assert sorted(element_lists[0]) == sorted(element_lists[1]), "Error: There is not a one to one mapping of the element sets: {0} and {1}".format(*map(sorted, element_lists))
    if has_bonds:
        assert( bonds[0].shape == tuple(map(len, point_lists)) ), "Error: Bonds array does have have the expected shape: {0} != {1}".format(bonds[0].shape, map(len, point_lists))
        assert( bonds[1].shape == tuple(map(len, point_lists)) ), "Error: Bonds array does have have the expected shape: {0} != {1}".format(bonds[1].shape, map(len, point_lists))

    point_arrays = map(np.array, point_lists)
    center_of_geometries = map(center_of_geometry, point_arrays)
    if has_elements:
        grouped_flavours_lists = map(lambda index:group_by(flavour_lists[index], lambda x:x ), ON_BOTH_LISTS)
        chemical_points_lists = get_chemical_points_lists(point_lists, element_lists, flavour_lists, has_flavours, grouped_flavours_lists)
        mask_array = np.zeros((len(flavour_lists[0]), len(flavour_lists[0])))
        dumb_array = np.chararray((len(flavour_lists[0]), len(flavour_lists[0])), itemsize=10)
        for i, chemical_point0 in enumerate(chemical_points_lists[0]):
            for j, chemical_point1 in enumerate(chemical_points_lists[1]):
                mask_array[i, j] = 0. if chemical_point0.canonical_rep == chemical_point1.canonical_rep else 1.0E5
                dumb_array[i, j] = "{0} {1}".format(chemical_point0.canonical_rep, chemical_point1.canonical_rep)
        if not silent: print chemical_points_lists
        if not silent: print dumb_array
    distance_function, distance_array_function = rmsd if not use_AD else ad, lambda *args, **kwargs: rmsd_array_for_loop(*args, mask_array=mask_array if mask_array is not None else None, **kwargs) if not use_AD else ad_array

    # First, remove translational part from both by putting the center of geometry in (0,0,0)
    centered_point_arrays = [point_arrays[0] - center_of_geometries[0], point_arrays[1] - center_of_geometries[1]]

    # Assert than the center of geometry of the translated point list are now on (0,0,0)
    [ assert_array_equal( center_of_geometry(a), np.array([0,0,0])) for a in centered_point_arrays ]

    # Break now if the molecule has less than 3 atoms
    if len(point_lists[0]) < 3 :
        return (centered_point_arrays[0]).tolist(), 0.0

    # Break now if there are no rotational component
    if distance_function(*centered_point_arrays) <= score_tolerance and ALLOW_SHORTCUTS:
        if not silent: print "Info: A simple translation was enough to match the two set of points. Exiting successfully."
        assert_found_permutation(*centered_point_arrays, silent=silent)
        return (centered_point_arrays[0] + center_of_geometries[1]).tolist(), distance_function(*centered_point_arrays)

    method_results = {}
    if not DISABLE_BRUTEFORCE_METHOD:
        # Try the bruteforce method first
        method_results['bruteforce'] = bruteforce_aligning_vectors_method(point_arrays, distance_array_function=distance_array_function, score_tolerance=score_tolerance, silent=silent and True)
        method_results['lucky_kabsch'] = lucky_kabsch_method(point_lists, element_lists, flavour_lists=flavour_lists, distance_array_function=distance_array_function, score_tolerance=score_tolerance, show_graph=show_graph, silent=silent)
        method_results['bruteforce_kabsch'] = bruteforce_kabsch_method(point_lists, element_lists, flavour_lists=flavour_lists, distance_array_function=distance_array_function, score_tolerance=score_tolerance, show_graph=show_graph, silent=silent)

    # Try the flavoured Kabsch method if we have elements
    if has_elements:
        method_results['kabsch'] = flavoured_kabsch_method(point_lists, element_lists, flavour_lists=flavour_lists, distance_array_function=distance_array_function, score_tolerance=score_tolerance, show_graph=show_graph, silent=silent)

    best_method = sorted(method_results.items(), key=lambda x:x[1]['score'] if 'score' in x[1] else 100.)[0][0]
    best_match = method_results[best_method]['array']
    if best_match == None:
        if not soft_fail: raise Exception("Best match is None. Something went wrong.")
        else: return None, float('inf')
    
    if not silent: print "Info: Scores of methods are: {0}".format(dict([ (k, v['score']) for (k,v) in method_results.items() if 'score' in v]))
    if not silent: print "Info: Best score was achieved with method: {0}".format(best_method)
    
    corrected_best_match = best_match - center_of_geometry(best_match) + center_of_geometries[1]
    assert_array_equal(*map(center_of_geometry, [corrected_best_match, point_arrays[1]]), message="{0} != {1}")
    assert_found_permutation(corrected_best_match, point_arrays[1], silent=silent)
    
    return corrected_best_match.tolist(), method_results[best_method]['score']

### METHODS ###

def bruteforce_aligning_vectors_method(centered_arrays, distance_array_function=rmsd_array, silent=True, score_tolerance=DEFAULT_SCORE_TOLERANCE):
    silent=True
    # First, select our first point on the translated structure; it is mandatory that this point is not on the center of geometry
    reference_vectors = [None, None]
    for point in centered_arrays[0][:,0:3]:
        reference_vectors[0] = Vector(point)
        break

    # Then try all the rotation that put one on the atom of the first set into one of the atoms of the second sets
    # There are N such rotations
    best_match, best_score = centered_arrays[0], distance_array_function(*centered_arrays, silent=silent)
    point_arrays = [None, None]
    for point_arrays[0] in centered_arrays[0][:,0:3]:
        for point_arrays[1] in centered_arrays[1][:,0:3]:

            reference_vectors[1] = Vector(point_arrays[1])

            r = rotmat(*reversed(reference_vectors))
            if not silent: print "    Info: Rotation parameters: {0} deg, axis {1}".format(m2rotaxis(r)[0]*180/np.pi, m2rotaxis(r)[1])
            assert m2rotaxis(r)[0] != 180., "Error: 180 degree rotation matrix currently not working"

            rotated_point_arrays = [np.dot(centered_arrays[0], r)]

            # If the norm of the vector are the same, check that the rotation effectively put p on q
            if all([ vect.norm() != 0. for vect in reference_vectors]):
                if reference_vectors[1].norm() == reference_vectors[0].norm():
                    assert_array_equal(rotated_point_arrays[0][0, 0:3], reference_vectors[1]._ar)
                # Else do the same operation on the normalized vectors
                else:
                    assert_array_equal(Vector(rotated_point_arrays[0][0, 0:3]).normalized()._ar, reference_vectors[1].normalized()._ar)

            current_score = distance_array_function(rotated_point_arrays[0], centered_arrays[1], silent=silent)
            if current_score <= best_score: 
                best_match, best_score = rotated_point_arrays[0], current_score

            if best_score <= score_tolerance and ALLOW_SHORTCUTS:
                if not silent: print "    Info: Found a really good match (Score={0}) worth aborting now. Exiting successfully.".format(best_score)
                break
        # Only iterate over the first point of centered_arrays[0]
        break
    
    if not silent: print "    Info: Minimum Score from bruteforce algorithm is: {0}".format(best_score)
    return {'array': best_match.tolist(), 'score': best_score, 'reference_array': centered_arrays[1]}


def get_chemical_points_lists(point_lists, element_lists, flavour_lists, has_flavours, grouped_flavours_lists):
    N_points = len(point_lists[0])
    chemical_points_lists = map(lambda index:[ChemicalPoint(*zipped_point, **{'grouped_flavours': grouped_flavours_lists[index]}) for zipped_point in zip(point_lists[index],
                                                                                                                                                   range(N_points),
                                                                                                                                                   element_lists[index],
                                                                                                                                                   flavour_lists[index] if has_flavours else [None] * N_points
                                                                                                                                                   )
                                              ],
                                ON_BOTH_LISTS)
    return chemical_points_lists

def flavoured_kabsch_method(point_lists, element_lists, silent=True, distance_array_function=rmsd_array, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE):
    point_arrays = map(np.array, point_lists)
    has_flavours= bool(flavour_lists)
    MIN_N_UNIQUE_POINTS = DEFAULT_MIN_N_UNIQUE_POINTS if len(point_lists[0]) >= DEFAULT_MIN_N_UNIQUE_POINTS else 3
    
    if not silent: print "    Info: Found element types. Trying flavoured {0}-point Kabsch algorithm on flavoured elements types ...".format(MIN_N_UNIQUE_POINTS)
    
    grouped_flavours_lists = map(lambda index:group_by(flavour_lists[index], lambda x:x ), ON_BOTH_LISTS)
    chemical_points_lists = get_chemical_points_lists(point_lists, element_lists, flavour_lists, has_flavours, grouped_flavours_lists)
    
    grouped_chemical_points_lists = map(lambda chemical_points:group_by(chemical_points, on_canonical_rep), chemical_points_lists)
    # Try to find MIN_N_UNIQUE_POINTS unique elements type points
    unique_points_lists = map(lambda grouped_chemical_points: [group[0] for group in grouped_chemical_points.values() if len(group) == 1], grouped_chemical_points_lists)
    # Order them by decreasing element type
    unique_points_lists = map(lambda index: sorted(unique_points_lists[index], key=lambda x: ELEMENT_NUMBERS[ x.element.upper() ], reverse=True),
                                     ON_BOTH_LISTS)

    assert len(unique_points_lists[0]) == len(unique_points_lists[1]), "Error: Non matching number of unique points in {0} and {1}".format(*unique_points_lists)

    if not silent: print "    Info: Unique groups found based on element types: {0}".format(unique_points_lists[0])

    if len(unique_points_lists[0]) < MIN_N_UNIQUE_POINTS:
        if not silent: print "    Warning: Unable to find at least {N} unique point with the elements provided. Trying to disambiguate enough points to make a fit.".format(N=MIN_N_UNIQUE_POINTS)

        missing_points = MIN_N_UNIQUE_POINTS - len(unique_points_lists[0])

        ambiguous_point_groups = map(lambda grouped_chemical_points: sorted([group for group in grouped_chemical_points.values() if 1 < len(group) <= MAX_N_COMPLEXITY ], key=len),
                                     grouped_chemical_points_lists )

        # Order them by number of heaviest atoms first
        ambiguous_point_groups = map(lambda index: sorted(ambiguous_point_groups[index], key=lambda x: ELEMENT_NUMBERS[ x[0].element.upper() ], reverse=True),
                                     ON_BOTH_LISTS)

        N_ambiguous_points = sum( map(len, ambiguous_point_groups[0]))

        if N_ambiguous_points < missing_points:
            if not silent: print "    Error: Couldn'd find enough point to disambiguate: {M} (unique points) + {P} (ambiguous points) < {N} (required points). Returning best found match ...".format(P=N_ambiguous_points, M=len(unique_points_lists[0]), N=MIN_N_UNIQUE_POINTS)
            return {'array': None, 'score': None, 'reference_array': point_arrays[1]}

        if not silent: print "    Info: Found enough point ({N}) to disambiguate. Trying kabsch algorithm ...".format(N=N_ambiguous_points)

        permutations_list = []
        atom_indexes = lambda chemical_points: map(lambda chemical_point: chemical_point.index, chemical_points)

        # For each ambiguous group
        ambiguous_points = 0
        N_list = []
        for group in ambiguous_point_groups[0]:
            if ambiguous_points == missing_points:
                N_list.append( 0 )
            else:
                N_list.append( min(len(group), missing_points-ambiguous_points) )
                ambiguous_points += N_list[-1]

        if not silent: print "    Info: Ambiguous groups are: {0} (number of points taken in each group: {1})".format(ambiguous_point_groups[0], N_list)

        permutation_lists = map(lambda group, N: permutations(atom_indexes(group), r=N),
                                ambiguous_point_groups[0],
                                N_list)

        complete_permutation_list = product(*permutation_lists)
        #print list(complete_permutation_list)

        for group, N in zip(ambiguous_point_groups[1], N_list):
            unique_points_lists[1] += group[0:N]

        best_match, best_score = None, None
        for group_permutations in complete_permutation_list:
            ambiguous_unique_points = [deepcopy(unique_points_lists[0])]
            for group in group_permutations:
                for index in group:
                    new_point = chemical_points_lists[0][index]
                    ambiguous_unique_points[0].append(new_point)

            if not silent: print '        Info: Attempting a fit between points {0} and {1}'.format(ambiguous_unique_points[0], unique_points_lists[1])
            assert( map(on_elements, ambiguous_unique_points[0]) == map(on_elements, unique_points_lists[1])), "Error: Trying to match points whose elements don't match: {0} != {1}".format(map(on_elements, ambiguous_unique_points[0]), map(on_elements, unique_points_lists[1]))
            P, Q = map(on_coords, ambiguous_unique_points[0]), map(on_coords, unique_points_lists[1])
            U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)
            kabsched_list1 = np.dot(point_arrays[0]-Pc, U) + Qc
            current_score = distance_array_function(kabsched_list1, point_arrays[1], silent=silent)
            if (not best_score) or current_score <= best_score:
                best_match, best_score = kabsched_list1, current_score
                if not silent: print "    Info: Best score so far with random {0}-point Kabsch fitting: {1}".format(MIN_N_UNIQUE_POINTS, best_score)
                if current_score <= score_tolerance: return {'array': best_match.tolist(), 'score': best_score, 'reference_array': point_arrays[1]}
            if show_graph: do_show_graph([(P-Pc,"P-Pc"), (Q-Qc, "Q-Qc"), (point_arrays[0] - Pc, "P1-Pc"), (point_arrays[1] - Qc, "P2-Qc")])

        if not silent: print "    Info: Returning best match with random {0}-point Kabsch fitting (Score: {1})".format(MIN_N_UNIQUE_POINTS, best_score)
        return {'array': best_match.tolist(), 'score': best_score, 'reference_array': point_arrays[1]}
    else:
        assert map(on_elements, unique_points_lists[0][0:MIN_N_UNIQUE_POINTS]) == map(on_elements, unique_points_lists[1][0:MIN_N_UNIQUE_POINTS]), "Error: Unique points have not been ordered properly: {0} and {1}".format(map(on_elements, unique_points_lists[0][0:MIN_N_UNIQUE_POINTS]), map(on_elements, unique_points_lists[1][0:MIN_N_UNIQUE_POINTS]))
        
        # Align those MIN_N_UNIQUE_POINTS points using Kabsch algorithm
        P, Q = map(on_coords, unique_points_lists[0]), map(on_coords, unique_points_lists[1])
        U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)
        kabsched_list1 = np.dot(point_arrays[0]-Pc, U) + Qc

        if show_graph: do_show_graph([(kabsched_list1, "P1_kabsch"), (point_arrays[1], "P2")])

        current_match, current_score = kabsched_list1, distance_array_function(kabsched_list1, point_arrays[1], silent=silent)
        
        if not silent: print "    Info: Klabsch algorithm on unique element types found a better match with a Score of {0}".format(current_score)
        return {'array': current_match.tolist(), 'score': current_score, 'reference_array': point_arrays[1]}

def lucky_kabsch_method(point_lists, element_lists, silent=True, distance_array_function=rmsd_array, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE):
    point_arrays = map(np.array, point_lists)

    P, Q = point_arrays
    U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)
    kabsched_list1 = np.dot(point_arrays[0]-Pc, U) + Qc

    current_match, current_score = kabsched_list1, distance_array_function(kabsched_list1, point_arrays[1], silent=silent)
    if not silent: print "    Info: Minimum Score from lucky Kabsch method is: {0}".format(current_score)
    return {'array': current_match.tolist(), 'score': current_score}

def bruteforce_kabsch_method(point_lists, element_lists, silent=True, distance_array_function=rmsd_array, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE):
    N_BRUTEFORCE_KABSCH = 4
    silent = True

    point_arrays = map(np.array, point_lists)

    unique_points = [None, point_arrays[1][0:N_BRUTEFORCE_KABSCH, 0:3] ]
    best_match, best_score = None, None

    for permutation in N_amongst_array(point_arrays[0], N_BRUTEFORCE_KABSCH):
        unique_points[0] = map(lambda index: point_arrays[0][index, 0:3], permutation)
        P, Q = unique_points
        U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)
        kabsched_list1 = np.dot(point_arrays[0]-Pc, U) + Qc

        current_match, current_score = kabsched_list1, distance_array_function(kabsched_list1, point_arrays[1], silent=silent)

        current_score = distance_array_function(kabsched_list1, point_arrays[1], silent=silent)
        if (not best_score) or current_score <= best_score:
            best_match, best_score = kabsched_list1, current_score
            if not silent: print "    Info: Best score so far with bruteforce {N}-point Kabsch fitting: {best_score}".format(best_score=best_score, N=N_BRUTEFORCE_KABSCH)
    if not silent: print "    Info: Minimum Score from bruteforce Kabsch method is: {0}".format(best_score)
    return {'array': current_match.tolist(), 'score': best_score}

#################
#### HELPERS ####
#################

def assert_array_equal(array1, array2, message="{0} and {1} are different"):
    assert np.allclose( array1, array2, atol=1e-5), message.format(array1, array2)

def assert_found_permutation(array1, array2, silent=True, hard_fail=False):
    perm_list = []
    for i, point1_array in enumerate(array1[:,0:3]):
        min_dist, min_index = None, None
        for j, point2_array in enumerate(array2[:,0:3]):
            distance = np.linalg.norm(point1_array-point2_array)
            min_dist = min(min_dist, distance) if min_dist else distance
            if distance == min_dist: min_index = j
        perm_list.append((i, min_index))

    offending_indexes = filter(lambda x: len(x[1])>=2, [ (value, list(group)) for value, group in groupby(perm_list, lambda x:x[1]) ])
    #ambiguous_indexes = list( set(zip(*perm_list)[0]) - set(zip(*perm_list)[1]) ) + [value for value, group in offending_indexes]

    # Assert that perm_list is a permutation, i.e. that every obj of the first list is assigned one and only once to an object of the second list
    if hard_fail: 
        assert sorted(zip(*perm_list)[1]) == list(zip(*perm_list)[0]), "Error: {0} is not a permutation of {1}, which means that the best fit does not allow an unambiguous one-on-one mapping of the atoms. The method failed.".format(sorted(zip(*perm_list)[1]), zip(*perm_list)[0])
        if not silent: print "Info: {0} is a permutation of {1}. This is a good indication the algorithm might have succeeded.".format(zip(*perm_list)[1], zip(*perm_list)[0])
    else:
        if not sorted(zip(*perm_list)[1]) == list(zip(*perm_list)[0]): 
            if not silent: print "Error: {0} is not a permutation of {1}, which means that the best fit does not allow an unambiguous one-on-one mapping of the atoms. The method failed.".format(sorted(zip(*perm_list)[1]), zip(*perm_list)[0])
        else:
            if not silent: print "Info: {0} is a permutation of {1}. This is a good indication the algorithm might have succeeded.".format(zip(*perm_list)[1], zip(*perm_list)[0])

def rotation_matrix_kabsch_on_points(points1, points2):
    # Align those points using Kabsch algorithm
    P, Q = np.array(points1), np.array(points2)
    #print P
    #print Q
    Pc, Qc = kabsch.centroid(P), kabsch.centroid(Q)
    P, Q = P - Pc, Q - Qc
    U = kabsch.kabsch(P, Q)
    return U, Pc, Qc

def do_show_graph(array_list):
    import plot3D as p
    symbol_list = ['x', 'o', '+', '^']
    colour_list = ['b', 'r', 'g', 'y']
    for symbol, colour, array in zip(symbol_list, colour_list, array_list):
        p.plotPoints(array[0], colour,  symbol, array[1])
    p.showGraph()

def distance(point1, point2):
    return np.linalg.norm(point1 - point2)

def center_of_geometry(point_array):
    return np.mean(point_array, axis=0)

# A reimplemetation of python crappy itertool's groupby method with dictionnaries and less BS
def group_by(iterable, key):
    group_dict = {}
    for obj in iterable:
        group_dict.setdefault( key(obj), [])
        group_dict[key(obj)].append(obj)
    return group_dict
