import numpy as np
from Vector import Vector, rotmat, m2rotaxis
import itertools
from functools import partial
from charnley_rmsd import kabsch
from copy import deepcopy
from scoring import rmsd_array, ad_array, rmsd, ad

on_self, on_first_element, on_second_element = lambda x:x, lambda x:x[0], lambda x:x[1]
on_second_element_and_flavour = lambda grouped_flavours, x: str(x[1]) + str(len(grouped_flavours[ x[2] ]))

MAX_N_COMPLEXITY = 3

ALLOW_SHORTCUTS = False

# Align points on points
def pointsOnPoints(point_lists, silent=True, use_AD=False, element_lists=None, flavour_lists=None, show_graph=False, score_tolerance=1E-3):

    # Initializers
    distance_function, distance_array_function = rmsd if not use_AD else ad, rmsd_array if not use_AD else ad_array
    has_flavours = True if all(flavour_lists) else False
    has_elements = True if all(element_lists) else False
    if not silent: print # Add a commencing newline

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

    point_arrays = map(np.array, point_lists)
    center_of_geometries = map(center_of_geometry, point_arrays)

    # First, remove translational part from both by putting the center of geometry in (0,0,0)
    centered_point_arrays = [point_arrays[0] - center_of_geometries[0], point_arrays[1] - center_of_geometries[1]]

    # Assert than the center of geometry of the translated point list are now on (0,0,0)
    [ assert_array_equal( center_of_geometry(a), np.array([0,0,0])) for a in centered_point_arrays ]

    # Break now if there are no rotational component
    if distance_function(*centered_point_arrays) <= score_tolerance and ALLOW_SHORTCUTS:
        if not silent: print "Info: A simple translation was enough to match the two set of points. Exiting successfully."
        assert_found_permutation(*centered_point_arrays, silent=silent)
        return (centered_point_arrays[0] + center_of_geometries[1]).tolist()

    method_results = {}
    # Try the bruteforce method first
    method_results['bruteforce'] = bruteforce_aligning_vectors_method(point_arrays, distance_array_function=distance_array_function, silent=silent and False)

    # Try the flavoured Kabsch method if we have elements
    if has_elements:
        method_results['kabsch'] = flavoured_kabsch_method(point_lists, element_lists, flavour_lists=flavour_lists, distance_array_function=distance_array_function, show_graph=show_graph, silent=silent)

    best_method = "kabsch" if method_results['kabsch']['score'] <= method_results['bruteforce']['score'] else "bruteforce"
    best_match = method_results[best_method]['array']
    
    if not silent: print "Info: Scores of methods are: {0}".format(dict([ (k, v['score']) for (k,v) in method_results.items() ]))
    if not silent: print "Info: Best score was achieved with method: {0}".format(best_method)
    
    corrected_best_match = best_match - center_of_geometry(best_match) + center_of_geometries[1]
    assert_array_equal(*map(center_of_geometry, [corrected_best_match, point_arrays[1]]), message="{0} != {1}")
    assert_found_permutation(corrected_best_match, point_arrays[1], silent=silent)
    
    return corrected_best_match.tolist()

### METHODS ###

def bruteforce_aligning_vectors_method(centered_arrays, distance_array_function=rmsd_array, silent=True, score_tolerance=0.01):
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
            if not silent: print "Info: Rotation parameters: {0} deg, axis {1}".format(m2rotaxis(r)[0]*180/np.pi, m2rotaxis(r)[1])
            assert m2rotaxis(r)[0] != 180., "Error: 180 degree rotation matrix currently not working"

            rotated_point_arrays = [np.dot(centered_arrays[0], r)]

            # If the norm of the vector are the same, check that the rotation effectively put p on q
            if reference_vectors[1].norm() == reference_vectors[0].norm():
                assert_array_equal(rotated_point_arrays[0][0, 0:3], reference_vectors[1]._ar)
            # Else do the same operation on the normalized vectors
            else:
                assert_array_equal(Vector(rotated_point_arrays[0][0, 0:3]).normalized()._ar, reference_vectors[1].normalized()._ar)

            current_score = distance_array_function(rotated_point_arrays[0], centered_arrays[1], silent=silent)
            if current_score <= best_score: 
                best_match, best_score = rotated_point_arrays[0], current_score

            if best_score <= score_tolerance and ALLOW_SHORTCUTS:
                if not silent: print "Info: Found a really good match (Score={0}) worth aborting now. Exiting successfully.".format(best_score)
                break
        # Only iterate over the first point of centered_arrays[0]
        break
    
    if not silent: print "Info: Minimum Score from unflavoured algorithm is: {0}".format(best_score)
    return {'array': best_match.tolist(), 'score': best_score}

def flavoured_kabsch_method(point_lists, element_lists, silent=True, distance_array_function=rmsd_array, flavour_lists=None, show_graph=False):
    has_flavours = True if flavour_lists else False
    point_arrays = map(np.array, point_lists)
    if not silent: print "Info: Found element types. Trying flavoured 3-point Kabsch algorithm on flavoured elements types ..."

    # Try to find three unique elements type points
    if has_flavours:
        element_points = map(lambda index: zip(point_lists[index], element_lists[index], flavour_lists[index]), [0,1])
        grouped_flavour_lists = [group_by(flavour_lists[0], on_self), group_by(flavour_lists[1], on_self)]
        grouping_functions = map( lambda index: partial(on_second_element_and_flavour, grouped_flavour_lists[index]) if has_flavours else on_second_element, [0,1])
    else:
        element_points = map(lambda index: zip(point_lists[index], element_lists[index]), [0,1])
        grouping_functions = [on_second_element, on_second_element]

    grouped_element_points = map( lambda index:group_by(element_points[index], grouping_functions[index]), [0,1])
    unique_points = map(lambda grouped_element_point: [group[0] for group in grouped_element_point.values() if len(group) == 1], grouped_element_points)

    assert len(unique_points[0]) == len(unique_points[1]), "Error: Non matching number of unique points in {0} and {1}".format(*unique_points)

    if not silent: print "    Info: Unique groups found based on element types: {0}".format(unique_points[0])

    if len(unique_points[0]) < 3:
        if not silent: print "    Warning: Unable to find at least three unique point with the elements provided. Trying to disambiguate enough points to make a fit."

        missing_points = 3 - len(unique_points[0])

        ambiguous_point_groups = map(lambda grouped_element_point: sorted([group for group in grouped_element_point.values() if 1 < len(group) <= MAX_N_COMPLEXITY ], key=len), grouped_element_points )

        if len(ambiguous_point_groups[0]) <= missing_points:
            if not silent: print "Error: Couldn'd find enough point to disambiguate. Returning best found match ..."
            return None

        if not silent: print "    Info: Found enough point to disambiguate. Trying kabsch algorithm ..."

        unique_points[1] += map(on_first_element, ambiguous_point_groups[1])[0:missing_points]
        permutations_list = itertools.product(*map(range, [len(group) for group in ambiguous_point_groups[0][0:missing_points] ]))
        
        best_match, best_score = None, None
        for permutation in permutations_list:

            ambiguous_unique_points = [deepcopy(unique_points[0])]
            for i, ambiguous_group in enumerate(ambiguous_point_groups[0]):
                ambiguous_unique_points[0].append(ambiguous_group[ permutation[i] ])
                if len(ambiguous_unique_points[0]) == 3: break
            
            # Align those three points using Kabsch algorithm
            print ambiguous_unique_points[0]
            print unique_points[1]
            P, Q = map(on_first_element, ambiguous_unique_points[0]), map(on_first_element, unique_points[1])
            U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)
            kabsched_list1 = np.dot(point_arrays[0]-Pc, U) + Qc
            current_score = distance_array_function(kabsched_list1, point_arrays[1], silent=silent)
            if (not best_score) or current_score <= best_score:
                best_match, best_score = kabsched_list1, current_score
                if not silent: print "    Info: Best RMSD so far with random 3-point Kabsch fitting: {0}".format(best_score)
            if show_graph: do_show_graph([(P-Pc,"P-Pc"), (Q-Qc, "Q-Qc"), (point_arrays[0] - Pc, "P1-Pc"), (point_arrays[1] - Qc, "P2-Qc")])
            
        if not silent: print "    Info: Returning best match with random 3-point Kabsch fitting (Score: {0})".format(best_score)
        return {'array': best_match.tolist(), 'score': best_score}
    else:
        assert map(on_second_element, unique_points[0][0:3]) == map(on_second_element, unique_points[1][0:3]), "Error: Unique points have not been ordered properly: {0} and {1}".format(map(on_second_element, unique_points[0][0:3]), map(on_second_element, unique_points[1][0:3]))
        
        # Align those three points using Kabsch algorithm
        P, Q = map(on_first_element, unique_points[0]), map(on_first_element, unique_points[1])
        U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)
        kabsched_list1 = np.dot(point_arrays[0]-Pc, U) + Qc

        if show_graph: do_show_graph([(kabsched_list1, "P1_kabsch"), (point_arrays[1], "P2")])

        current_match, current_score = kabsched_list1, distance_array_function(kabsched_list1, point_arrays[1], silent=silent)
        
        if not silent: print "    Info: Klabsch algorithm on unique element types found a better match with a Score of {0}".format(current_score)
        return {'array': current_match.tolist(), 'score': current_score}
    
#################
#### HELPERS ####
#################

def assert_array_equal(array1, array2, message="{0} and {1} are different"):
    assert np.allclose( array1, array2), message.format(array1, array2)

def assert_found_permutation(array1, array2, silent=True, hard_fail=False):
    perm_list = []
    for i, point1_array in enumerate(array1[:,0:3]):
        min_dist, min_index = None, None
        for j, point2_array in enumerate(array2[:,0:3]):
            distance = np.linalg.norm(point1_array-point2_array)
            min_dist = min(min_dist, distance) if min_dist else distance
            if distance == min_dist: min_index = j
        perm_list.append((i, min_index))

    offending_indexes = filter(lambda x: len(x[1])>=2, [ (value, list(group)) for value, group in itertools.groupby(perm_list, lambda x:x[1]) ])
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
    # Align those three points using Kabsch algorithm
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
