from numpy import sqrt, mean, square, minimum
import numpy as np
from scipy.spatial.distance import cdist
from Vector import Vector, rotmat, m2rotaxis
import math
import itertools
from functools import partial
import charnley_rmsd.kabsch as kabsch
from copy import deepcopy

on_self, on_first_element, on_second_element = lambda x:x, lambda x:x[0], lambda x:x[1]
on_second_element_and_flavour = lambda grouped_flavours, x: str(x[1]) + str(len(grouped_flavours[ x[2] ]))

N_COMPLEXITY = 3

ALLOW_SHORTCUTS = False

def alignPointsOnPoints(point_list1, point_list2, silent=True, use_AD=False, element_list1=None, element_list2=None, flavour_list1=None, flavour_list2=None, show_graph=False, rmsd_tolerance=1E-3):
    point_lists = [point_list1, point_list2]
    element_lists = [element_list1, element_list2]
    flavour_lists = [flavour_list1, flavour_list2]

    # Initilaizers
    distance_function = rmsd_array if not use_AD else ad_array
    has_flavours = True if flavour_list1 and flavour_list2 else False
    has_elements = True if element_list1 and element_list2 else False
    if not silent: print # Add a commencing newline

    # Assert that the fitting make sense
    assert len(point_lists[0]) == len(point_lists[1]), "Error: Size of point lists doesn't match: {0} and {1}".format(len(point_lists[0]), len(point_lists[1]))
    if has_flavours:
        assert len(flavour_lists[0]) == len(flavour_lists[1]), "Error: Size of flavour lists doesn't match: {0} and {1}".format(len(flavour_lists[0]), len(flavour_lists[1]))
        assert len(flavour_lists[0]) == len(point_lists[0]), "Error: Size of flavour lists doesn't match size of point lists: {0} and {1}".format(len(flavour_lists[0]), len(point_lists[1]))
        get_sorted_eqgroup_lengths = lambda flavour_list: sorted(map(len, group_by(flavour_list, on_self).values()))
        grouped_flavour_lists[0], grouped_flavour_lists[1] = group_by(flavour_lists[0], on_self), group_by(flavour_lists[1], on_self)
        assert get_sorted_eqgroup_lengths(flavour_lists[0]) == get_sorted_eqgroup_lengths(flavour_lists[1]), "Error: There is not a one to one mapping between the lengths of the flavour sets: {0} and {1}".format(get_sorted_eqgroup_lengths(flavour_lists[0]), get_sorted_eqgroup_lengths(flavour_lists[1]))
    if has_elements:
        assert len(element_list1) == len(element_list2), "Error: Size of element lists doesn't match: {0} and {1}".format(len(element_list1), len(element_list2))
        assert len(element_list1) == len(point_lists[0]), "Error: Size of element lists doesn't match size of point lists: {0} and {1}".format(len(element_list1), len(point_lists[1]))
        assert sorted(element_list1) == sorted(element_list2), "Error: There is not a one to one mapping of the element sets: {0} and {1}".format(sorted(element_list1), sorted(element_list2))

    point_arrays = map(np.array, point_lists)
    center_of_geometries = map(center_of_geometry, point_arrays)

    # First, remove tranlational part from both by putting the cog in (0,0,0)
    centered_point_arrays = [point_arrays[0] - center_of_geometries[0], point_arrays[1] - center_of_geometries[1]]

    # Assert than the center of geometry of the translated point list are now on (0,0,0)
    [ assert_array_equal(a, np.array([0,0,0])) for a in centered_point_arrays ]

    # Break now if there are no rotationnal component
    if rmsd(*centered_point_arrays) <= rmsd_tolerance and ALLOW_SHORTCUTS:
        return translated_point_array1 + cog2

    bruteforce_aligning_vectors_method()

    if has_elements: # Additional method if we have elements
        flavoured_kabsch_method()

    assert_found_permutation(best_aligned_point_array1, point_array2, silent=silent)
    return best_aligned_point_array1.tolist()

### METHODS ###

def bruteforce_aligning_vectors_method(centered_arrays):
    # First, select our first point on the translated structure; it is mandatory that this point is not on the center of geometry
    for point in centered_arrays[0][:,0:3]:
        reference_vectors[0] = Vector(point)
        break

    # Break now if there are no rotationnal component
    if minimum_rmsd <= rmsd_tolerance and ALLOW_SHORTCUTS:
        assert_found_permutation(best_aligned_point_array1, point_array2, silent=silent)
        if not silent: print "Info: A simple translation was enough to match the two set of points. Exiting successfully."
        return best_aligned_point_array1.tolist()

    # Then try all the rotation that put one on the atom of the first set into one of the atoms of the second sets
    # There are N such rotations
    best_aligned_array, best_rmsd = centered_arrays[0], distance_function(*centered_arrays, silent=silent)
    for i, point1_array in enumerate(centered_arrays[0][:,0:3]):
        for j, point2_array in enumerate(centered_arrays[1][:,0:3]):

            reference_vectors[1] = Vector(point2_array)

            r = rotmat(reversed(*reference_vectors))
            if not silent: print "Info: Rotation parameters: {0} deg, axis {1}".format(m2rotaxis(r)[0]*180/np.pi, m2rotaxis(r)[1])
            assert m2rotaxis(r)[0] != 180., "Error: 180 degree rotation matrix currently not working"

            rotated_point_array1 = np.dot(centered_arrays[0], r)

            # If the norm of the vector are the same, check that the rotation effectively put p on q
            if reference_vectors[1].norm() == reference_vectors[0].norm():
                assert_array_equal(rotated_point_array1[0, 0:3], reference_vectors[1]._ar)
            # Else do the same operation on the normalized vectors
            else:
                assert_array_equal(Vector(rotated_point_array1[0, 0:3]).normalized()._ar, reference_vectors[1].normalized()._ar)

            current_rmsd = distance_function(rotated_point_array1, centered_arrays[1], silent=silent)
            minimum_rmsd = minimum(minimum_rmsd, current_rmsd)
            if current_rmsd == minimum_rmsd: best_aligned_point_array1 = rotated_point_array1 + cog2

            if current_rmsd <= rmsd_tolerance:
                break
        # Only iterate over the first point of centered_arrays[0]
        break
    if not silent: print "Info: Minimum RMSD from unflavoured algorithm is: {0}".format(minimum_rmsd)

    # Break now if there if the unflavoured algorithm succeeded
    if minimum_rmsd <= rmsd_tolerance and ALLOW_SHORTCUTS:
        assert_found_permutation(best_aligned_point_array1, point_array2, silent=silent)
        if not silent: print "Info: Unflavoured algorithm succeeded in finding a match within the given tolerance. Exiting successfully."
        return best_aligned_point_array1.tolist()

def flavoured_kabsch_method():
    if not silent: print "Info: Found element types. Trying flavoured 3-point Kabsch algorithm on flavoured elements types ..."

    # Try to find three unique elements type points
    if has_flavours:
        element_points = map(lambda index: zip(point_lists[index], element_list[index], flavour_lists[index]), [0,1])
        grouping_functions = map( lambda index: partial(on_second_element_and_flavour, grouped_flavour_lists[0]) if has_flavours else on_second_element
    else:
        element_points = map(lambda index: zip(point_lists[index], element_list[index]), [0,1])
        grouping_functions = [on_second_element, on_second_element]

    grouped_element_points = map( lambda index:group_by(element_points[index], grouping_functions[index]), [0,1])
    unique_points = map(lambda grouped_element_point: [group[0] for group in grouped_element_point.values() if len(group) == 1], grouped_element_points)

    assert len(unique_points[0]) == len(unique_points[1]), "Error: Non matching number of unique points in {0} and {1}".format(*unique_points)

    if not silent: print "    Info: Unique groups found based on element types: {0}".format(unique_points1)

    if len(unique_points[0]) < 3:
        if not silent: print "    Warning: Unable to find at least three unique point with the elements provided. Trying to disambiguate enough points to make a fit."

        missing_points = 3 - len(unique_points[0])

        ambiguous_point_groups = map(lambda grouped_element_point: sorted([group for group in grouped_element_points1.values() if 1 < len(group) <= N_COMPLEXITY ], key=len), grouped_element_points]

        if len(ambiguous_point_groups[0]) <= missing_points:
            if not silent: print "Error: Couldn'd find enough point to disambiguate. Returning best found match ..."
            return None

        if not silent: print "    Info: Found enough point to disambiguate. Trying kabsch algorithm ..."

        unique_points[1] += map(on_first_element, ambiguous_point_groups[1])[0:missing_points]
        permutations_list = itertools.product(*map(range, [len(group) for group in ambiguous_point_groups[0][0:missing_points] ]))

        best_kabsched_list1, best_rmsd = None, None
        for permutation in permutations_list:

            ambiguous_unique_points1 = deepcopy(unique_points1)
            for i, ambiguous_group in enumerate(ambiguous_point_groups1):
                ambiguous_unique_points1.append(ambiguous_group[ permutation[i] ])
                if len(ambiguous_unique_points1) == 3: break
            # Align those three points using Kabsch algorithm
            print ambiguous_unique_points1
            print unique_points2
            P, Q = map(on_first_element, ambiguous_unique_points1), map(on_first_element, unique_points2)
            U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)
            kabsched_list1 = np.dot(point_array1-Pc, U) + Qc
            current_rmsd = distance_function(kabsched_list1, point_array2, silent=silent)
            if (not best_rmsd) or current_rmsd <= best_rmsd:
                best_kabsched_list1, best_rmsd = kabsched_list1, current_rmsd
                if not silent: print "    Info: Best RMSD so far with random 3-point Kabsch fitting: {0}".format(best_rmsd)
            if show_graph: do_show_graph([(P-Pc,"P-Pc"), (Q-Qc, "Q-Qc"), (point_array1-Pc,"P1-Pc"), (point_array2-Qc,"P2-Qc")])
        if not silent: print "    Info: Returning best match with random 3-point Kabsch fitting (RMSD: {0})".format(best_rmsd)
        return best_kabsched_list1.tolist()
    else:
        assert map(on_second_element, unique_points1[0:3]) == map(on_second_element, unique_points2[0:3]), "Error: Unique points have not been ordered properly: {0} and {1}".format(map(on_second_element, unique_points1[0:3]), map(on_second_element, unique_points2[0:3]))
        # Align those three points using Kabsch algorithm
        P, Q = map(on_first_element, unique_points1), map(on_first_element, unique_points2)
        U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)
        kabsched_list1 = np.dot(point_array1-Pc, U) + Qc

        if show_graph: do_show_graph([(point_array2, "P1"), (kabsched_list1, "P2")])

        current_match = kabsched_list1
        assert_array_equal( center_of_geometry(best_aligned_point_array1), cog2, "Error: Center of geometry of fitted list1 doesn't match center of geometry of list2 ({0} != {1})")
        current_rmsd = distance_function(current_match, point_array2, silent=silent)
        if current_rmsd <= minimum_rmsd:
            minimum_rmsd = current_rmsd
            best_aligned_point_array1 = current_match
            if not silent: print "    Info: Klabsch algorithm on unique element types found a better match with a RMSD of {0}".format(current_rmsd)
        else:
            if not silent: print "    Info: Klabsch algorithm on unique element types could not found a better match with (best found RMSD: {0})".format(current_rmsd)


#################
#### HELPERS ####
#################

def assert_array_equal(array1, array2, message="{0} and {1} are different"):
    assert np.allclose( array1, array2), message.format(array1, array2)

def assert_found_permutation(array1, array2, silent=True):
    perm_list = []
    for i, point1_array in enumerate(array1[:,0:3]):
        min_dist, min_index = None, None
        for j, point2_array in enumerate(array2[:,0:3]):
            distance = np.linalg.norm(point1_array-point2_array)
            min_dist = min(min_dist, distance) if min_dist else distance
            if distance == min_dist: min_index = j
        perm_list.append((i, min_index))

    offending_indexes = filter(lambda x: len(x[1])>=2, [ (value, list(group)) for value, group in itertools.groupby(perm_list, lambda x:x[1]) ])
    ambiguous_indexes = list( set(zip(*perm_list)[0]) - set(zip(*perm_list)[1]) ) + [value for value, group in offending_indexes]

    # Assert that perm_list is a permutation, i.e. that every obj of the first list is assigned one and only once to an object of the second list
    assert sorted(zip(*perm_list)[1]) == list(zip(*perm_list)[0]), "Error: {0} is not a permutation of {1}, which means that the best fit does not allow an unambiguous one-on-one mapping of the atoms. The method failed.".format(sorted(zip(*perm_list)[1]), zip(*perm_list)[0])
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

def rmsd(point_list1, point_list2):
    point_array1 = np.array(point_list1)
    point_array2 = np.array(point_list2)
    return rmsd_array(point_array1, point_array2)

def rmsd_array(point_array1, point_array2, silent=True):
    distance_matrix = get_distance_matrix(point_array1, point_array2)
    if not silent: print "    Info: Number of contact points: {0}/{1}".format(count_contact_points(distance_matrix), point_array1.shape[0])
    
    # Do you like my lisp skills?
    # This convoluted one-liner computes the square (R)oot of the (M)ean (S)quared (M)inimum (D)istances
    # We should call it the RMSMD :).
    # I think this is my favourite one-liner ever! (Probably because it look me 1 hour to construct and it's still beautiful)
    rmsd = sqrt( mean( square( np.min( distance_matrix, axis=0 ) ) ) )
    if not silent: print "    Info: New RMSD: {0}".format(rmsd)
    return rmsd

# Absolute Deviation
def ad(point_list1, point_list2):
    point_array1 = np.array(point_list1)
    point_array2 = np.array(point_list2)
    return ad_array(point_array1, point_array2)

def ad_array(point_array1, point_array2, silent=True):
    distance_matrix = get_distance_matrix(point_array1, point_array2)
    ad = max( np.min( distance_matrix, axis=0 ) )
    if not silent: print "    Info: New AD: {0}".format(ad)
    return ad

def distance(point1, point2):
    return np.linalg.norm(point1 - point2)

def center_of_geometry(point_array):
    return np.mean(point_array, axis=0)

def get_distance_matrix(x, y):
    return cdist(x, y, metric='euclidean')

# Return how many point of list1 are within 0.1 Angstrom to another point of list2
# Throws an error if several points are considered in contact to the same one
CONTACT_THRESHOLD = 0.15
def count_contact_points(distance_matrix):
    size = distance_matrix.shape[0]
    contacts = 0
    for line in distance_matrix[:,0:size]:
        new_contacts = sum([int(dist <= CONTACT_THRESHOLD) for dist in line])
        if new_contacts in [0,1]: contacts += new_contacts
        else: raise Exception("Error: Several points are in contact with the same one: {0}".format(line))
    return contacts

# A reimplemetation of python crappy itertool's broupby method with dictionnaries and less BS
def group_by(iterable, key):
    group_dict = {}
    for obj in iterable:
        group_dict.setdefault( key(obj), [])
        group_dict[key(obj)].append(obj)
    return group_dict
