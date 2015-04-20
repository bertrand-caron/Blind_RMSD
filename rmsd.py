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

def alignPointsOnPoints(point_list1, point_list2, silent=True, use_AD=False, element_list1=None, element_list2=None, flavour_list1=None, flavour_list2=None, show_graph=False, rmsd_tolerance=1E-3):
    distance_function = rmsd_array if not use_AD else ad_array
    has_flavours = True if flavour_list1 and flavour_list2 else False
    has_elements = True if element_list1 and element_list2 else False
    assert len(point_list1) == len(point_list2), "Error: Size of point lists doesn't match: {0} and {1}".format(len(point_list1), len(point_list2))
    if has_flavours:
        assert len(flavour_list1) == len(flavour_list2), "Error: Size of flavour lists doesn't match: {0} and {1}".format(len(flavour_list1), len(flavour_list2))
        assert len(flavour_list1) == len(point_list1), "Error: Size of flavour lists doesn't match size of point lists: {0} and {1}".format(len(flavour_list1), len(point_list2))
        get_sorted_eqgroup_lengths = lambda flavour_list: sorted(map(len, group_by(flavour_list, on_self).values()))
        grouped_flavour_list1, grouped_flavour_list2 = group_by(flavour_list1, on_self), group_by(flavour_list2, on_self)
        assert get_sorted_eqgroup_lengths(flavour_list1) == get_sorted_eqgroup_lengths(flavour_list2), "Error: There is not a one to one mapping between the lengths of the flavour sets: {0} and {1}".format(get_sorted_eqgroup_lengths(flavour_list1), get_sorted_eqgroup_lengths(flavour_list2))
    if has_elements:
        assert len(element_list1) == len(element_list2), "Error: Size of element lists doesn't match: {0} and {1}".format(len(element_list1), len(element_list2))
        assert len(element_list1) == len(point_list1), "Error: Size of element lists doesn't match size of point lists: {0} and {1}".format(len(element_list1), len(point_list2))
        assert sorted(element_list1) == sorted(element_list2), "Error: There is not a one to one mapping of the element sets: {0} and {1}".format(sorted(element_list1), sorted(element_list2))
    if not silent: print # Add a commencing newline

    point_array1, point_array2 = map(np.array, (point_list1, point_list2))
    cog1, cog2 = map(center_of_geometry, (point_array1, point_array2))

    # First, remove tranlational part from both by putting the cog in (0,0,0)
    translated_point_array1, translated_point_array2 = point_array1 - cog1, point_array2 - cog2

    # Assert than the center of geometry of the translated point list are now on (0,0,0)
    assert_array_equal(center_of_geometry(translated_point_array1), np.array([0,0,0]))
    assert_array_equal(center_of_geometry(translated_point_array2), np.array([0,0,0]))

    # Break now if there are no rotationnal component
    if rmsd(translated_point_array1, translated_point_array2) <= rmsd_tolerance and ALLOW_SHORTCUTS: 
        return translated_point_array1 + cog2
    # First, select our first point on the translated structure; it is mandatory that this point is not on the center of geometry
    for point in translated_point_array1[:,0:3]:
        #print "point: {0}".format(point)
        pass

    point1_vector = Vector(translated_point_array1[0,0:3])

    # Assert than the first vector is effectively the first point of the point list translated by minus the center of mass
    assert_array_equal(point1_vector._ar, np.array(point_list1[0]) - cog1)

    minimum_rmsd = distance_function(translated_point_array1, translated_point_array2, silent=silent)
    best_aligned_point_array1 = translated_point_array1 + cog2

    # Break now if there are no rotationnal component
    if minimum_rmsd <= rmsd_tolerance and ALLOW_SHORTCUTS:
        assert_found_permutation(best_aligned_point_array1, point_array2, silent=silent)
        if not silent: print "Info: A simple translation was enough to match the two set of points. Exiting successfully."
        return best_aligned_point_array1.tolist()

    # Then try all the rotation that put one on the atom of the first set into one of the atoms of the second sets
    # There are N such rotations
    for i, point1_array in enumerate(translated_point_array1[:,0:3]):
        for j, point2_array in enumerate(translated_point_array2[:,0:3]):

            point2_vector = Vector(point2_array)
            #print "\nVector 2 is: {0}".format(point2_vector)

            # If the points are already superimposed, continue as the rotation matrix would be [[Nan, Nan, Nan], ...
            if point1_vector == point2_vector:
                if not silent: print "Info: {0} and {1} are superimposed".format(point1_vector, point2_vector)
                continue

            r = rotmat(point2_vector, point1_vector)
            if not silent: print "Info: Rotation parameters: {0} deg, axis {1}".format(m2rotaxis(r)[0]*180/np.pi, m2rotaxis(r)[1])
            assert m2rotaxis(r)[0] != 180., "Error: 180 degree rotation matrix currently not working"

            rotated_point_array1 = np.dot(translated_point_array1, r)

            # If the norm of the vector are the same, check that the rotation effectively put p on q
            if point2_vector.norm() == point1_vector.norm():
                assert_array_equal(rotated_point_array1[0, 0:3], point2_vector._ar)
            # Else do the same operation on the normalized vectors
            else:
                assert_array_equal(Vector(rotated_point_array1[0, 0:3]).normalized()._ar, point2_vector.normalized()._ar)

            current_rmsd = distance_function(rotated_point_array1, translated_point_array2, silent=silent)
            minimum_rmsd = minimum(minimum_rmsd, current_rmsd)
            if current_rmsd == minimum_rmsd: best_aligned_point_array1 = rotated_point_array1 + cog2

            if current_rmsd <= rmsd_tolerance:
                break
        # Only iterate over the first point of translated_point_array1
        break
    if not silent: print "Info: Minimum RMSD from unflavoured algorithm is: {0}".format(minimum_rmsd)

    # Break now if there if the unflavoured algorithm succeeded
    if minimum_rmsd <= rmsd_tolerance and ALLOW_SHORTCUTS:
        assert_found_permutation(best_aligned_point_array1, point_array2, silent=silent)
        if not silent: print "Info: Unflavoured algorithm succeeded in finding a match within the given tolerance. Exiting successfully."
        return best_aligned_point_array1.tolist()

    #has_flavours = False
    if has_elements: # Additional method if we have elements
        if not silent: print "Info: Found element types. Trying Kabsch algorithm on unique elements ..."
        # Try to find three unique elements type points
        element_points1, element_points2 = zip(point_list1, element_list1, flavour_list1) if has_flavours else zip(point_list1, element_list1), zip(point_list2, element_list2, flavour_list2) if has_flavours else zip(point_list2, element_list2)
        grouping_function1, grouping_function2 = partial(on_second_element_and_flavour, grouped_flavour_list1) if has_flavours else on_second_element, partial(on_second_element_and_flavour, grouped_flavour_list2) if has_flavours else on_second_element
        grouped_element_points1, grouped_element_points2 = group_by(element_points1, grouping_function1), group_by(element_points2, grouping_function2)
        unique_points1, unique_points2 = [group[0] for group in grouped_element_points1.values() if len(group) == 1], [group[0] for group in grouped_element_points2.values() if len(group) == 1 ]
        assert len(unique_points1) == len(unique_points2), "Error: Non matching number of unique points in {0} and {1}".format(unique_points1, unique_points2)
        if not silent: print "    Info: Unique groups found based on element types: {0}".format(unique_points1)
        if len(unique_points1) < 3:
            if not silent: print "    Warning: Unable to find at least three unique point with the elements provided. Trying to disambiguate enough points to make a fit."
            missing_points = 3 - len(unique_points1)
            ambiguous_point_groups1, ambiguous_point_groups2 = [group for group in grouped_element_points1.values() if 1 < len(group) <= N_COMPLEXITY ], [group for group in grouped_element_points2.values() if 1 < len(group) <= N_COMPLEXITY]
            ambiguous_point_groups1, ambiguous_point_groups2 = sorted(ambiguous_point_groups1, key=len), sorted(ambiguous_point_groups2, key=len)
            if len(ambiguous_point_groups1) <= missing_points:
                if not silent: print "Error: Couldn'd find enough point to disambiguate. Returning best found match ..."
                return best_aligned_point_array1.tolist()
            if not silent: print "    Info: Found enough point to disambiguate. Trying kabsch algorithm ..."
            unique_points2 = unique_points2 + map(on_first_element, ambiguous_point_groups2)[0:missing_points]
            permutations_list1 = itertools.product(*map(range, [len(group) for group in ambiguous_point_groups1[0:missing_points] ]))
            best_kabsched_list1, best_rmsd = None, None
            for permutation in permutations_list1:
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

    assert_found_permutation(best_aligned_point_array1, point_array2, silent=silent)
    return best_aligned_point_array1.tolist()

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
