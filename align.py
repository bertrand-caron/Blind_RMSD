import numpy as np
np.set_printoptions(precision=1, linewidth=300)
from itertools import product, groupby, permutations
from functools import partial
from copy import deepcopy
from pprint import PrettyPrinter
from collections import namedtuple

from Blind_RMSD.helpers.Vector import Vector, rotmat, m2rotaxis
from Blind_RMSD.helpers.ChemicalPoint import ChemicalPoint, on_elements, on_coords, on_canonical_rep, ELEMENT_NUMBERS
from Blind_RMSD.helpers.moldata import group_by
from Blind_RMSD.helpers.permutations import N_amongst_array
from Blind_RMSD.helpers.scoring import rmsd_array, ad_array, rmsd, ad, rmsd_array_for_loop, get_distance_matrix

from Blind_RMSD.lib.charnley_rmsd import kabsch

pp = PrettyPrinter(indent=2)

on_self, on_first_element, on_second_element = lambda x:x, lambda x:x[0], lambda x:x[1]
on_third_element, on_fourth_element = lambda x: x[2], lambda x: x[3]
on_second_element_and_flavour = lambda grouped_flavours, x: str(x[1]) + str(len(grouped_flavours[ x[2] ]))

# Kabsch Algorithm options
DEFAULT_MIN_N_UNIQUE_POINTS = 5
MAX_N_COMPLEXITY = 10 # Maximum number of permutations is MAX_N_COMPLEXITY^(N_UNIQUE_POINTS - MIN_N_UNIQUE_POINTS)

ALLOW_SHORTCUTS = False
DEFAULT_SCORE_TOLERANCE = 0.01

DISABLE_BRUTEFORCE_METHOD = True

BYPASS_SILENT = False

ORIGIN, ZERO_VECTOR = np.array([0.,0.,0.]), np.array([0.,0.,0.])

DEFAULT_MATRICES = (np.identity(3), ZERO_VECTOR, ZERO_VECTOR)
NO_TRANSFORM = lambda point_array: point_array

FIRST_STRUCTURE, SECOND_STRUCTURE = (0, 1)
UNTIL_SECOND_STRUCTURE = SECOND_STRUCTURE + 1
EXTRA_POINTS = 2
ON_BOTH_LISTS = [FIRST_STRUCTURE, SECOND_STRUCTURE]

NULL_RMSD = 0.0
INFINITE_RMSD = float('inf')

Alignment = namedtuple('Alignment', 'aligned_points, score , extra_points, final_permutation')

FAILED_ALIGNMENT = Alignment(None, INFINITE_RMSD, None, None)

def do_assert(something, error_msg):
    assert something, error_msg

def transform_mapping(point_array_1, point_array_2):
    P, Q = point_array_1, point_array_2
    U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)
    transform = lambda point_array: np.dot(point_array - Pc, U) + Qc
    return transform

# Align points on points
def pointsOnPoints(point_lists, silent=True, use_AD=False, element_lists=None, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE, soft_fail=False, extra_points=[]):
    '''
    '''

    # Initializers
    has_elements = True if element_lists and all(element_lists) else False
    has_flavours = True if flavour_lists and all(flavour_lists) else False
    if not has_flavours:
        flavour_lists = map(
            lambda a_list: range(len(a_list)),
            element_lists,
        )
        has_flavours = True

    has_extra_points = (extra_points != [])

    # Assert that the fitting make sense
    do_assert(
        len(point_lists[FIRST_STRUCTURE]) == len(point_lists[SECOND_STRUCTURE]),
        "Error: Size of point lists doesn't match: {0} and {1}".format(*map(len, point_lists)),
    )
    if has_flavours:
        do_assert(
            len(flavour_lists[FIRST_STRUCTURE]) == len(flavour_lists[SECOND_STRUCTURE]),
            "Error: Size of flavour lists doesn't match: {0} and {1}".format(*map(len, flavour_lists)),
        )
        do_assert(
            len(flavour_lists[FIRST_STRUCTURE]) == len(point_lists[SECOND_STRUCTURE]),
            "Error: Size of flavour lists doesn't match size of point lists: {0} and {1}".format(*map(len, [flavour_lists[FIRST_STRUCTURE], point_lists[SECOND_STRUCTURE]])),
        )
        do_assert(
            sorted(flavour_lists[FIRST_STRUCTURE]) == sorted(flavour_lists[SECOND_STRUCTURE]),
            "Error: There is not a one to one mapping between the sorted elements of the flavour sets: {0} and {1}".format(*map(sorted, flavour_lists)),
        )
    if has_elements:
        do_assert(
            len(element_lists[FIRST_STRUCTURE]) == len(element_lists[SECOND_STRUCTURE]),
            "Error: Size of element lists doesn't match: {0} and {1}".format(*map(len, element_lists)),
        )
        do_assert(
            len(element_lists[FIRST_STRUCTURE]) == len(point_lists[SECOND_STRUCTURE]),
            "Error: Size of element lists doesn't match size of point lists: {0} and {1}".format(*map(len, [element_lists[FIRST_STRUCTURE], point_lists[SECOND_STRUCTURE]])),
        )
        do_assert(
            sorted(element_lists[FIRST_STRUCTURE]) == sorted(element_lists[SECOND_STRUCTURE]),
            "Error: There is not a one to one mapping of the element sets: {0} and {1}".format(*map(sorted, element_lists)),
        )

    point_arrays = map(
        np.array,
        point_lists + [extra_points],
    )
    center_of_geometries = map(
        center_of_geometry,
        point_arrays,
    )

    if has_elements:
        grouped_flavours_lists = map(
            lambda index: group_by(flavour_lists[index], lambda x:x ),
            ON_BOTH_LISTS,
        )
        chemical_points_lists = get_chemical_points_lists(
            point_lists,
            element_lists,
            flavour_lists,
            has_flavours,
            grouped_flavours_lists,
        )

        mask_array = np.zeros(
            (
                len(flavour_lists[FIRST_STRUCTURE]),
                len(flavour_lists[SECOND_STRUCTURE]),
            ),
        )
        dumb_array = np.chararray((len(flavour_lists[FIRST_STRUCTURE]), len(flavour_lists[FIRST_STRUCTURE])), itemsize=10)
        for i, chemical_point0 in enumerate(chemical_points_lists[FIRST_STRUCTURE]):
            for j, chemical_point1 in enumerate(chemical_points_lists[SECOND_STRUCTURE]):
                mask_array[i, j] = 0. if chemical_point0.canonical_rep == chemical_point1.canonical_rep else float('inf')
                dumb_array[i, j] = "{0} {1}".format(chemical_point0.canonical_rep, chemical_point1.canonical_rep)
        if not silent:
            print chemical_points_lists
        if not silent:
            print dumb_array

    distance_function, distance_array_function = rmsd if not use_AD else ad, lambda *args, **kwargs: rmsd_array_for_loop(*args, mask_array=mask_array if mask_array is not None else None, **kwargs) if not use_AD else ad_array

    # First, remove translational part from both by putting the center of geometry in (0,0,0)
    centered_point_arrays = [
        (a_point_array - a_center_of_geometry)
        for a_point_array, a_center_of_geometry in zip(
            point_arrays,
            center_of_geometries[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE] + [center_of_geometries[FIRST_STRUCTURE]],
        )
    ]

    # Assert than the center of geometry of the translated point list are now on (0,0,0)
    [ assert_array_equal( center_of_geometry(array), ORIGIN) for array in centered_point_arrays[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE] ]

    # Break now if the molecule has less than 3 atoms
    if len(point_lists[FIRST_STRUCTURE]) < 3 :
        return Alignment(
            (centered_point_arrays[FIRST_STRUCTURE]).tolist(),
            NULL_RMSD,
            (centered_point_arrays[EXTRA_POINTS].tolist() if has_extra_points else None),
            None,
        )

    # Break now if there are no rotational component
    if distance_function(*centered_point_arrays[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE]) <= score_tolerance and ALLOW_SHORTCUTS:
        if not silent: print "Info: A simple translation was enough to match the two set of points. Exiting successfully."
        assert_found_permutation_array(*centered_point_arrays[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE], silent=silent)

        def center_on_second_structure(point_array):
            return point_array + center_of_geometries[SECOND_STRUCTURE]

        return Alignment(
            center_on_second_structure(centered_point_arrays[FIRST_STRUCTURE]).tolist(),
            distance_function(*centered_point_arrays[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE]),
            center_on_second_structure(centered_point_arrays[EXTRA_POINTS]).tolist(),
            None,
        )

    method_results = {}
    if not DISABLE_BRUTEFORCE_METHOD:
        # Try the bruteforce method first
        method_results['bruteforce'] = bruteforce_aligning_vectors_method(
            point_arrays[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE],
            distance_array_function=distance_array_function,
            score_tolerance=score_tolerance,
            silent=silent and True,
        )
        method_results['lucky_kabsch'] = lucky_kabsch_method(
            point_lists,
            element_lists,
            flavour_lists=flavour_lists,
            distance_array_function=distance_array_function,
            score_tolerance=score_tolerance,
            show_graph=show_graph,
            silent=silent,
        )
        method_results['bruteforce_kabsch'] = bruteforce_kabsch_method(
            point_lists,
            element_lists,
            flavour_lists=flavour_lists,
            distance_array_function=distance_array_function,
            score_tolerance=score_tolerance,
            show_graph=show_graph,
            silent=silent,
        )

    # Try the flavoured Kabsch method if we have elements
    if has_elements:
        method_results['kabsch'] = flavoured_kabsch_method(
            point_lists,
            element_lists,
            flavour_lists=flavour_lists,
            distance_array_function=distance_array_function,
            score_tolerance=score_tolerance,
            show_graph=show_graph,
            silent=silent,
        )

    best_method = sorted(
        method_results.items(),
        key=lambda x:x[1]['score'] if 'score' in x[1] else 100.,
    )[0][0]
    best_match, best_score = method_results[best_method]['array'], method_results[best_method]['score']

    if has_extra_points:
        transform_function = method_results[best_method]['transform']
        aligned_extra_points_array = transform_function(centered_point_arrays[EXTRA_POINTS])
    else:
        aligned_extra_points_array = None

    if best_match == None:
        if not soft_fail:
            raise Exception("Best match is None. Something went wrong.")
        else:
            return FAILED_ALIGNMENT

    if not silent: print "Info: Scores of methods are: {0}".format(dict([ (k, v['score']) for (k,v) in method_results.items() if 'score' in v]))
    if not silent: print "Info: Best score was achieved with method: {0}".format(best_method)

    def corrected(points_array):
        return points_array - center_of_geometry(best_match) + center_of_geometries[SECOND_STRUCTURE]

    corrected_best_match = corrected(best_match)
    if has_extra_points:
        corrected_extra_points = corrected(aligned_extra_points_array)
    else:
        corrected_extra_points = None

    assert_array_equal(*map(center_of_geometry, [corrected_best_match, point_arrays[SECOND_STRUCTURE]]), message="{0} != {1}")

    def assert_found_permutation_array(array1, array2, mask_array=None, silent=True, hard_fail=False):
        perm_list = []
        masked_rmsd_array = get_distance_matrix(array1, array2) + mask_array
        dim = masked_rmsd_array.shape
        assert dim[0] == dim[1]
        for i in range(dim[0]):
            min_dist, min_index = None, None
            for j in range(dim[0]):
                distance = masked_rmsd_array[i,j]
                min_dist = min(min_dist, distance) if min_dist else distance
                if distance == min_dist: min_index = j
            perm_list.append((i, min_index))

        offending_indexes = [ (value, map(lambda x: x[0], group)) for value, group in group_by(perm_list, lambda x:x[1]).items() if len(group) >= 2]
        misdefined_indexes = list( set(zip(*perm_list)[0]) - set(zip(*perm_list)[1]) ) + [value for value, group in offending_indexes]

        if not silent:
            print zip(misdefined_indexes, map(lambda x: (element_lists[SECOND_STRUCTURE][x], flavour_lists[SECOND_STRUCTURE][x]), misdefined_indexes))

        # Assert that perm_list is a permutation, i.e. that every obj of the first list is assigned one and only once to an object of the second list
        if hard_fail:
            do_assert(
                sorted(zip(*perm_list)[1]) == list(zip(*perm_list)[0]),
                "Error: {0} is not a permutation of {1}, which means that the best fit does not allow an unambiguous one-on-one mapping of the atoms. The method failed.".format(sorted(zip(*perm_list)[1]), list(zip(*perm_list)[0])),
            )
            if not silent:
                print "Info: {0} is a permutation of {1}. This is a good indication the algorithm might have succeeded.".format(zip(*perm_list)[1], zip(*perm_list)[0])
        else:
            if not sorted(zip(*perm_list)[1]) == list(zip(*perm_list)[0]):
                if not silent or BYPASS_SILENT:
                    print "Error: {0} is not a permutation of {1}, which means that the best fit does not allow an unambiguous one-on-one mapping of the atoms. The method failed.".format(sorted(zip(*perm_list)[1]), list(zip(*perm_list)[0]))
                    print "Error: Troublesome indexes are {0}".format(misdefined_indexes)
                final_permutation = None
            else:
                if not silent or BYPASS_SILENT:
                    print "Info: {0} is a permutation of {1}. This is a good indication the algorithm might have succeeded.".format(zip(*perm_list)[1], zip(*perm_list)[0])
                final_permutation = perm_list
        return final_permutation

    final_permutation = assert_found_permutation_array(
        corrected_best_match,
        point_arrays[SECOND_STRUCTURE],
        mask_array=mask_array if mask_array is not None else None,
        silent=silent,
    )

    return Alignment(
        corrected_best_match.tolist(),
        method_results[best_method]['score'],
        corrected_extra_points.tolist() if has_extra_points else [],
        final_permutation,
    )

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
    for point_arrays[FIRST_STRUCTURE] in centered_arrays[0][:,0:3]:
        for point_arrays[SECOND_STRUCTURE] in centered_arrays[SECOND_STRUCTURE][:,0:3]:

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
                if not silent:
                    print "    Info: Found a really good match (Score={0}) worth aborting now. Exiting successfully.".format(best_score)
                break
        # Only iterate over the first point of centered_arrays[0]
        break

    if not silent:
        print "    Info: Minimum Score from bruteforce algorithm is: {0}".format(best_score)

    return {
        'array': best_match.tolist(),
        'score': best_score,
        'reference_array': centered_arrays[SECOND_STRUCTURE],
    }

def get_chemical_points_lists(point_lists, element_lists, flavour_lists, has_flavours, grouped_flavours_lists):
    N_points = len(point_lists[FIRST_STRUCTURE])
    chemical_points_lists = map(
        lambda index: [
            ChemicalPoint(*zipped_point, **{'grouped_flavours': grouped_flavours_lists[index]})
            for zipped_point in zip(
                point_lists[index],
                range(N_points),
                element_lists[index],
                flavour_lists[index] if has_flavours else [None] * N_points,
            )
        ],
        ON_BOTH_LISTS,
    )
    return chemical_points_lists

def flavoured_kabsch_method(point_lists, element_lists, silent=True, distance_array_function=rmsd_array, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE):
    point_arrays = map(
        np.array,
        point_lists,
    )
    has_flavours= bool(flavour_lists)
    MIN_N_UNIQUE_POINTS = DEFAULT_MIN_N_UNIQUE_POINTS if len(point_lists[0]) >= DEFAULT_MIN_N_UNIQUE_POINTS else 3

    if not silent:
        print "    Info: Found element types. Trying flavoured {0}-point Kabsch algorithm on flavoured elements types ...".format(MIN_N_UNIQUE_POINTS)

    grouped_flavours_lists = map(
        lambda index:group_by(flavour_lists[index], lambda x:x ),
        ON_BOTH_LISTS,
    )
    chemical_points_lists = get_chemical_points_lists(
        point_lists,
        element_lists,
        flavour_lists,
        has_flavours,
        grouped_flavours_lists,
    )

    grouped_chemical_points_lists = map(
        lambda chemical_points:group_by(chemical_points, on_canonical_rep),
        chemical_points_lists,
    )
    # Try to find MIN_N_UNIQUE_POINTS unique elements type points
    unique_points_lists = map(
        lambda grouped_chemical_points: [group[0] for group in grouped_chemical_points.values() if len(group) == 1],
        grouped_chemical_points_lists,
    )
    # Order them by decreasing element type
    unique_points_lists = map(
        lambda index: sorted(unique_points_lists[index], key=lambda x: ELEMENT_NUMBERS[ x.element.upper() ], reverse=True),
        ON_BOTH_LISTS,
    )

    if not len(unique_points_lists[FIRST_STRUCTURE]) == len(unique_points_lists[SECOND_STRUCTURE]):
        import yaml
        canonical_reps = [map(on_canonical_rep, unique_points_list) for unique_points_list in unique_points_lists]
        total_canonical_reps = set(canonical_reps[0] + canonical_reps[1])
        raise Exception( '''
Error: Non matching number of unique points in
{0}
(len={1})

and

{2}
(len={3})

Non-matching unique points are:
{4}

and

{5}
'''.format(
        yaml.dump(sorted(unique_points_lists[FIRST_STRUCTURE], key=lambda x: x.canonical_rep)),
        len(unique_points_lists[FIRST_STRUCTURE]),
        yaml.dump(sorted(unique_points_lists[SECOND_STRUCTURE], key=lambda x: x.canonical_rep)),
        len(unique_points_lists[SECOND_STRUCTURE]),
        '\n'.join([ str(chem_point) for chem_point in unique_points_lists[FIRST_STRUCTURE] if on_canonical_rep(chem_point) not in canonical_reps[1] ]),
        '\n'.join([ str(chem_point) for chem_point in unique_points_lists[SECOND_STRUCTURE] if on_canonical_rep(chem_point) not in canonical_reps[0] ]),
    ))

    if not silent:
        print "    Info: Unique groups found based on element types: {0}".format(unique_points_lists[0])

    if len(unique_points_lists[FIRST_STRUCTURE]) < MIN_N_UNIQUE_POINTS:
        if not silent:
            print "    Warning: Unable to find at least {N} unique point with the elements provided. Trying to disambiguate enough points to make a fit.".format(N=MIN_N_UNIQUE_POINTS)

        missing_points = MIN_N_UNIQUE_POINTS - len(unique_points_lists[FIRST_STRUCTURE])

        ambiguous_point_groups = map(
            lambda grouped_chemical_points: sorted(
                [group for group in grouped_chemical_points.values() if 1 < len(group) <= MAX_N_COMPLEXITY ],
                key=len,
            ),
            grouped_chemical_points_lists,
        )

        # Order them by number of heaviest atoms first
        ambiguous_point_groups = map(
            lambda index: sorted(
                ambiguous_point_groups[index],
                key=lambda x: ELEMENT_NUMBERS[ x[0].element.upper() ],
                reverse=True,
            ),
            ON_BOTH_LISTS,
        )

        N_ambiguous_points = sum(
            map(
                len,
                ambiguous_point_groups[FIRST_STRUCTURE],
            ),
        )

        if N_ambiguous_points < missing_points:
            if not silent: print "    Error: Couldn'd find enough point to disambiguate: {M} (unique points) + {P} (ambiguous points) < {N} (required points). Returning best found match ...".format(P=N_ambiguous_points, M=len(unique_points_lists[0]), N=MIN_N_UNIQUE_POINTS)
            return {
                'array': None,
                'score': None,
                'reference_array': point_arrays[SECOND_STRUCTURE],
                'transform': NO_TRANSFORM,
            }

        if not silent:
            print "    Info: Found enough point ({N}) to disambiguate. Trying kabsch algorithm ...".format(N=N_ambiguous_points)

        permutations_list = []
        atom_indexes = lambda chemical_points: map(lambda chemical_point: chemical_point.index, chemical_points)

        # For each ambiguous group
        ambiguous_points = 0
        N_list = []
        for group in ambiguous_point_groups[FIRST_STRUCTURE]:
            if ambiguous_points == missing_points:
                N_list.append( 0 )
            else:
                N_list.append( min(len(group), missing_points-ambiguous_points) )
                ambiguous_points += N_list[-1]

        if not silent: print "    Info: Ambiguous groups are: {0} (number of points taken in each group: {1})".format(ambiguous_point_groups[0], N_list)

        permutation_lists = map(
            lambda group, N: permutations(atom_indexes(group), r=N),
            ambiguous_point_groups[FIRST_STRUCTURE],
            N_list,
        )

        complete_permutation_list = product(*permutation_lists)
        #print list(complete_permutation_list)

        for group, N in zip(ambiguous_point_groups[SECOND_STRUCTURE], N_list):
            unique_points_lists[SECOND_STRUCTURE] += group[0:N]

        best_match, best_score, best_matrices = [None]*3
        for group_permutations in complete_permutation_list:
            ambiguous_unique_points = [deepcopy(unique_points_lists[FIRST_STRUCTURE])]
            for group in group_permutations:
                for index in group:
                    new_point = chemical_points_lists[0][index]
                    ambiguous_unique_points[FIRST_STRUCTURE].append(new_point)

            if not silent:
                print '        Info: Attempting a fit between points {0} and {1}'.format(ambiguous_unique_points[FIRST_STRUCTURE], unique_points_lists[SECOND_STRUCTURE])

            do_assert(
                map(on_elements, ambiguous_unique_points[0]) == map(on_elements, unique_points_lists[SECOND_STRUCTURE]),
                "Error: Trying to match points whose elements don't match: {0} != {1}".format(map(on_elements, ambiguous_unique_points[FIRST_STRUCTURE]), map(on_elements, unique_points_lists[SECOND_STRUCTURE])),
            )

            transform = transform_mapping(
                map(on_coords, ambiguous_unique_points[FIRST_STRUCTURE]),
                map(on_coords, unique_points_lists[SECOND_STRUCTURE]),
            )
            kabsched_list1 = transform(point_arrays[FIRST_STRUCTURE])
            current_score = distance_array_function(kabsched_list1, point_arrays[SECOND_STRUCTURE], silent=silent)

            if (not best_score) or current_score <= best_score:
                best_match, best_score, best_transform = kabsched_list1, current_score, transform

                if not silent:
                    print "    Info: Best score so far with random {0}-point Kabsch fitting: {1}".format(MIN_N_UNIQUE_POINTS, best_score)

                if current_score <= score_tolerance:
                    return {
                        'array': best_match.tolist(),
                        'score': best_score,
                        'reference_array': point_arrays[SECOND_STRUCTURE],
                        'transform': best_transform,
                    }

#            if show_graph:
#                do_show_graph(
#                    [
#                        (P-Pc,"P-Pc"),
#                        (Q-Qc, "Q-Qc"),
#                        (point_arrays[FIRST_STRUCTURE] - Pc, "P1-Pc"),
#                        (point_arrays[SECOND_STRUCTURE] - Qc, "P2-Qc"),
#                    ],
#                )

        if not silent:
            print "    Info: Returning best match with random {0}-point Kabsch fitting (Score: {1})".format(
                MIN_N_UNIQUE_POINTS,
                best_score,
            )

        return {
            'array': best_match.tolist(),
            'score': best_score,
            'reference_array': point_arrays[SECOND_STRUCTURE],
            'transform': best_transform,
        }
    else:
        do_assert(
            map(on_elements, unique_points_lists[FIRST_STRUCTURE][0:MIN_N_UNIQUE_POINTS]) == map(on_elements, unique_points_lists[SECOND_STRUCTURE][0:MIN_N_UNIQUE_POINTS]),
            "Error: Unique points have not been ordered properly: {0} and {1}".format(map(on_elements, unique_points_lists[0][0:MIN_N_UNIQUE_POINTS]), map(on_elements, unique_points_lists[1][0:MIN_N_UNIQUE_POINTS])),
        )

        # Align those MIN_N_UNIQUE_POINTS points using Kabsch algorithm
        current_transform = transform_mapping(
            map(on_coords, unique_points_lists[FIRST_STRUCTURE]),
            map(on_coords, unique_points_lists[SECOND_STRUCTURE]),
        )
        kabsched_list1 = current_transform(point_arrays[FIRST_STRUCTURE])

#        if show_graph:
#            do_show_graph(
#                [
#                    (kabsched_list1, "P1_kabsch"),
#                    (point_arrays[SECOND_STRUCTURE], "P2"),
#                ],
#            )

        current_match = kabsched_list1
        current_score= distance_array_function(
            kabsched_list1,
            point_arrays[SECOND_STRUCTURE],
            silent=silent,
        )

        if not silent:
            print "    Info: Klabsch algorithm on unique element types found a better match with a Score of {0}".format(current_score)

        return {
            'array': current_match.tolist(),
            'score': current_score,
            'reference_array': point_arrays[SECOND_STRUCTURE],
            'transform': current_transform,
        }

def lucky_kabsch_method(point_lists, element_lists, silent=True, distance_array_function=rmsd_array, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE):
    point_arrays = map(
        np.array,
        point_lists,
    )

    transform = transform_mapping(*point_arrays)
    kabsched_list1 = transform(point_arrays[FIRST_STRUCTURE])

    current_match = kabsched_list1
    current_score = distance_array_function(
        kabsched_list1,
        point_arrays[FIRST_STRUCTURE],
        silent=silent,
    )

    if not silent:
        print "    Info: Minimum Score from lucky Kabsch method is: {0}".format(current_score)

    return {
        'array': current_match.tolist(),
        'score': current_score,
    }

def bruteforce_kabsch_method(point_lists, element_lists, silent=True, distance_array_function=rmsd_array, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE):
    N_BRUTEFORCE_KABSCH = 4
    silent = True

    point_arrays = map(
        np.array,
        point_lists,
    )

    unique_points = [None, point_arrays[SECOND_STRUCTURE][0:N_BRUTEFORCE_KABSCH, 0:3] ]
    best_match, best_score = None, None

    for permutation in N_amongst_array(point_arrays[FIRST_STRUCTURE], N_BRUTEFORCE_KABSCH):
        unique_points[FIRST_STRUCTURE] = map(
            lambda index: point_arrays[FIRST_STRUCTURE][index, 0:3],
            permutation,
        )
        transform = transform_mapping(*unique_points)
        kabsched_list1 = transform(point_arrays[FIRST_STRUCTURE])

        current_match = kabsched_list1
        current_score = distance_array_function(
            kabsched_list1,
            point_arrays[SECOND_STRUCTURE],
            silent=silent,
        )

        if (not best_score) or current_score <= best_score:
            best_match, best_score = kabsched_list1, current_score
            if not silent:
                print "    Info: Best score so far with bruteforce {N}-point Kabsch fitting: {best_score}".format(
                    best_score=best_score,
                    N=N_BRUTEFORCE_KABSCH,
                )

    if not silent:
        print "    Info: Minimum Score from bruteforce Kabsch method is: {0}".format(best_score)

    return {
        'array': current_match.tolist(),
        'score': best_score,
    }

#################
#### HELPERS ####
#################

def assert_array_equal(array1, array2, message="{0} and {1} are different"):
    assert np.allclose( array1, array2, atol=1e-5), message.format(array1, array2)

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
