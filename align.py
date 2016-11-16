from Blind_RMSD.helpers.numpy_helpers import *
from itertools import product, groupby, permutations
from functools import partial
from copy import deepcopy
from pprint import PrettyPrinter
from collections import namedtuple

from Blind_RMSD.helpers.Vector import Vector, rotmat, m2rotaxis
from Blind_RMSD.helpers.ChemicalPoint import ChemicalPoint, on_coords, on_flavour, ELEMENT_NUMBERS
from Blind_RMSD.helpers.moldata import group_by
from Blind_RMSD.helpers.permutations import N_amongst_array
from Blind_RMSD.helpers.scoring import rmsd_array, ad_array, rmsd, ad, rmsd_array_for_loop, NULL_RMSD, INFINITE_RMSD
from Blind_RMSD.helpers.assertions import do_assert, assert_array_equal, assert_found_permutation_array, do_assert_is_isometry, distance_matrix, pdist, is_close, assert_blind_rmsd_symmetry
from Blind_RMSD.helpers.exceptions import Topology_Error
from Blind_RMSD.helpers.kabsch import kabsch, centroid, Kabsch_Error

pp = PrettyPrinter(indent=2).pprint

on_self, on_first_object, on_second_object = lambda x: x, lambda x: x[0], lambda x: x[1]
on_third_object, on_fourth_object = lambda x: x[2], lambda x: x[3]
on_second_object_and_flavour = lambda grouped_flavours, x: str(x[1]) + str(len(grouped_flavours[ x[2] ]))

# Kabsch Algorithm options
DEFAULT_MIN_N_UNIQUE_POINTS = 5
MAX_N_COMPLEXITY = 10 # Maximum number of permutations is MAX_N_COMPLEXITY^(N_UNIQUE_POINTS - MIN_N_UNIQUE_POINTS)

ALLOW_SHORTCUTS = False
DEFAULT_SCORE_TOLERANCE = 0.01

DISABLE_BRUTEFORCE_METHOD = True

ORIGIN, ZERO_VECTOR = np.array([0.,0.,0.]), np.array([0.,0.,0.])

DEFAULT_MATRICES = (np.identity(3), ZERO_VECTOR, ZERO_VECTOR)
NO_TRANSFORM = lambda point_array: point_array

FIRST_STRUCTURE, SECOND_STRUCTURE = (0, 1)
UNTIL_SECOND_STRUCTURE = SECOND_STRUCTURE + 1
EXTRA_POINTS = 2
ON_BOTH_LISTS = [FIRST_STRUCTURE, SECOND_STRUCTURE]

Alignment = namedtuple('Alignment', 'aligned_points, score , extra_points, final_permutation')

Alignment_Method_Result = namedtuple('Alignment_Method_Result', 'method_name, method_result')

FAILED_ALIGNMENT = Alignment(None, INFINITE_RMSD, None, None)

NULL_PDB_WRITING_FCT = lambda alignment, file_name: None

DUMMY_DUMP_PDB = lambda point_list, transform, file_name: None

def transform_mapping(P, Q, verbosity=0):
    assert len(P) == len(Q)

    old_P, old_Q = list(map(
        deepcopy,
        (P, Q),
    ))

    U, Pc, Qc = rotation_matrix_kabsch_on_points(P, Q)

    def transform(point_array):
        new_point_array = np.dot(point_array - Pc, U) + Qc

        if False:
            try:
                assert_array_equal(
                    center_of_geometry(point_array),
                    center_of_geometry(new_point_array),
                    'ERROR: {0} != {1}'.format(
                        center_of_geometry(point_array),
                        center_of_geometry(new_point_array),
                    ),
                    atol=0.01,
                    rtol=0.0,
                )
            except:
                if verbosity >= 5:
                    print(center_of_geometry(point_array))
                    print(center_of_geometry(new_point_array))
                    print(Pc)
                    print(Qc)
                raise

        return new_point_array

    # Make sure the arrays were not modified
    assert_array_equal(P, old_P)
    assert_array_equal(Q, old_Q)

    return transform

# Align points on points
def pointsOnPoints(point_lists, use_AD: bool = False, flavour_lists=None, show_graph: bool = False, score_tolerance=DEFAULT_SCORE_TOLERANCE, soft_fail: bool = False, extra_points=[], assert_is_isometry: bool = False, verbosity: int = 0, pdb_writing_fct=None):
    '''
    '''

    # Initializers
    has_flavours = True if flavour_lists and all(flavour_lists) else False

    if verbosity >= 1:
        print()

    if verbosity >= 3:
        print('INFO: {0}'.format(
            dict(has_flavours=has_flavours),
        ))

    if not has_flavours:
        flavour_lists = list(map(
            lambda a_list: list(range(len(a_list))),
            point_lists,
        ))
        has_flavours = True

    has_extra_points = (extra_points != [])

    # Assert that the fitting make sense
    do_assert(
        len(point_lists[FIRST_STRUCTURE]) == len(point_lists[SECOND_STRUCTURE]),
        "ERROR: Size of point lists doesn't match: {0} and {1}".format(*list(map(len, point_lists))),
        exception_type=Topology_Error,
    )
    if has_flavours:
        do_assert(
            len(flavour_lists[FIRST_STRUCTURE]) == len(flavour_lists[SECOND_STRUCTURE]),
            "ERROR: Size of flavour lists doesn't match: {0} and {1}".format(*list(map(len, flavour_lists))),
            exception_type=Topology_Error,
        )
        do_assert(
            len(flavour_lists[FIRST_STRUCTURE]) == len(point_lists[SECOND_STRUCTURE]),
            "ERROR: Size of flavour lists doesn't match size of point lists: {0} and {1}".format(*list(map(len, [flavour_lists[FIRST_STRUCTURE], point_lists[SECOND_STRUCTURE]]))),
        )
        do_assert(
            sorted(flavour_lists[FIRST_STRUCTURE]) == sorted(flavour_lists[SECOND_STRUCTURE]),
            "ERROR: There is not a one to one mapping between the sorted flavour of the sets: {0} and {1}".format(*list(map(sorted, flavour_lists))),
            exception_type=Topology_Error,
        )

    point_arrays = list(map(
        np.array,
        point_lists + [extra_points],
    ))
    center_of_geometries = list(map(
        center_of_geometry,
        point_arrays,
    ))

    if verbosity >= 4:
        print(point_arrays[0])
        print(point_arrays[1])
        print(point_arrays[2])

    if has_flavours:
        grouped_flavours_lists = list(map(
            lambda index: group_by(flavour_lists[index], lambda x:x ),
            ON_BOTH_LISTS,
        ))
        chemical_points_lists = get_chemical_points_lists(
            point_lists,
            flavour_lists,
            has_flavours,
        )

        mask_array = np.zeros(
            (
                len(flavour_lists[FIRST_STRUCTURE]),
                len(flavour_lists[SECOND_STRUCTURE]),
            ),
        )

        dumb_array = np.chararray(
            (len(flavour_lists[FIRST_STRUCTURE]), len(flavour_lists[FIRST_STRUCTURE])),
            itemsize=60,
        )

        for i, chemical_point_0 in enumerate(chemical_points_lists[FIRST_STRUCTURE]):
            for j, chemical_point_1 in enumerate(chemical_points_lists[SECOND_STRUCTURE]):
                mask_array[i, j] = 0. if chemical_point_0.flavour == chemical_point_1.flavour else float('inf')
                dumb_array[i, j] = "{0} =?= {1}".format(chemical_point_0.flavour, chemical_point_1.flavour)

        if verbosity >= 5:
            print('INFO: chemical_points_lists:')
            print(chemical_points_lists)
        if verbosity >= 5:
            print('INFO: dumb_array:')
            print(dumb_array)
    else:
        mask_array = None
        chemical_points_lists = None

    if not use_AD:
        distance_array_function = lambda array_1, array_2, transpose_mask_array=False: rmsd_array_for_loop(
            array_1,
            array_2,
            mask_array=(np.transpose(mask_array) if transpose_mask_array else mask_array) if mask_array is not None else None,
            verbosity=verbosity,
        )
    else:
        raise AssertionError('This has not been implemented yet.')

    if pdb_writing_fct:
        def dump_pdb(point_array, transform, file_name):
            pdb_writing_fct(
                Alignment(
                    point_array.tolist(),
                    distance_array_function(
                        point_array,
                        point_arrays[SECOND_STRUCTURE],
                    ),
                    transform(point_arrays[EXTRA_POINTS]).tolist() if has_extra_points else [],
                    None,
                ),
                file_name,
            )
    else:
        dump_pdb = DUMMY_DUMP_PDB

    dump_pdb(
        point_arrays[FIRST_STRUCTURE],
        NO_TRANSFORM,
        'init.pdb',
    )

    # First, remove translational part from both by putting the center of geometry in (0,0,0)
    centered_point_arrays = [
        (a_point_array - a_center_of_geometry)
        for a_point_array, a_center_of_geometry in zip(
            point_arrays,
            center_of_geometries[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE] + ([center_of_geometries[FIRST_STRUCTURE]] if has_extra_points else []),
        )
    ]

    # Assert than the center of geometry of the translated point list are now on (0,0,0)
    [ assert_array_equal( center_of_geometry(array), ORIGIN) for array in centered_point_arrays[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE] ]

    def center_on_second_structure(point_array):
        return point_array + center_of_geometries[SECOND_STRUCTURE]

    # Break now if the molecule has less than 3 atoms
    if len(point_lists[FIRST_STRUCTURE]) < 3 :
        return formatted_and_validated_Aligment(
            center_on_second_structure(centered_point_arrays[FIRST_STRUCTURE]),
            point_arrays[SECOND_STRUCTURE],
            distance_array_function,
            aligned_extra_points=centered_point_arrays[EXTRA_POINTS] if has_extra_points else None,
            chemical_points_lists=chemical_points_lists,
            mask_array=mask_array,
            verbosity=verbosity,
        )

    # Break now if there are no rotational component
    current_score = distance_array_function(*centered_point_arrays[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE])

    if verbosity >= 2:
        print('INFO: {0}'.format(dict(current_score=current_score)))

    if current_score <= score_tolerance and ALLOW_SHORTCUTS:
        if verbosity >= 1:
            print("INFO: A simple translation was enough to match the two set of points. Exiting successfully.")

        return formatted_and_validated_Aligment(
            center_on_second_structure(centered_point_arrays[FIRST_STRUCTURE]),
            point_arrays[SECOND_STRUCTURE],
            distance_array_function,
            aligned_extra_points=center_on_second_structure(centered_point_arrays[EXTRA_POINTS]),
            chemical_points_lists=chemical_points_lists,
            mask_array=mask_array,
            verbosity=verbosity,
        )

    method_results = {}

    def add_method_result(method_result):
        method_results[method_result.method_name] = method_result.method_result

    add_method_result(
        Alignment_Method_Result(
            'translation',
            {
                'array': point_arrays[FIRST_STRUCTURE],
                'score': current_score,
                'reference_array': centered_point_arrays[SECOND_STRUCTURE],
                'transform': NO_TRANSFORM,
            },
        ),
    )

    if not DISABLE_BRUTEFORCE_METHOD:
        add_method_result(
            bruteforce_aligning_vectors_method(
                point_arrays[FIRST_STRUCTURE:UNTIL_SECOND_STRUCTURE],
                distance_array_function,
                score_tolerance=score_tolerance,
                verbosity=verbosity,
            ),
        )

        add_method_result(
            lucky_kabsch_method(
                point_lists,
                distance_array_function,
                flavour_lists=flavour_lists,
                score_tolerance=score_tolerance,
                show_graph=show_graph,
                verbosity=verbosity,
            ),
        )

        add_method_result(
            bruteforce_kabsch_method(
                point_lists,
                distance_array_function,
                flavour_lists=flavour_lists,
                score_tolerance=score_tolerance,
                show_graph=show_graph,
                verbosity=verbosity,
            ),
        )

    # Try the flavoured Kabsch method if we have flavours
    if has_flavours:
        add_method_result(
            flavoured_kabsch_method(
                point_lists,
                distance_array_function,
                flavour_lists=flavour_lists,
                score_tolerance=score_tolerance,
                show_graph=show_graph,
                verbosity=verbosity,
                dump_pdb=dump_pdb,
            ),
        )

    best_method = sorted(
        list(method_results.items()),
        key=lambda x:x[1]['score'] if 'score' in x[1] else 100.,
    )[0][0]
    best_match, best_score = method_results[best_method]['array'], method_results[best_method]['score']

    if verbosity >= 1:
        print("INFO: Scores of methods are: {0}".format(
            dict([ (k, v['score']) for (k,v) in list(method_results.items()) if 'score' in v]),
        ))
    if verbosity >= 1:
        print("INFO: Best score was achieved with method: {0}".format(best_method))

    if has_extra_points:
        transform_function = method_results[best_method]['transform']
        aligned_extra_points_array = transform_function(point_arrays[EXTRA_POINTS])

        try:
            do_assert(
                is_close(
                    distance(center_of_geometries[FIRST_STRUCTURE], center_of_geometries[EXTRA_POINTS]),
                    distance(center_of_geometry(best_match), center_of_geometry(aligned_extra_points_array)),
                ),
                'ERROR: EXTRA POINTS WERE SHUFFLED AROUND {0} ({1}) {2} ({3})'.format(
                    center_of_geometries[FIRST_STRUCTURE] - center_of_geometries[EXTRA_POINTS],
                    distance(center_of_geometries[FIRST_STRUCTURE], center_of_geometries[EXTRA_POINTS]),
                    center_of_geometry(best_match) - center_of_geometry(aligned_extra_points_array),
                    distance(center_of_geometry(best_match), center_of_geometry(aligned_extra_points_array)),
                ),
            )
        except:
            if verbosity >= 3:
                print(center_of_geometries[FIRST_STRUCTURE], center_of_geometries[EXTRA_POINTS], distance(center_of_geometries[FIRST_STRUCTURE], center_of_geometries[EXTRA_POINTS]))
                print(center_of_geometry(best_match), center_of_geometry(aligned_extra_points_array), distance(center_of_geometry(best_match), center_of_geometry(aligned_extra_points_array)))
            raise
    else:
        aligned_extra_points_array = None

    if best_match is None:
        if not soft_fail:
            raise AssertionError("Best match is None. Something went wrong.")
        else:
            return FAILED_ALIGNMENT

    def corrected(points_array):
        return points_array - center_of_geometry(best_match) + center_of_geometries[SECOND_STRUCTURE]

    corrected_best_match = corrected(best_match)# + center_of_geometry(best_match)
    if has_extra_points:
        corrected_extra_points = corrected(aligned_extra_points_array)# + center_of_geometry(best_match)

        # Make sure the correction has not distorted the inter-distances between EXTRA_POINTS and FIRST_STRUCTURE
        if not soft_fail:
            if verbosity >= 5:
                print(extra_points)
                print()
                print(corrected_extra_points)
                print()
                print(corrected_extra_points - extra_points)
                print()
                print(center_of_geometry(best_match))

        complete_molecule_before = np.concatenate((point_arrays[FIRST_STRUCTURE], point_arrays[EXTRA_POINTS]))
        complete_molecule_after = np.concatenate((corrected_best_match, corrected_extra_points))

        if assert_is_isometry:

            do_assert_is_isometry(
                point_arrays[FIRST_STRUCTURE],
                corrected_best_match,
                verbosity=verbosity,
                success_msg='INFO: aligned_points is isometric',
            )

            do_assert_is_isometry(
                point_arrays[EXTRA_POINTS],
                corrected_extra_points,
                verbosity=verbosity,
                success_msg='INFO: extra_points is isometric',
            )

            do_assert_is_isometry(
                complete_molecule_after,
                complete_molecule_before,
                verbosity=verbosity,
                success_msg='INFO: complete_molecule is isometric',
            )
    else:
        corrected_extra_points = None

    return formatted_and_validated_Aligment(
        corrected_best_match,
        point_arrays[SECOND_STRUCTURE],
        distance_array_function,
        aligned_extra_points=corrected_extra_points,
        chemical_points_lists=chemical_points_lists,
        mask_array=mask_array,
        verbosity=verbosity,
        dump_pdb=dump_pdb,
    )

def formatted_and_validated_Aligment(aligned_point_array, reference_point_array, distance_array_function, chemical_points_lists=None, aligned_extra_points=None, mask_array=None, verbosity=0, dump_pdb=DUMMY_DUMP_PDB):
    assert_array_equal(*
        list(map(center_of_geometry, (aligned_point_array, reference_point_array,))),
        message="{0} != {1}"
    )

    dump_pdb(
        aligned_point_array,
        NO_TRANSFORM,
        'last.pdb'
    )

    final_permutation = assert_found_permutation_array(
        aligned_point_array,
        reference_point_array,
        chemical_points_lists=chemical_points_lists,
        mask_array=mask_array if mask_array is not None else None,
        verbosity=verbosity,
    )

    if final_permutation is not None:
        assert_blind_rmsd_symmetry(
            aligned_point_array,
            reference_point_array,
            distance_array_function,
            verbosity=verbosity,
        )

    return Alignment(
        aligned_point_array.tolist(),
        distance_array_function( #FIXME
            aligned_point_array,
            reference_point_array,
        ),
        aligned_extra_points.tolist() if aligned_extra_points is not None else [],
        final_permutation,
    )

### METHODS ###

def bruteforce_aligning_vectors_method(centered_arrays, distance_array_function, score_tolerance=DEFAULT_SCORE_TOLERANCE, verbosity=0):
    # First, select our first point on the translated structure; it is mandatory that this point is not on the center of geometry
    reference_vectors = [None, None]
    for point in centered_arrays[0][:,0:3]:
        reference_vectors[0] = Vector(point)
        break

    # Then try all the rotation that put one on the atom of the first set into one of the atoms of the second sets
    # There are N such rotations
    best_match, best_score = centered_arrays[0], distance_array_function(*centered_arrays)
    point_arrays = [None, None]
    for point_arrays[FIRST_STRUCTURE] in centered_arrays[0][:,0:3]:
        for point_arrays[SECOND_STRUCTURE] in centered_arrays[SECOND_STRUCTURE][:,0:3]:

            reference_vectors[1] = Vector(point_arrays[1])

            r = rotmat(*reversed(reference_vectors))

            if verbosity >= 4:
                print("    INFO: Rotation parameters: {0} deg, axis {1}".format(
                    m2rotaxis(r)[0] * 180 / np.pi,
                    m2rotaxis(r)[1],
                ))

            assert m2rotaxis(r)[0] != 180., "ERROR: 180 degree rotation matrix currently not working"

            rotated_point_arrays = [np.dot(centered_arrays[0], r)]

            # If the norm of the vector are the same, check that the rotation effectively put p on q
            if all([ vect.norm() != 0. for vect in reference_vectors]):
                if reference_vectors[1].norm() == reference_vectors[0].norm():
                    assert_array_equal(rotated_point_arrays[0][0, 0:3], reference_vectors[1]._ar)
                # Else do the same operation on the normalized vectors
                else:
                    assert_array_equal(Vector(rotated_point_arrays[0][0, 0:3]).normalized()._ar, reference_vectors[1].normalized()._ar)

            current_score = distance_array_function(
                rotated_point_arrays[0],
                centered_arrays[1],
            )

            if current_score <= best_score:
                best_match, best_score = rotated_point_arrays[0], current_score

            if best_score <= score_tolerance and ALLOW_SHORTCUTS:
                if verbosity >= 1:
                    print("    INFO: Found a really good match (Score={0}) worth aborting now. Exiting successfully.".format(
                        best_score,
                    ))
                break
        # Only iterate over the first point of centered_arrays[0]
        break

    if verbosity >= 2:
        print("    INFO: Minimum Score from bruteforce algorithm is: {0}".format(best_score))

    return Alignment_Method_Result(
        'bruteforce_aligning_vectors',
        {
            'array': best_match.tolist(),
            'score': best_score,
            'reference_array': centered_arrays[SECOND_STRUCTURE],
        },
    )

def get_chemical_points_lists(point_lists, flavour_lists, has_flavours):
    N_points = len(point_lists[FIRST_STRUCTURE])
    chemical_points_lists = list(map(
        lambda index: [
            ChemicalPoint(*zipped_point)
            for zipped_point in zip(
                point_lists[index],
                list(range(N_points)),
                flavour_lists[index] if has_flavours else [None] * N_points,
            )
        ],
        ON_BOTH_LISTS,
    ))
    return chemical_points_lists

def flavoured_kabsch_method(point_lists, distance_array_function, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE, extra_points=[], verbosity=0, dump_pdb=DUMMY_DUMP_PDB):
    point_arrays = list(map(
        np.array,
        point_lists,
    ))
    has_flavours= bool(flavour_lists)

    old_point_arrays = deepcopy(point_arrays)

    def assert_constant_point_arrays():
        for an_array, old_array in zip(point_arrays, old_point_arrays):
            assert_array_equal(
                an_array,
                old_array,
            )
        if verbosity >= 5:
            print('INFO: POINT ARRAYS WERE CONSTANT')

    if verbosity >= 5:
        print('    INFO: {0}'.format(dict(has_flavours=has_flavours)))

    MIN_N_UNIQUE_POINTS = min(
        DEFAULT_MIN_N_UNIQUE_POINTS,
        len(point_lists[0]),
    )

    if verbosity >= 3:
        print("    INFO: Found flavours. Trying flavoured {0}-point Kabsch algorithm on flavoured types ...".format(MIN_N_UNIQUE_POINTS))

    grouped_flavours_lists = list(map(
        lambda index: group_by(flavour_lists[index], lambda x:x ),
        ON_BOTH_LISTS,
    ))
    chemical_points_lists = get_chemical_points_lists(
        point_lists,
        flavour_lists,
        has_flavours,
    )

    grouped_chemical_points_lists = list(map(
        lambda chemical_points: group_by(chemical_points, on_flavour),
        chemical_points_lists,
    ))
    # Try to find MIN_N_UNIQUE_POINTS uniquely-flavoured type points
    unique_points_lists = list(map(
        lambda grouped_chemical_points: [group[0] for group in list(grouped_chemical_points.values()) if len(group) == 1],
        grouped_chemical_points_lists,
    ))
    # Order them by (decreasing flavour)
    unique_points_lists = list(map(
        lambda index: sorted(
            unique_points_lists[index],
            key=lambda x: (x.flavour,),
            reverse=True,
        ),
        ON_BOTH_LISTS,
    ))

    if not len(unique_points_lists[FIRST_STRUCTURE]) == len(unique_points_lists[SECOND_STRUCTURE]):
        import yaml
        flavours = [list(map(on_flavour, unique_points_list)) for unique_points_list in unique_points_lists]
        total_flavours = set(flavours[0] + flavours[1])
        raise AssertionError( '''
ERROR: Non matching number of unique points in
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
        yaml.dump(sorted(unique_points_lists[FIRST_STRUCTURE], key=lambda x: x.flavour)),
        len(unique_points_lists[FIRST_STRUCTURE]),
        yaml.dump(sorted(unique_points_lists[SECOND_STRUCTURE], key=lambda x: x.flavour)),
        len(unique_points_lists[SECOND_STRUCTURE]),
        '\n'.join([ str(chem_point) for chem_point in unique_points_lists[FIRST_STRUCTURE] if on_flavour(chem_point) not in flavours[1] ]),
        '\n'.join([ str(chem_point) for chem_point in unique_points_lists[SECOND_STRUCTURE] if on_flavour(chem_point) not in flavours[0] ]),
    ))

    if verbosity >= 3:
        print("    INFO: Unique groups found based on flavours: {0}".format(unique_points_lists[0]))

    if len(unique_points_lists[FIRST_STRUCTURE]) < MIN_N_UNIQUE_POINTS:
        if verbosity >= 3:
            print("    Warning: Unable to find at least {N} unique point with the flavoured points provided. Trying to disambiguate enough points to make a fit.".format(N=MIN_N_UNIQUE_POINTS))

        missing_points = MIN_N_UNIQUE_POINTS - len(unique_points_lists[FIRST_STRUCTURE])

        # Order groups by length, and then flavour
        ambiguous_point_groups = list(map(
            lambda grouped_chemical_points: sorted(
                [group for group in list(grouped_chemical_points.values()) if 1 < len(group) <= MAX_N_COMPLEXITY ],
                key= lambda group: (len(group), group[0].flavour),
            ),
            grouped_chemical_points_lists,
        ))

        N_ambiguous_points = sum(
            list(map(
                len,
                ambiguous_point_groups[FIRST_STRUCTURE],
            )),
        )

        if N_ambiguous_points < missing_points:
            if verbosity >= 1:
                print("    ERROR: Couldn'd find enough point to disambiguate: {M} (unique points) + {P} (ambiguous points) < {N} (required points). Returning best found match ...".format(
                    P=N_ambiguous_points,
                    M=len(unique_points_lists[0]),
                    N=MIN_N_UNIQUE_POINTS,
                ))

            return Alignment_Method_Result(
                'flavoured_kabsch_failed',
                {
                    'array': None,
                    'score': INFINITE_RMSD,
                    'reference_array': point_arrays[SECOND_STRUCTURE],
                    'transform': NO_TRANSFORM,
                },
            )

        if verbosity >= 2:
            print("    INFO: Found enough point ({N}) to disambiguate. Trying kabsch algorithm ...".format(N=N_ambiguous_points))

        permutations_list = []
        atom_indexes = lambda chemical_points: [chemical_point.index for chemical_point in chemical_points]

        # For each ambiguous group
        ambiguous_points = 0
        N_list = []
        for group in ambiguous_point_groups[FIRST_STRUCTURE]:
            if ambiguous_points == missing_points:
                N_list.append( 0 )
            else:
                N_list.append(
                    min(
                        len(group),
                        missing_points-ambiguous_points,
                    ),
                )
                ambiguous_points += N_list[-1]

        if verbosity >= 3:
            print("    INFO: Ambiguous groups are: {0} (number of points taken in each group: {1})".format(
                ambiguous_point_groups[0],
                N_list,
            ))

        permutation_lists = list(map(
            lambda group, N: permutations(atom_indexes(group), r=N),
            ambiguous_point_groups[FIRST_STRUCTURE],
            N_list,
        ))

        complete_permutation_list = product(*permutation_lists)

        for (group, N) in zip(ambiguous_point_groups[SECOND_STRUCTURE], N_list):
            unique_points_lists[SECOND_STRUCTURE] += group[0:N]

        best_match, best_score = [None]*2
        for (i, group_permutations) in enumerate(complete_permutation_list):
            ambiguous_unique_points = [deepcopy(unique_points_lists[FIRST_STRUCTURE])]

            for group in group_permutations:
                for index in group:
                    new_point = chemical_points_lists[0][index]
                    ambiguous_unique_points[FIRST_STRUCTURE].append(new_point)

            if verbosity >= 5:
                print('        INFO: Attempting a fit between points\n{0}\nand\n{1}'.format(
                    ambiguous_unique_points[FIRST_STRUCTURE],
                    unique_points_lists[SECOND_STRUCTURE],
                ))

            do_assert(
                list(map(on_flavour, ambiguous_unique_points[0])) == list(map(on_flavour, unique_points_lists[SECOND_STRUCTURE])),
                "ERROR: Trying to match points whose flavours don't match: {0} != {1}".format(
                    list(map(
                        on_flavour,
                        ambiguous_unique_points[FIRST_STRUCTURE],
                    )),
                    list(map(
                        on_flavour,
                        unique_points_lists[SECOND_STRUCTURE],
                    )),
                ),
            )

            try:
                transform = transform_mapping(
                    list(map(on_coords, ambiguous_unique_points[FIRST_STRUCTURE])),
                    list(map(on_coords, unique_points_lists[SECOND_STRUCTURE])),
                )
            except Kabsch_Error as e:
                if verbosity >= 1:
                    print(e)
                return Alignment_Method_Result('flavoured_kabsch', FAILED_ALIGNMENT)

            kabsched_list1 = transform(point_arrays[FIRST_STRUCTURE])

            current_score = distance_array_function(
                kabsched_list1,
                point_arrays[SECOND_STRUCTURE],
            )

            dump_pdb(
                kabsched_list1,
                transform,
                'kabsch_{0}.pdb'.format(i),
            )

            if (best_score is None) or current_score <= best_score:
                best_match, best_score, best_transform = kabsched_list1, current_score, transform

                if verbosity >= 5:
                    print("    INFO: Best score so far with random {0}-point Kabsch fitting: {1}".format(MIN_N_UNIQUE_POINTS, best_score))

                if current_score <= score_tolerance:
                    return Alignment_Method_Result(
                        'flavoured_kabsch_ambiguous_early_success',
                        {
                            'array': best_match.tolist(),
                            'score': best_score,
                            'reference_array': point_arrays[SECOND_STRUCTURE],
                            'transform': best_transform,
                        },
                    )

        if verbosity >= 5:
            print("    INFO: Returning best match with random {0}-point Kabsch fitting (Score: {1})".format(
                MIN_N_UNIQUE_POINTS,
                best_score,
            ))

        assert_constant_point_arrays()
        return Alignment_Method_Result(
            'flavoured_kabsch_ambiguous',
            {
                'array': best_match.tolist(),
                'score': best_score,
                'reference_array': point_arrays[SECOND_STRUCTURE],
                'transform': best_transform,
            },
        )
    else:
        do_assert(
            list(map(on_flavour, unique_points_lists[FIRST_STRUCTURE])) == list(map(on_flavour, unique_points_lists[SECOND_STRUCTURE])),
            "ERROR: Unique points have not been ordered properly: {0} and {1}".format(
                list(map(on_flavour, unique_points_lists[0])),
                list(map(on_flavour, unique_points_lists[1])),
            ),
        )

        do_assert(
            list(map(on_flavour, unique_points_lists[FIRST_STRUCTURE])) == list(map(on_flavour, unique_points_lists[SECOND_STRUCTURE])),
            'ERROR: Canonical representation of unique points do not match: {0} != {1}'.format(
                list(map(on_flavour, unique_points_lists[FIRST_STRUCTURE])),
                list(map(on_flavour, unique_points_lists[SECOND_STRUCTURE])),
            ),
        )

        # Align all unique_points using Kabsch algorithm
        try:
            current_transform = transform_mapping(
                list(map(on_coords, unique_points_lists[FIRST_STRUCTURE])),
                list(map(on_coords, unique_points_lists[SECOND_STRUCTURE])),
            )
        except Kabsch_Error as e:
            if verbosity >= 1:
                print(e)
            return Alignment_Method_Result('flavoured_kabsch', FAILED_ALIGNMENT)

        kabsched_list1 = current_transform(point_arrays[FIRST_STRUCTURE])

        if verbosity >= 5:
            print(kabsched_list1, point_arrays[SECOND_STRUCTURE])

        current_match = kabsched_list1

        current_score= distance_array_function(
            kabsched_list1,
            point_arrays[SECOND_STRUCTURE],
        )

        dump_pdb(
            kabsched_list1,
            current_transform,
            'kabsch_unique.pdb',
        )

        if verbosity >= 2:
            print("    INFO: Klabsch algorithm on unique flavours found a better match with a Score of {0}".format(
                current_score,
            ))

        assert_constant_point_arrays()

        return Alignment_Method_Result(
            'flavoured_kabsch_unique',
            {
                'array': current_match.tolist(),
                'score': current_score,
                'reference_array': point_arrays[SECOND_STRUCTURE],
                'transform': current_transform,
            }
        )

def lucky_kabsch_method(point_lists, distance_array_function, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE, verbosity=0):
    point_arrays = list(map(
        np.array,
        point_lists,
    ))

    try:
        transform = transform_mapping(*point_arrays)
    except Kabsch_Error as e:
        if verbosity >= 1:
            print(e)
        return Alignment_Method_Result('lucky_kabsch', FAILED_ALIGNMENT)

    kabsched_list1 = transform(point_arrays[FIRST_STRUCTURE])

    current_match = kabsched_list1
    current_score = distance_array_function(
        kabsched_list1,
        point_arrays[FIRST_STRUCTURE],
    )

    if verbosity >= 1:
        print("    INFO: Minimum Score from lucky Kabsch method is: {0}".format(
            current_score,
        ))

    return Alignment_Method_Result(
        'lucky_kabsch',
        {
            'array': current_match.tolist(),
            'score': current_score,
        },
    )

def bruteforce_kabsch_method(point_lists, distance_array_function, flavour_lists=None, show_graph=False, score_tolerance=DEFAULT_SCORE_TOLERANCE, verbosity=0):
    N_BRUTEFORCE_KABSCH = 4

    point_arrays = list(map(
        np.array,
        point_lists,
    ))

    unique_points = list(map(
        deepcopy,
        [None, point_arrays[SECOND_STRUCTURE][0:N_BRUTEFORCE_KABSCH, 0:3] ],
    ))

    best_match, best_score = None, None

    for permutation in N_amongst_array(point_arrays[FIRST_STRUCTURE], N_BRUTEFORCE_KABSCH):
        unique_points[FIRST_STRUCTURE] = list(map(
            lambda index: point_arrays[FIRST_STRUCTURE][index, 0:3],
            permutation,
        ))

        try:
            transform = transform_mapping(*unique_points)
        except Kabsch_Error as e:
            if verbosity >= 1:
                print(e)
            return Alignment_Method_Result('bruteforce_kabsch_method', FAILED_ALIGNMENT)

        kabsched_list1 = transform(point_arrays[FIRST_STRUCTURE])

        current_match = kabsched_list1

        current_score = distance_array_function(
            kabsched_list1,
            point_arrays[SECOND_STRUCTURE],
        )

        if (best_score is None) or current_score <= best_score:
            best_match, best_score = kabsched_list1, current_score
            if verbosity >= 3:
                print("    INFO: Best score so far with bruteforce {N}-point Kabsch fitting: {best_score}".format(
                    best_score=best_score,
                    N=N_BRUTEFORCE_KABSCH,
                ))

    if verbosity >= 3:
        print("    INFO: Minimum Score from bruteforce Kabsch method is: {0}".format(
            best_score,
        ))

    return Alignment_Method_Result(
        'bruteforce_kabsch',
        {
            'array': current_match.tolist(),
            'score': best_score,
        },
    )

#################
#### HELPERS ####
#################

def rotation_matrix_kabsch_on_points(points1, points2):
    # Align those points using Kabsch algorithm
    P, Q = np.array(points1), np.array(points2)
    #print P
    #print Q
    Pc, Qc = centroid(P), centroid(Q)
    P, Q = P - Pc, Q - Qc
    U = kabsch(P, Q)
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
