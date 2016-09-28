from scipy.spatial.distance import cdist, pdist

from Blind_RMSD.helpers.numpy_helpers import *
from Blind_RMSD.helpers.moldata import group_by

BYPASS_SILENT = False

INF = float('inf')

def do_assert(something, error_msg, exception_type=None, verbosity=0):
    try:
        assert something, error_msg
    except AssertionError:
        if exception_type is None:
            raise
        else:
            raise exception_type(error_msg)

def do_assert_equal(a, b, error_msg, exception_type=None, soft_fail=False, verbosity=0):
    if verbosity >= 3:
        print('INFO: Asserting that {0} == {1}'.format(a, b))
    try:
        do_assert(
            (a == b if type(a) == float else is_close(a, b)),
            error_msg.format(a, b),
            exception_type=exception_type,
            verbosity=verbosity,
        )
    except:
        if not soft_fail:
            raise
        else:
            'ERROR: Failed assertion'

def assert_array_equal(array1, array2, message="\n{0}\n\nand\n\n{1}\nare different", rtol=1e-5, atol=1e-08):
    assert np.allclose( array1, array2, rtol=rtol, atol=atol), message.format(array1, array2)

def is_close(a, b, atol=1E-8, rtol=1E-5):
    return abs(a - b) <= (atol + rtol * max(a, b))

def distance_matrix(array):
    if False:
        return cdist(array, array)
    else:
        return pdist(array)

def do_assert_is_isometry(array_1, array_2, atol=1E-8, rtol=1E-5, success_msg='Success', verbosity=0):
    try:
        assert_array_equal(
            distance_matrix(array_1),
            distance_matrix(array_2),
            atol=atol, # Because of the rotation, the coordinates get truncated quite a bit
            rtol=rtol,
        )
        if verbosity >= 3:
            print(success_msg)
    except:
        print('ERROR: Transformation was not isometric')
        if verbosity >= 5:
            print(distance_matrix(array_1))
            print(distance_matrix(array_2))

            for i, (d_1, d_2) in enumerate(zip(distance_matrix(array_1), distance_matrix(array_2))):
                if not is_close(d_1, d_2, atol=atol, rtol=rtol):
                    print(i, d_1 - d_2)

def assert_blind_rmsd_symmetry(array_1, array_2, distance_function, verbosity=1):
    do_assert_equal(
        distance_function(array_1, array_2),
        distance_function(array_2, array_1, transpose_mask_array=True),
        'ERROR: {0} != {1}',
        verbosity=verbosity,
        soft_fail=True,
    )

def assert_found_permutation_array(array1, array2, chemical_points_lists=None, mask_array=None, hard_fail=False, verbosity=0):
    from Blind_RMSD.align import FIRST_STRUCTURE, SECOND_STRUCTURE

    distance_matrix = get_distance_matrix(array1, array2)

    if mask_array is not None:
        distance_matrix += mask_array

    if verbosity >=5:
        print('Distance matrix for permutation array:')
        print(distance_matrix)

    dim = distance_matrix.shape
    assert dim[0] == dim[1]

    def closest_point_in_second_structure_to(i):
        min_j = 0
        min_dist = distance_matrix[i, min_j]
        for j in range(1, dim[1]):
            distance = distance_matrix[i, j]
            if distance < min_dist:
                min_dist, min_j = distance, j
            elif distance == min_dist:
                if False:
                    if distance != INF:
                        raise Exception('ERROR: equal distance clash (i={0}, j={1})'.format(i, j))
        return (
            chemical_points_lists[FIRST_STRUCTURE][i],
            chemical_points_lists[SECOND_STRUCTURE][min_j],
            min_dist,
        )

    point_mapping = [closest_point_in_second_structure_to(i) for i in range(dim[0])]

    def mapped_with_condition(condition = lambda group: len(group) == 1):
        return set(
            reduce(
                lambda acc, e: acc + e,
                [
                    group
                    for (key, group) in list(group_by(
                        [point_2 for (_, point_2, _) in point_mapping],
                        key=lambda x: x,
                    ).items())
                    if condition(group)
                ],
                [],
            ),
        )

    mapped_several_times = mapped_with_condition(condition=lambda group: (len(group) >= 2))
    mapped_exactly_once = mapped_with_condition(condition=lambda group: (len(group) == 1))

    not_mapped = set(chemical_points_lists[SECOND_STRUCTURE]) - mapped_several_times - mapped_exactly_once

    print('ERROR: Points of reference structure mapped several times: {0}'.format(
        mapped_several_times,
    ))

    for point in mapped_several_times:
        print('    ERROR: Point {0} from reference structure was mapped to the following points of the aligned structure: {1}'.format(
            point,
            [point_1 for (point_1, point_2, _) in point_mapping if point_2 == point],
        ))

    print('ERROR: Points of reference structure not mapped at all: {0}'.format(
        not_mapped,
    ))

    try:
        assert mapped_several_times == set()
    except AssertionError:
        if hard_fail:
            raise

    index_permutation = [(point_1.index, point_2.index) for (point_1, point_2, _) in point_mapping]

    if verbosity >= 3:
        print('INFO: Found a permutation between points of the aligned structure and points of the reference structure: {0} and {1}'.format(*list(zip(*index_permutation))))

    return index_permutation
