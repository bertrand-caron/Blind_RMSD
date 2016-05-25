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

def do_assert_equal(a, b, error_msg, exception_type=None, verbosity=0):
    if verbosity >= 3:
        print 'INFO: Asserting that {0} == {1}'.format(a, b)
    do_assert(
        (a == b if type(a) == float else is_close(a, b)),
        error_msg.format(a, b),
        exception_type=exception_type,
        verbosity=verbosity,
    )

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
            print success_msg
    except:
        print 'ERROR: Transformation was not isometric'
        if verbosity >= 5:
            print distance_matrix(array_1)
            print distance_matrix(array_2)

            for i, (d_1, d_2) in enumerate(zip(distance_matrix(array_1), distance_matrix(array_2))):
                if not is_close(d_1, d_2, atol=atol, rtol=rtol):
                    print i, d_1 - d_2

def assert_blind_rmsd_symmetry(array_1, array_2, distance_function, verbosity=1):
    do_assert_equal(
        distance_function(array_1, array_2),
        distance_function(array_2, array_1, transpose_mask_array=True),
        'ERROR: {0} != {1}',
        verbosity=verbosity,
    )

def assert_found_permutation_array(array1, array2, flavour_lists=None, mask_array=None, hard_fail=False, verbosity=0):
    from Blind_RMSD.align import FIRST_STRUCTURE, SECOND_STRUCTURE

    distance_matrix = get_distance_matrix(array1, array2)

    if mask_array is not None:
        distance_matrix += mask_array

    if verbosity >=5:
        print 'Distance matrix for permutation array:'
        print distance_matrix

    dim = distance_matrix.shape
    assert dim[0] == dim[1]

    perm_list = []
    for j in range(dim[0]):
        min_index = 0
        min_dist = distance_matrix[min_index, j]
        for i in range(1, dim[0]):
            distance = distance_matrix[i,j]
            if distance < min_dist:
                min_dist, min_index = distance, i
            elif distance == min_dist:
                if distance != INF:
                    raise Exception('ERROR: equal distance clash (i={0}, j={1})'.format(i, j))
        perm_list.append((j, min_index))

    offending_indexes = [
        (value, map(lambda x: x[0], group))
        for value, group in group_by(
            perm_list,
            lambda x: x[1],
        ).items()
        if len(group) >= 2
    ]

    misdefined_indexes = list( set(zip(*perm_list)[0]) - set(zip(*perm_list)[1]) ) + [value for value, group in offending_indexes]

    if flavour_lists is not None:
        if verbosity >= 4:
            print zip(
                misdefined_indexes,
                map(
                    lambda x: (flavour_lists[SECOND_STRUCTURE][x],),
                    misdefined_indexes,
                ),
            )

    # Assert that perm_list is a permutation, i.e. that every obj of the first list is assigned one and only once to an object of the second list
    if hard_fail:
        do_assert(
            sorted(zip(*perm_list)[1]) == list(zip(*perm_list)[0]),
            "ERROR: {0} is not a permutation of {1}, which means that the best fit does not allow an unambiguous one-on-one mapping of the atoms. The method failed.".format(
                sorted(zip(*perm_list)[1]),
                list(zip(*perm_list)[0]),
            ),
        )
        if verbosity >= 1:
            print "INFO: {0} is a permutation of {1}. This is a good indication the algorithm might have succeeded.".format(
                zip(*perm_list)[1],
                zip(*perm_list)[0],
            )
    else:
        if not sorted(zip(*perm_list)[1]) == list(zip(*perm_list)[0]):
            if (verbosity >= 1) or BYPASS_SILENT:
                print "ERROR: {0} is not a permutation of {1}, which means that the best fit does not allow an unambiguous one-on-one mapping of the atoms. The method failed.".format(
                    sorted(zip(*perm_list)[1]),
                    list(zip(*perm_list)[0]),
                )
                print "ERROR: Troublesome indexes are {0}".format(misdefined_indexes)
            final_permutation = None
        else:
            if (verbosity >= 1) or BYPASS_SILENT:
                print "INFO: {0} is a permutation of {1}. This is a good indication the algorithm might have succeeded.".format(
                    zip(*perm_list)[1],
                    zip(*perm_list)[0],
                )
            final_permutation = perm_list
    return final_permutation
