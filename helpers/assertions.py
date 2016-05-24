from scipy.spatial.distance import cdist, pdist

from Blind_RMSD.helpers.numpy_helpers import *
from Blind_RMSD.helpers.moldata import group_by

BYPASS_SILENT = False

def do_assert(something, error_msg, exception_type=None):
    try:
        assert something, error_msg
    except:
        if exception_type is None:
            raise
        else:
            raise exception_type(error_msg)

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

def assert_found_permutation_array(array1, array2, element_lists=None, flavour_lists=None, mask_array=None, hard_fail=False, verbosity=0):
    from Blind_RMSD.align import FIRST_STRUCTURE, SECOND_STRUCTURE

    perm_list = []
    masked_rmsd_array = get_distance_matrix(array1, array2) + mask_array
    dim = masked_rmsd_array.shape
    assert dim[0] == dim[1]
    for i in range(dim[0]):
        min_dist, min_index = None, None
        for j in range(dim[0]):
            distance = masked_rmsd_array[i,j]
            min_dist = min(min_dist, distance) if min_dist is not None else distance
            if distance == min_dist:
                min_index = j
        perm_list.append((i, min_index))

    offending_indexes = [
        (value, map(lambda x: x[0], group))
        for value, group in group_by(
            perm_list,
            lambda x:x[1],
        ).items()
        if len(group) >= 2
    ]

    misdefined_indexes = list( set(zip(*perm_list)[0]) - set(zip(*perm_list)[1]) ) + [value for value, group in offending_indexes]

    if verbosity >= 4:
        print zip(
            misdefined_indexes,
            map(
                lambda x: (element_lists[SECOND_STRUCTURE][x], flavour_lists[SECOND_STRUCTURE][x]),
                misdefined_indexes,
            ),
        )

    # Assert that perm_list is a permutation, i.e. that every obj of the first list is assigned one and only once to an object of the second list
    if hard_fail:
        do_assert(
            sorted(zip(*perm_list)[1]) == list(zip(*perm_list)[0]),
            "Error: {0} is not a permutation of {1}, which means that the best fit does not allow an unambiguous one-on-one mapping of the atoms. The method failed.".format(sorted(zip(*perm_list)[1]), list(zip(*perm_list)[0])),
        )
        if verbosity >= 1:
            print "INFO: {0} is a permutation of {1}. This is a good indication the algorithm might have succeeded.".format(zip(*perm_list)[1], zip(*perm_list)[0])
    else:
        if not sorted(zip(*perm_list)[1]) == list(zip(*perm_list)[0]):
            if (verbosity >= 1) or BYPASS_SILENT:
                print "ERROR: {0} is not a permutation of {1}, which means that the best fit does not allow an unambiguous one-on-one mapping of the atoms. The method failed.".format(sorted(zip(*perm_list)[1]), list(zip(*perm_list)[0]))
                print "ERROR: Troublesome indexes are {0}".format(misdefined_indexes)
            final_permutation = None
        else:
            if (verbosity >= 1) or BYPASS_SILENT:
                print "INFO: {0} is a permutation of {1}. This is a good indication the algorithm might have succeeded.".format(zip(*perm_list)[1], zip(*perm_list)[0])
            final_permutation = perm_list
    return final_permutation
