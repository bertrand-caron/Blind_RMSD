from Blind_RMSD.helpers.log import log
from Blind_RMSD.helpers.assertions import do_assert
from Blind_RMSD.helpers.numpy_helpers import *

NULL_RMSD = 0.0
INFINITE_RMSD = float('inf')

def rmsd(point_list1, point_list2, mask_array=None):
    point_array1 = np.array(point_list1)
    point_array2 = np.array(point_list2)
    return rmsd_array_for_loop(point_array1, point_array2, mask_array=mask_array)

def rmsd_array(point_array1, point_array2, mask_array = None, verbosity=0):
    do_assert(
        point_array1.shape == point_array2.shape,
        "Won't compute RMSD on arrays with different sizes: {0} and {1}".format(*[x.shape for x in [point_array1, point_array2]]),
    )

    distance_matrix = get_distance_matrix(point_array1, point_array2)

    if mask_array is not None:
        do_assert(
            mask_array.shape == distance_matrix.shape,
            'Shapes of mask arrays do not match: {0} != {1}'.format(
                mask_array.shape,
                distance_matrix.shape,
            ),
        )
    else:
        mask_array = np.zeros((point_array1.shape[0], point_array1.shape[0]))

    if verbosity >= 3:
        log.debug("Number of contact points: {0}/{1}".format(count_contact_points(distance_matrix), point_array1.shape[0]))

    rmsd = sqrt( mean( square( np.min( distance_matrix + np.transpose(mask_array), axis=0 ) ) ) )

    if verbosity >= 4:
        log.debug("New RMSD: {0}".format(rmsd))

    raise Exception() # Wrong until proven otherwise

    return rmsd

def rmsd_array_for_loop(point_array1, point_array2, mask_array = None, verbosity=0):
    assert point_array1.shape == point_array2.shape, "Error: Won't compute RMSD on arrays with different sizes: {0} and {1}".format(*[x.shape for x in [point_array1, point_array2]])

    distance_matrix = get_distance_matrix(point_array1, point_array2)

    if mask_array is not None:
        assert mask_array.shape == distance_matrix.shape, 'Error: Shapes of mask arrays do not match: {0} != {1}'.format(mask_array.shape, distance_matrix.shape)
    else:
        mask_array = np.zeros((point_array1.shape[0], point_array1.shape[0]))

    if verbosity >= 5:
        log.debug("Number of contact points: {0}/{1}".format(count_contact_points(distance_matrix), point_array1.shape[0]))
    #mask_array = np.transpose(mask_array)
    if verbosity >= 500:
        log.debug('mask_array:')
        log.debug(mask_array)
        log.debug('distance_matrix:')
        log.debug(distance_matrix)
    if verbosity >= 5:
        log.debug('masked distance_matrix:')
        log.debug(distance_matrix + mask_array)

    distances = []
    for point1 in range(mask_array.shape[0]):
        closest_distance = None
        for point2 in range(mask_array.shape[0]):
            current_distance = distance_matrix[point1, point2] + mask_array[point1, point2]
            if (closest_distance is None) or current_distance <= closest_distance:
                closest_distance = current_distance
        distances.append(closest_distance)
    rmsd = sqrt( mean( square( distances ) ) )

    assert rmsd != INFINITE_RMSD

    if verbosity >= 4:
        log.debug("New RMSD: {0}".format(rmsd))

    return rmsd


# Absolute Deviation
def ad(point_list1, point_list2, mask_array=None):
    point_array1 = np.array(point_list1)
    point_array2 = np.array(point_list2)
    return ad_array(point_array1, point_array2, mask_array=mask_array)

def ad_array(point_array1, point_array2, mask_array=None, verbosity=0):
    if mask_array:
        assert mask_array.shape == point_array1.shape, 'Error: Shapes of mask arrays do not match'
    else:
        mask_array = np.zeros((point_array1.shape[0], point_array1.shape[0]))
    assert point_array1.shape == point_array2.shape, "Error: Won't compute AD on arrays with diference sizes: {0} and {1}".format(*[x.shape for x in [point_array1, point_array2]])
    distance_matrix = get_distance_matrix(point_array1, point_array2)
    ad = max( np.min( distance_matrix, axis=0 ) )

    if verbosity >= 4:
        log.debug("New AD: {0}".format(ad))

    return ad

# Return how many point of list1 are within 0.1 Angstrom to another point of list2
# Throws an error if several points are considered in contact to the same one
CONTACT_THRESHOLD = 0.15
def count_contact_points(distance_matrix):
    size = distance_matrix.shape[0]
    contacts = 0
    for line in distance_matrix[:,0:size]:
        new_contacts = sum([int(dist <= CONTACT_THRESHOLD) for dist in line])
        if new_contacts in [0,1]:
            contacts += new_contacts
        else:
            raise Exception("Error: Several points are in contact with the same one: {0}".format(line))
    return contacts
