from numpy import sqrt, mean, square, min, minimum
import numpy as np
from scipy.spatial.distance import cdist
from Vector import Vector, rotmat
import math

RMSD_TOLERANCE = 1E-3

def assert_array_equal(array1, array2):
    assert np.allclose( array1, array2), "{0} and {1} are different".format(array1, array2)


def alignPointsOnPoints(point_list1, point_list2):
    point_array1 = np.array(point_list1)
    point_array2 = np.array(point_list2)
    cog1 = center_of_geometry(point_array1)
    cog2 = center_of_geometry(point_array2)

    print "cog1: {0}".format(cog1)
    print "cog2: {0}".format(cog2)

    # First, remove tranlational part from both by putting the cog in (0,0,0)
    translated_point_array1 = point_array1 - cog1
    translated_point_array2 = point_array2 - cog2

    print "translated_point_array1:\n{0}\n".format(translated_point_array1)
    print "translated_point_array2:\n{0}\n".format(translated_point_array2)

    # Assert than the center of geometry of the translated point list are now on (0,0,0)
    assert_array_equal(center_of_geometry(translated_point_array1), np.array([0,0,0]))
    assert_array_equal(center_of_geometry(translated_point_array2), np.array([0,0,0]))

    # Break now if there are no rotationnal component
    if rmsd(translated_point_array1, translated_point_array2) <= RMSD_TOLERANCE: return translated_point_array1 + cog2

    # Them try all the rotation that put one on the atom of the first set into one of the atoms of the second sets
    # There are N such rotations

    # First, select our first point on the translated structure; it is mandatory that this point is not on the center of geometry
    for point in translated_point_array1[:,0:3]:
        print "point: {0}".format(point)

    point1_vector = Vector(translated_point_array1[0,0:3])

    # Assert than the first vector is effectively the first point of the point list translated by minus the center of mass
    assert_array_equal(point1_vector._ar, np.array(point_list1[0]) - cog1)

    minimum_rmsd = 100.
    best_aligned_point_array1 = translated_point_array1
    #print "Vector 1 is: {0}".format(point1_vector)

    for i, point2_array in enumerate(translated_point_array2[:,0:3]):

        point2_vector = Vector(point2_array)
        print "\nVector 2 is: {0}".format(point2_vector)

        # If the points are already superimposed, continue as the rotation matrix would be [[Nan, Nan, Nan], ...
        if point1_vector == point2_vector:
            print "{0} and {1} are superimposed".format(point1_vector, point2_vector)
            continue

        r = rotmat(point2_vector, point1_vector)
        print "\nMatrix rotating {0} onto {1}:".format(point1_vector, point2_vector)
        print r

        rotated_point_array1 = np.dot(translated_point_array1, r)

        print "\nCoordinate of first point array before rotation:"
        print translated_point_array1
        print "Coordinate after rotation:"
        print rotated_point_array1

        print "Coordinate after rotation in the other direction:"
        print np.dot(translated_point_array1, np.transpose(r))

        print "\nCoordinate of second point array:"
        print translated_point_array2

        # If the norm of the vector are the same, check that the rotation effectively put p on q
        if point2_vector.norm() == point1_vector.norm():
            assert_array_equal(rotated_point_array1[0, 0:3], point2_vector._ar)
        # Else do the same operation on the normalized vectors
        else:
            assert_array_equal(Vector(rotated_point_array1[0, 0:3]).normalized()._ar, point2_vector.normalized()._ar)

        current_rmsd = rmsd_array(rotated_point_array1, translated_point_array2)
        minimum_rmsd = minimum(minimum_rmsd, current_rmsd)
        if current_rmsd == minimum_rmsd: best_aligned_point_array1 = rotated_point_array1 + cog2

        print "    New RMSD after rotation: {0}".format(current_rmsd)

        print "    New AD after rotation: {0}".format(ad_array(rotated_point_array1, translated_point_array2))
        #print "{0} has been rotated to {1}".format(point1_vector, rotated_point_array1)

        if current_rmsd <= RMSD_TOLERANCE:
            break

    return best_aligned_point_array1

def rmsd(point_list1, point_list2):
    point_array1 = np.array(point_list1)
    point_array2 = np.array(point_list2)
    return rmsd_array(point_array1, point_array2)

def rmsd_array(point_array1, point_array2):
    distance_matrix = get_distance_matrix(point_array1, point_array2)
    
    # Do you like my lisp skills?
    # This convoluted one-liner computes the square (R)oot of the (M)ean (S)quared (M)inimum (D)istances
    # We should call it the RMSMD :).
    # I think this is my favourite one-liner ever! (Probably because it look me 1 hour to construct and it's still beautiful)
    rmsd = sqrt( mean( square( min( distance_matrix, axis=0 ) ) ) )
    return rmsd

# Absolute Deviation
def ad(point_list1, point_list2):
    point_array1 = np.array(point_list1)
    point_array2 = np.array(point_list2)
    return ad_array(point_array1, point_array2)

def ad_array(point_array1, point_array2):
    distance_matrix = get_distance_matrix(point_array1, point_array2)
    ad = max( min( distance_matrix, axis=0 ) )
    return ad

def distance(point1, point2):
    return np.linalg.norm(point1 - point2)

def center_of_geometry(point_array):
    return np.mean(point_array, axis=0)

def get_distance_matrix(x, y):
    return cdist(x, y, metric='euclidean')
    
