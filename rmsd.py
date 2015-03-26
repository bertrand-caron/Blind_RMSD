from numpy import sqrt, mean, square, min, minimum
import numpy as np
from scipy.spatial.distance import cdist
from Vector import Vector, rotmat
import math

RMSD_TOLERANCE = 1E-3

def alignPointsOnPoints(point_list1, point_list2):
    point_array1 = np.array(point_list1)
    point_array2 = np.array(point_list2)
    cog1 = center_of_geometry(point_array1)
    cog2 = center_of_geometry(point_array2)

    # First, remove tranlational part
    translate_cog = cog2 - cog1
    print "Translational component: {0}".format(translate_cog)
    translated_point_array1 = point_array1 + translate_cog

    # Break now if there are no rotationnal component
    if rmsd(translated_point_array1, point_array2) <= RMSD_TOLERANCE: return translated_point_array1

    # Them try all the rotation that put one on the atom of the first set into one of the atoms of the second sets
    # There are N such rotations

    # First, select our first point on the translated structure; it is mandatory that this point is not on the center of geometry
    for poin1 in point_list1[0]:
        pass

    point1_vector = Vector(translated_point_array1[0])
    minimum_rmsd = 100.
    best_aligned_point_array1 = translated_point_array1
    #print "Vector 1 is: {0}".format(point1_vector)

    for point2 in point_list2:

        point2_vector = Vector(point2)
        #print "\nVector 2 is: {0}".format(point2_vector)

        # If the points are already superimposed, continue as the rotation matrix would be [[Nan, Nan, Nan], ...
        if point1_vector == point2_vector:
            print "{0} and {1} are superimposed".format(point1_vector, point2_vector)
            continue

        r = rotmat(point2_vector, point1_vector)
        print "\nMatrix rotating {0} onto {1}:".format(point1_vector, point2_vector)
        print r
        rotated_point1_vector = point1_vector.left_multiply(r)
        rotated_point_array1 = np.dot(translated_point_array1, r)
        print "Coordinate before rotation:"
        print translated_point_array1
        print "Coordinate after rotation:"
        print rotated_point_array1

        #print rotated_point_array1
        #print point_array2

        current_rmsd = rmsd_array(rotated_point_array1, point_array2)
        minimum_rmsd = minimum(minimum_rmsd, current_rmsd)
        if current_rmsd == minimum_rmsd: best_aligned_point_array1 = rotated_point_array1

        print "    New RMSD after rotation: {0}".format(current_rmsd)
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

def distance(point1, point2):
    return np.linalg.norm(point1 - point2)

def center_of_geometry(point_array):
    return np.mean(point_array, axis=0)

def get_distance_matrix(x, y):
    return cdist(x, y, metric='euclidean')
    
