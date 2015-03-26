from numpy import sqrt, mean, square, min
import numpy as np
from scipy.spatial.distance import cdist

def alignPointsOnPoints(point_list1, point_list2):
    point_array1 = np.array(point_list1)
    point_array2 = np.array(point_list2)
    cog1 = center_of_geometry(point_array1)
    cog2 = center_of_geometry(point_array2)
    translate_cog = cog2 - cog1
    aligned_point_array = point_array1 + translate_cog
    return aligned_point_array

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
    
