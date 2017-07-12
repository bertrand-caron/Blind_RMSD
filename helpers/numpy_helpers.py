import numpy as np
from numpy import sqrt, mean, square
from scipy.spatial.distance import cdist
np.seterr(all='raise')
np.set_printoptions(precision=3, linewidth=300)

def get_distance_matrix(x, y):
    return cdist(x, y, metric='euclidean')

Array = np.ndarray

array = np.array
