from Blind_RMSD.helpers.numpy_helpers import *
from Blind_RMSD.helpers.assertions import is_close

class Kabsch_Error(Exception):
    pass

def kabsch(P, Q):
    """
    Copied from https://github.com/charnley/rmsd/blob/master/rmsd/calculate_rmsd.py
    distributed under the BSD 2-Clause "Simplified" License (see LICENSE).

    Modified to raise a Kabsch_Error for linear or planar molecules.

    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)

    if sum([1 for x in S if is_close(x, 0.0)]) >= 1:
        raise Kabsch_Error("ERROR: Kabsch points are either coplanar or colinear. Algorithm won't work (SVD was: {0})".format(S))

    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U

def centroid(X):
    """
    Copied from https://github.com/charnley/rmsd/blob/master/rmsd/calculate_rmsd.py
    distributed under the BSD 2-Clause "Simplified" License (see LICENSE).

    Calculate the centroid from a vectorset X
    """
    return X.mean(axis=0)
