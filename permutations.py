from itertools import product

def N_amongst_array(point_array, N=3):
    def all_different(a, b, c):
        return a != b and a != c and b != c
    N_points = point_array.shape[0]
    return [perm for perm in product(*map(lambda x:range(x), range(N_points, N_points - N, -1))) if all_different(*perm)] # Selecting three points at random amongst N
