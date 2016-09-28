from itertools import product

def N_amongst_array(point_array, N=3):
    def all_different(the_list):
        for a in the_list:
            if the_list.count(a) > 1 : return False
        return True
    N_points = point_array.shape[0]
    return [perm for perm in product(*[list(range(x)) for x in range(N_points, N_points - N, -1)]) if all_different(perm)] # Selecting three points at random amongst N
