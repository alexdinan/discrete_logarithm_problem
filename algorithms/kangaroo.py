from math import log, ceil, sqrt
from random import randint


def dlp(p, g, h, n, step_num_const=2):
    # random map of elements to bound distances
    rand_bound_dist = lambda x, k: 2 ** (x % k)

    # define order / width of search interval
    width = n - 1

    # define the bound-distances {1, 2, ..., 2^k-1}
    k = ceil(log(sqrt(width), 2)) + 3

    # how many tame kangaroo bounds
    N = step_num_const * int(sqrt(width))

    # compute tame position after N steps
    s = randint(0, width)
    tame, tot_tame_bound = pow(g, s, p), 0
    for _ in range(N):
        # compute bound distance
        curr_tame_bound = rand_bound_dist(tame, k)
        tot_tame_bound += curr_tame_bound
        # update tame kangaroo position
        tame = (tame * pow(g, curr_tame_bound, p)) % p

    # move wild kangaroo
    wild, tot_wild_bound = h, 0
    while tot_wild_bound <= tot_tame_bound + width:
        # check for collision
        if wild == tame:
            return (tot_tame_bound + s - tot_wild_bound) % n
        
        # compute bound distance
        curr_wild_bound = rand_bound_dist(wild, k)
        tot_wild_bound += curr_wild_bound

        # update wild kangaroo position
        wild = (wild * pow(g, curr_wild_bound, p)) % p
    
    # No solution found
    return None
