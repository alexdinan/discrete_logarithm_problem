from math import sqrt, ceil, log, e
from random import randint


def dlp(p, g, h, n, c=1, step_num_const=2):
    # random map of elements to bound distances
    rand_bound_dist = lambda x, k: 2 ** (x % k)

    # define group order + width of search interval for x
    width = n - 1

    # k specifies the set of bound-sizes, S = {1,2,...,2^k-1}
    k = ceil(log(sqrt(width), 2)) + 3

    # how many tame kangaroo bounds 
    N = step_num_const * int(sqrt(width))

    # init. distinguished points
    f = max(0, int((0.5 * log(n, 2)) - (c * log(log(n, e), 2))))
    mask = (1 << f) - 1
    dist_pts = {}

    # compute tame kangaroo position after N steps
    s = randint(0, width)
    tame, total_tame_bound = pow(g, s, p), 0
    for _ in range(N):
        # check if curr position is a D.P
        if tame & mask == 0:
            dist_pts[tame] = total_tame_bound

        # next bound distance - based on curr position
        curr_tame_bound = rand_bound_dist(tame, k)
        total_tame_bound += curr_tame_bound
        # move kangaroo
        tame = (tame * pow(g, curr_tame_bound, p)) % p

    # add final position (trap) to set of D.Ps
    dist_pts[tame] = total_tame_bound
    
    # begin moving wild kangaroo
    wild, total_wild_bound = h, 0
    while total_wild_bound <= total_tame_bound + width:
        # check for collision with set of D.Ps
        if wild in dist_pts:
            return (dist_pts[wild] + s - total_wild_bound) % n
        
        # calculate next bound distance
        curr_wild_bound = rand_bound_dist(wild, k)
        total_wild_bound += curr_wild_bound
        # move wild kangaroo
        wild = (wild * pow(g, curr_wild_bound, p)) % p

    # no solution found
    return None
