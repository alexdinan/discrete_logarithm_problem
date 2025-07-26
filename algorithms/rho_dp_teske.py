from math import gcd, log, sqrt, pi, e
from random import randint


def dlp(p, g, h, n, c=1.5, r=20):
    # r-adding iterating function
    rand_walk = []
    for _ in range(r):
        ms, ns = randint(1, n), randint(1, n)
        Ms = (pow(g, ms, p) * pow(h, ns, p)) % p
        rand_walk.append(make_walk_func(Ms, ms, ns, p, n))

    # distinguished points init.
    f = max(0, int((0.5 * log(n, 2)) - (c * log(log(n, e), 2))))
    dist_pts = {}
    mask = (1 << f) - 1
    i = 0
    max_iters = sqrt((pi * n) / 2)

    # main loop
    x, a, b = 1, 0, 0
    while True:
        # move forward in walk
        x, a, b = rand_walk[x % r](x, a, b)

        # check for collision with d.ps
        if x in dist_pts:
            break
        
        # is x a d.p
        if not (x & mask):
            dist_pts[x] = (a, b)
        
        # keep track of number of iterations
        i += 1
        if i > max_iters:
            # decrement f -> double |dist_pts|
            f = max(0, f - 1)
            mask = (1 << f) - 1
            # reset iteration counter
            i = 0

    # collision found - solve
    sols = solve_cong((b - dist_pts[x][1]) % n, (dist_pts[x][0] - a) % n, n)
    for x in sols:
        if pow(g, x, p) == h:
            return x
    # no solution found
    return None


def make_walk_func(Ms, ms, ns, p, n):
    return lambda x, a, b: (x * Ms % p, (a + ms) % n, (b + ns) % n)


def solve_cong(lhs, rhs, n):
    d = gcd(lhs, n)
    # no soutions exist
    if rhs % d != 0:
        return []
    
    # reduce congruence
    lhs_red, rhs_red, n_red = lhs // d, rhs // d, n // d

    # return list of possible solutions
    x0 = (rhs_red * pow(lhs_red, -1, n_red)) % n_red
    return [x0 + (k * n_red) for k in range(d)]
