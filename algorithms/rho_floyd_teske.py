from random import randint
from math import gcd


def dlp(p, g, h, n, r=20):
    # r-adding iterating function
    rand_walk = []
    for _ in range(r):
        ms, ns = randint(1, n), randint(1, n)
        Ms = (pow(g, ms, p) * pow(h, ns, p)) % p
        rand_walk.append(make_walk_func(Ms, ms, ns, p, n))

    # init. hare/tortoise positions
    x, a, b = 1, 0, 0
    x2, a2, b2 = rand_walk[x % r](x, a, b)

    while x != x2:
        # move tortoise
        x, a, b = rand_walk[x % r](x, a, b)

        # move hare
        x2, a2, b2 = rand_walk[x2 % r](x2, a2, b2)
        x2, a2, b2 = rand_walk[x2 % r](x2, a2, b2)

    # collision found - solve
    sols = solve_cong((b2 - b) % n, (a - a2) % n, n)
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
