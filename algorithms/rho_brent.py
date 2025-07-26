from math import gcd


def dlp(p, g, h, n):
    # define random walk function
    rand_walk = [
        lambda x, a, b: (pow(x, 2, p), a * 2 % n, b * 2 % n),
        lambda x, a, b: (x * h % p, a, (b + 1) % n),
        lambda x, a, b: (x * g % p, (a + 1) % n, b)
    ]

    # set initial values
    x, a, b = 1, 0, 0
    x2, a2, b2 = rand_walk[x % 3](x, a, b)
    i, t = 1, 1

    while x != x2:
        if i == t:
            # teleport tortoise to hare
            x, a, b = x2, a2, b2
            t *= 2
        # move hare
        x2, a2, b2 = rand_walk[x2 % 3](x2, a2, b2)
        i += 1

    # collision found, solve
    sols = solve_cong((b2 - b) % n, (a - a2) % n, n)
    for x in sols:
        if pow(g, x, p) == h:
            return x
    # no solution found
    return None


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
