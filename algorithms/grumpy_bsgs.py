# import libraries
from math import gcd


def dlp(p, g, h, n):
    # determine step sizes
    m = int((n / 2) ** 0.5) + 1

    # initialisation
    baby, giant1, giant2 = {}, {}, {}

    baby_val = 1
    giant1_val, giant2_val = h, pow(h, 2, p)

    giant1_step = pow(g, m, p)
    giant2_step = pow(g, -(m + 1), p)

    # main loop
    for i in range(10000000):
        # store baby step
        baby[baby_val] = i
        # check for collisions
        if baby_val in giant1:
            return (i - (m * giant1[baby_val])) % n
        if baby_val in giant2:
            rhs = i + (giant2[baby_val] * (m + 1))
            return solve_cong(rhs, p, g, h, n)
        # compute next baby val
        baby_val = (baby_val * g) % p

        # store giant1 step
        giant1[giant1_val] = i
        # check for collisions
        if giant1_val in baby:
            return (baby[giant1_val] - (m * i)) % n
        if giant1_val in giant2:
            return ((m * i) + (giant2[giant1_val] * (m + 1))) % n
        # compute next giant1 val
        giant1_val = (giant1_val * giant1_step) % p

        # store giant2 step
        giant2[giant2_val] = i
        # check for collisions
        if giant2_val in baby:
            rhs = baby[giant2_val] + (i * (m + 1))
            return solve_cong(rhs, p, g, h, n)
        if giant2_val in giant1:
            return ((m * giant1[giant2_val]) + (i * (m + 1))) % n
        # compute next giant2 val
        giant2_val = (giant2_val * giant2_step) % p
    
    # no solution found
    return None


def solve_cong(rhs, p, g, h, order):
    """ Solves 2x ≡ i + j(m + 1) (mod order) """
    # if order is odd 2^-1 exists -> solve directly
    if gcd(2, order) == 1:
        return (rhs * pow(2, -1, order)) % order
    
    # reduce congruence
    rhs_red = rhs // 2
    order_red = order // 2

    # solve reduced congruence
    sol = rhs_red % order_red
    
    # verify sol.
    if h == pow(g, sol, p):
        return sol
    elif h == pow(g, sol + order_red, p):
        return sol + order_red
    else:
        return None
