
def dlp(p, g, h, n):
    # determine search space size
    m = int(n ** 0.5) + 1

    # initialisation
    baby, giant = {}, {}
    g_m_inv = pow(g, m, p)
    baby_val, giant_val = 1, h

    for i in range(m):
        # store baby step
        baby[baby_val] = i
        # look for collision
        if baby_val in giant:
            return ((-giant[baby_val] * m) + i) % n
        # compute next baby step
        baby_val = (baby_val * g) % p

        # store giant step
        giant[giant_val] = i
        # look for collision
        if giant_val in baby:
            return ((-i * m) + baby[giant_val]) % n
        # compute next giant step
        giant_val = (giant_val * g_m_inv) % p
    
    # no solution found
    return None
