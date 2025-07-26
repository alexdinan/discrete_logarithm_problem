
def dlp(p, g, h, n):
    # determines no. baby steps
    m = int((n / 2) ** 0.5) + 1

    # compute baby steps
    baby = {}
    val = 1
    for i in range(m):
        # store value, exponent pair in hash table
        baby[val] = i
        # update value
        val = (val * g) % p
    
    # precompute g^-m
    g_m_inv = pow(g, -m, p)

    # compute giant steps
    val = h
    for j in range((n // m) + 1):
        # check for collision
        if val in baby:
            return ((j * m) + baby[val]) % n
        val = (val * g_m_inv) % p

    # no solution found
    return None
