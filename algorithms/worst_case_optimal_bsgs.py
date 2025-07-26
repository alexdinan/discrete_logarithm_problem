
def dlp(p, g, h, n):
    # determine the search space size
    m = int(n ** 0.5) + 1

    # compute baby steps - {1,g^1,...,g^n-1}
    baby = {}
    val = 1
    for i in range(m):
        # store pair (g^i (mod p) : i) in hash table
        baby[val] = i
        # compute g^i+1 (mod p)
        val = (val * g) % p
    
    # find g^-m (mod p)
    g_n_inv = pow(g, -m, p)
    
    # compute giant steps - {h,hg^-m,...,hg^-m^2}
    val = h
    for j in range(m):
        # check for collision
        if val in baby:
            return (j * m) + baby[val]
        val = (val * g_n_inv) % p
    
    # no solution found
    return None
