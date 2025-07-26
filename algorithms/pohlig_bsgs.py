from sympy.ntheory import factorint


def dlp(p, g, h, n, prime_factors=None):
    if prime_factors is None:
        prime_factors = factorint(n)

    moduli, remainders = [], []
    # solve smaller DLP instances independently
    for prime, power in prime_factors.items():
        # solve using two-stage pohlig-hellman
        sol = pohlig_second_stage(p, g, h, n, prime, power)
        moduli.append(pow(prime, power))
        remainders.append(sol)
    
    # Combine solutions using Chinese Remainder Theorem
    return crt(moduli, remainders)


def pohlig_second_stage(p, g, h, n, prime, power):
    # initialise variables
    x = 0
    g_mod = pow(g, n // prime, p)
    h_i = h
    curr_mod = prime

    # iteratively compute x0,x1,...
    for i in range(power):
        # update the target element
        h_mod = pow(h_i, n // curr_mod, p)
        
        # solve the DLP instance using BSGS
        x_i = bsgs(p, g_mod, h_mod, prime)
       
        # update the solution, x by adding xi*prime^i
        x += x_i * pow(prime, i)

        # update h for next iteration
        h_i = (h_i * pow(g, -x_i * pow(prime, i), p)) % p

        # update current modulus for next iteration
        curr_mod *= prime
    # return solution
    return x % pow(prime, power)


def crt(moduli, remainders):
    # calculate product of moduli
    product = 1
    for modulus in moduli:
        product *= modulus

    # apply the CRT
    sol = 0
    for remainder, modulus in zip(remainders, moduli):
        partial_product = product // modulus
        inverse = pow(partial_product, -1, modulus)
        sol += remainder * inverse * partial_product
    return sol % product


def bsgs(p, g, h, n):
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
