# import modules
from sympy.ntheory import factorint


def find_generator(p, prime_fact=None):
    """ Finds a generator, g of Zp* """
    # get prime factorisation of p - 1
    if prime_fact is None:
        prime_fact = factorint(p - 1)
    
    # try candidate values for g
    for g in range(2, p):
        is_gen = True

        for prime in prime_fact.keys():
            if pow(g, (p-1)//prime, p) == 1:
                # g is not a generator - exit
                is_gen = False
                break
        if is_gen:
            return g
    

def phi(n, prime_fact=None):
    """ Computes euler's totient function for n """
    # get prime factorisation of n
    if prime_fact is None:
        prime_fact = factorint(n)

    # compute phi(n)
    phi_n = 1
    for prime, power in prime_fact.items():
        phi_n *= (prime ** (power - 1)) * (prime - 1)
    
    # return 0 < phi(n)/n-1 <= 1
    return phi_n / (n - 1)


def smoothness(n, prime_fact=None):
    """ Finds the smoothness, B of an integer, n """
    # get prime factorisation of n
    if prime_fact is None:
        prime_fact = factorint(n)
    # return largest prime factor
    return max(prime_fact.keys())


def instance_from_p(p, prime_fact=None):
    """ generates a dlp instance from a prime, p """
    # get prime factorisation
    if prime_fact is None:
        prime_fact = factorint(p - 1)

    # find generator
    g = find_generator(p, prime_fact)

    # find smoothness
    B = smoothness(p - 1, prime_fact)

    # find euler's totient ratio
    phi_n = phi(p - 1, prime_fact)
    return g, B, phi_n


def write_csv(data, dst_path, mode='w', delim=','):
    """ write data to a csv file """
    # convert data to correct string format
    output = ""
    for row in data:
        output += delim.join(map(str, row)) + '\n'

    # write to file
    with open(dst_path, mode) as f:
        f.write(output)
    print('SUCCESS: data written to file')
