# import modules
from instance_gen_helper import find_generator, write_csv
from sys import argv
from sympy import randprime, isprime


def gen_safe_sub_dlp_instances(nMin, nMax, dstPath, file_mode='w'):
    # generate safe primes
    safe_primes = gen_n_bit_safe_primes(nMin, nMax+1)

    instances = []
    for p in safe_primes:
        # sophie-germain prime / subgroup order
        q = (p - 1) // 2

        # q is prime => smoothness(q) = q
        B = q

        # find primitive root
        g_ = find_generator(p)
        # g has order q = (p-1)/2
        g = pow(g_, 2, p)

        # n is prime => phi(n)/n-1 = 1
        phi_ratio = 1

        instances.append([p, g, q, phi_ratio, B])
    # write to file
    write_csv(instances, dstPath, mode=file_mode, delim=',')


def gen_safe_prime(pMin, pMax, attempts=10**7):
    """ Generates a prime p: p = 2q + 1, pMin < p < pMax, q prime """
    qMin, qMax = (pMin - 1) / 2, (pMax -1) / 2

    for _ in range(attempts):
        # generate a random prime, q in range
        q = randprime(qMin, qMax)

        # test if p = 2q + 1 prime also
        p = (2 * q) + 1
        if isprime(p):
            return p
    
    # no safe prime found within attempts
    print('ERROR: failed to generate safe prime')
    return None


def gen_n_bit_safe_primes(nMin, nMax):
    """ Generates one n-bit safe prime for each nMin <= n < nMax """
    safePrimes = []

    for n in range(nMin, nMax):
        # compute range for n-bit numbers
        pMin, pMax = 2 ** (n - 1), (2 ** n) - 1

        # find safe prime
        p = gen_safe_prime(pMin, pMax)

        if p is not None:
            # add safe prime to list
            safePrimes.append(p)
            print(f'SUCCESS: found {n}-bit safe prime')
        else:
            # report error
            print(f'ERROR: failed to find {n}-bit safe prime')
    return safePrimes


if __name__ == "__main__":
    # get command line args
    if len(argv) != 5:
        print("Usage: python script.py <nMin> <nMax> <output_file_path> <file_mode>")
    else:
        try:
            nMin, nMax, dstPath, file_mode = argv[1:]
            gen_safe_sub_dlp_instances(int(nMin), int(nMax), dstPath, file_mode)
            print(f"Instances generated and saved to {dstPath}")
        except ValueError:
            print("Error: nMin and nMax should be integers.")
