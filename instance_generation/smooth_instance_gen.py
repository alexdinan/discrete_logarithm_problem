# import modules
from sympy import randprime, isprime
from math import log, e
from instance_gen_helper import write_csv, instance_from_p
from sys import argv


def gen_smooth_num(nMin, nMax, bMin, bMax, attempts=100):
    for _ in range(attempts):
        # initialise
        success = True
        n, fact = 1, {}

        # get first factor
        B = randprime(bMin, bMax)
        n *= B
        fact[B] = 1

        while not nMin < n < nMax:
            # upper limit to ensure n remains in range
            lim = min(bMax, nMax / n)

            # pick next factor
            if lim <= 2:
                success = False
                break
            elif lim < 3:
                p = 2
            else:
                p = randprime(2, lim)
            
            # add factor
            n *= p
            if p in fact:
                fact[p] += 1
            else:
                fact[p] = 1
        
        if success:
            return n, fact
    # failure
    return None, None


def gen_smooth_prime(pMin, pMax, bMin, bMax):
    # max_attempts derived from prime number theorem
    max_attempts = int(50 * log(pMax, e))

    for _ in range(max_attempts):
        p, fact = gen_smooth_num(pMin, pMax, bMin, bMax)

        if p is not None and isprime(p + 1):
            return p + 1, fact
    # failure
    return None, None


def gen_n_bit_smooth_primes(pMin, pMax, bMinBits, bMaxBits):
    smooth_primes = {}

    for i in range(bMinBits, bMaxBits+1):
        bMin, bMax = 2 ** i, (2 ** (i + 1)) - 1

        p, fact = gen_smooth_prime(pMin, pMax, bMin, bMax)

        if p is not None:
            smooth_primes[p] = fact
            print(f'SUCCESS: {i}-bit smooth prime found in range')
        else:
            print(f'ERROR: {i}-bit smooth prime not found')
    return smooth_primes


def gen_smooth_instances(pMin, pMax, bMinBits, bMaxBits, dst, file_mode='w'):
    # get smooth primes
    smooth_primes = gen_n_bit_smooth_primes(pMin, pMax, bMinBits, bMaxBits)

    instances = []
    for p, fact in smooth_primes.items():
        # get instance info.
        g, B, phi = instance_from_p(p, fact)
        n = p - 1
        instances.append([p, g, n, phi, B])
    write_csv(instances, dst, mode=file_mode, delim=',')


if __name__ == "__main__":
    # get command line args
    if len(argv) != 7:
        print("Usage: python script.py <pMin> <pMax> <bMinBits> <bMaxBits> <output_file_path> <file_mode>")
    else:
        try:
            pMin, pMax, bMin, bMax, dstPath, file_mode = argv[1:]
            gen_smooth_instances(int(pMin), int(pMax), int(bMin), int(bMax), dstPath, file_mode)
            print(f"Instances generated and saved to {dstPath}")
        except ValueError:
            print("Error: pMin, pMax, bMinBits, bMaxBits should be integers.")
