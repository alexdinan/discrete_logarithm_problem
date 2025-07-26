from math import e, sqrt, log, prod, gcd
from sympy.ntheory import factorint
from random import randint


def dlp(p, g, h, n, prime_fact=None):
    if n == p - 1:
        g_ = g
    else:
        g_ = find_generator(p, prime_fact)
    
    # init. factor base
    B = int(e ** (0.6 * sqrt(log(p) * log(log(p))))) + 10
    fb = primes_up_to(B)
    fb_size = len(fb)
    
    # collect fb_size + 1 relations
    rc = RelationCollector(p, g_, B, fb)
    rels = rc.collect_rels((fb_size + 1))
    
    while True:
        # row reduce mod n (in-place)
        full_rank = row_reduce(rels, n)
        if full_rank:
            break

        # collect more relations
        rels += rc.collect_rels(int(fb_size * 0.2))

    # extract discrete logs from matrix
    fb_dls = [int(rels[i][-1]) for i in range(fb_size)]

    # compute x: g_^x ≡ g mod p
    log_gprime_g = compute_dlog(p, g_, g, fb, fb_dls, n)
    # compute x: g_^x ≡ h mod p
    log_gprime_h = compute_dlog(p, g_, h, fb, fb_dls, n)
    # combine for x: g^x ≡ h mod p
    return (log_gprime_h * pow(log_gprime_g, -1, n)) % n


def find_generator(p, prime_fact=None):
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


def primes_up_to(n):
    # not sieve[i] <=> i is composite
    sieve = [True] * n
    sieve[0] = sieve[1] = False
    p = 2

    # perform the sieve
    while (p * p) <= n:
        if sieve[p]:
            # all multiples of p are composite
            for i in range(p * p, n, p):
                sieve[i] = False
        p += 1

    # extract primes
    return [i for i, is_prime in enumerate(sieve) if is_prime]


class RelationCollector:
    def __init__(self, p, g, B, fb):
        self.p = p
        self.g = g
        self.B = B
        self.fb = fb
        self.k, self.gk = 1, g
        self.batch_size = 4 * B
        self.z = prod(fb)
    

    def collect_rels(self, num_rels):
        # initialise relation matrix
        rels = []

        while len(rels) < num_rels:
            # create batch
            batch = []
            for _ in range(self.batch_size):
                batch.append((self.gk, self.k))
                # increment
                self.k += 1
                self.gk = (self.gk * self.g) % self.p
            
            # collect rels with bernstein batch smoothness
            rels += (self.batch_smooth(batch))
        return rels
    

    def batch_smooth(self, batch):
        rels = []
        # compute z % xi using product tree
        rems = self.rems_from_tree(self.build_prod_tree([n for n, _ in batch]))

        # x = g^k mod p, r = z mod x
        for (x, k), r in zip(batch, rems):
            e2, y = 2, r
            while e2 < x and y != 0:
                # repeated squaring
                e2 **= 2
                y = (y * y) % x
            
            # x is smooth over B
            if y == 0:
                # factorise and add to relations
                rels.append(fact_over_list(x, self.fb) + [k])
        return rels
    

    def rems_from_tree(self, T):
        # computes n % v for each leaf node of T
        rems = [self.z]
        for level in T:
            rems = [rems[i // 2] % node for i, node in enumerate(level)]
        return rems
    

    def build_prod_tree(self, arr):
        # build product tree given leaf nodes
        tree = [arr]
        while len(arr) > 1:
            # build next tree level
            arr = [prod(arr[i * 2 : (i + 1) * 2]) for i in range((len(arr) + 1) // 2)]
            # prepend to list
            tree = [arr] + tree
        return tree


def row_reduce(M, n):
    # in-place row reduction of 2D matrix mod n
    rows, cols = len(M), len(M[0])
    
    curr_col = -1
    for curr_row in range(cols - 1):
        curr_col += 1

        # search for pivot element
        pivot = False
        for i in range(curr_row, rows):
            if gcd(M[i][curr_col], n) == 1:
                pivot = True
                break
        
        # no pivot found in col
        if not pivot:
            return False
                
        # swap current row and pivot row
        M[curr_row], M[i] = M[i], M[curr_row]

        # normalise curr row
        mult_inv = pow(M[curr_row][curr_col], -1, n)
        if mult_inv != 1:
            for j in range(curr_col, cols):
                M[curr_row][j] = M[curr_row][j] * mult_inv % n
        
        # eliminate non-zero entries in other rows
        for k in range(rows):
            if M[k][curr_col] != 0 and k != curr_row:
                # reduce with additive inverse
                add_inv = n - M[k][curr_col]
                M[k][curr_col] = 0
                for l in range(curr_col+1, cols):
                    M[k][l] = (M[k][l] + (M[curr_row][l] * add_inv)) % n
    # return if matrix is full rank
    return True


def compute_dlog(p, g, h, fb, fb_dls, order):
    # compute x: g^x ≡ h mod p
    while True:
        # pick a random s, compute h.g^-s
        s = randint(0, order - 1)
        g_s_h = (h * pow(g, s, p)) % p
        #print(g_s_h, h, g, s, p)

        # factorise h.g^s over the factor base
        fact = fact_over_list(g_s_h, fb)

        if fact is not None:
            # solve for x = log(h)

            x = sum(pow_f * log_f for pow_f, log_f in zip(fact, fb_dls)) - s
            return x % order
        

def fact_over_list(x, l):
    fact = [0] * len(l)
    for i, elem in enumerate(l):
        # take out all factors of each element
        while x % elem == 0:
            x //= elem
            fact[i] += 1
        
        # factorisation complete
        if x == 1:
            return fact
    # x cannot be factorised over l
    return None
