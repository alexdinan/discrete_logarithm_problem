# This code is adpated from the implementation by gilcu3
# https://github.com/gilcu3/discretelog/blob/main/discretelog/linear_sieve_index_calculus.py

# import modules
from math import isqrt, log, exp, gcd
from sympy.ntheory import factorint
from functools import lru_cache, partial
from primefac import isprime
import sys


def L(p, alpha, c):
    return int(exp(c * log(p) ** alpha * log(log(p)) ** (1 - alpha)))


def default_climit(p):
    return int(3 * L(p, 0.5, 0.45)) + 10


def default_qlimit(p):
    return int(2 * L(p, 0.5, 0.475)) + 10


@lru_cache(maxsize=None)
def smooth_primes(b):
    primes = [True] * b
    pp = [2]
    for i in range(2, b, 2):
        primes[i] = False
    for i in range(3, b, 2):
        if primes[i]:
            pp += [i]
            for j in range(i**2, b, 2 * i):
                primes[j] = False
    return pp


@lru_cache(maxsize=None)
def primitive_root(q):
    o = phi(q)
    for g in range(2, q):
        if multiplicative_order(g, q) == o:
            return g


@lru_cache(maxsize=None)
def phi(q):
    ans = 1
    for p, e in factorint(q).items():
        ans *= p ** (e - 1) * (p - 1)
    return ans


@lru_cache(maxsize=None)
def multiplicative_order(g, q):
    assert gcd(g, q) == 1
    d = phi(q)
    qs = factorint(d).keys()
    while True:
        found = True
        assert pow(g, d, q) == 1
        for p in qs:
            if d % p == 0 and pow(g, d // p, q) == 1:
                d //= p
                found = False
                break
        if found:
            break
    assert pow(g, d, q) == 1, f"order failed: g={g} d={d} q={q}"
    return d


class CongruenceFinder:
    def __init__(self, H, p, opq, fbq, fbqlogs, ifbq, qlimit, np):
        self.H, self.p, self.np = H, p, np
        self.fbq, self.fbqlogs, self.qlimit = fbq, fbqlogs, qlimit
        self.opq = opq
        self.c2 = 0
        self.fb = []
        self.infb = []
        self.clog = []
        self.ifbq = ifbq
        self.eps = 1e-9

    def sieve_values(self, n):
        for c2 in range(self.c2, self.c2 + n):
            for c1, fn in linear_sieve(
                self.p,
                self.H + c2,
                self.H,
                c2 + 1,
                self.clog,
                self.eps,
                self.qlimit,
                self.fbq,
                self.fbqlogs,
            ):
                yield c2, c1, fn

    def sieve_values_parallel(self, n):
        if True:
            yield from self.sieve_values(n)
            return


    def get(self, n):
        rels, relsex = [], []
        self.infb += [None] * n
        self.clog += [
            log(c) if c > 0 else 0 for c in range(3 * self.c2, 3 * (self.c2 + n))
        ]

        for c2, c1, fn in self.sieve_values_parallel(n):
            if self.infb[c1] is None:
                self.infb[c1] = len(self.fb)
                self.fb += [c1]
            if self.infb[c2] is None:
                self.infb[c2] = len(self.fb)
                self.fb += [c2]
            cr = [0] * self.np
            if c1 == c2:
                crex = [self.infb[c1]]
                for q, e in fn.items():
                    cr[self.ifbq[q]] = e
            else:
                for q, e in fn.items():
                    cr[self.ifbq[q]] = 2 * e
                crex = [self.infb[c1], self.infb[c2]]
            rels += [cr]
            relsex += [crex]
        self.c2 += n
        return rels, relsex
    

def linear_sieve(
    p, A, B, climit=None, clog=None, eps=None, qlimit=None, fbq=None, fbqlogs=None
):
    if climit is None:
        climit = default_climit(p)
    M = A * B // p
    base = A * (B - 1) - M * p
    logbase = base // A + 1
    if clog is None:
        clog = [log(logbase + c) if c > 0 else 0 for c in range(climit)]
    if eps is None:
        eps = 1e-9
    if qlimit is None:
        qlimit = default_qlimit(p)
    if fbq is None:
        fbq = smooth_primes(qlimit)
    if fbqlogs is None:
        fbqlogs = [log(q) for q in fbq]

    sieve = [0] * climit
    num = M * p - A * B
    topp = (M + 1) * p
    topc2 = min(climit, (topp + A - 1) // A - B)
    for q, logq in zip(fbq, fbqlogs):
        if A % q == 0:
            continue
        for yi, qpow in inverse_seq(A, q):
            # assert Hc1 * y % qpow == 1
            c2 = num * yi % qpow
            if A * (B + c2) >= topp or c2 >= climit:
                break
            while c2 < topc2:
                sieve[c2] += logq
                c2 += qpow
    logA = log(A)
    cur = base
    for c2 in range(topc2):
        cur += A
        loglowerbound = logA + clog[c2]
        if loglowerbound < sieve[c2] and abs(sieve[c2] - log(cur)) < eps:
            iss, fn = is_Bsmooth(qlimit, cur)
            if iss:
                yield c2, fn


def inverse_seq(x, q):
    qpow = 1
    y, yq = 0, pow(x % q, -1, q)
    while True:
        s = (x % (qpow * q) * y - 1) // (qpow)
        k = (-s * yq) % q
        y += qpow * k
        qpow *= q
        yield y, qpow


def structured_gaussian_elimination(rels, relsex, k):
    cols = [set() for _ in range(k)]
    kp = len(rels[0])
    n = len(rels)
    colsw = [set() for _ in range(n)]
    singlerow = []
    doublerow = []
    for i, crex in enumerate(relsex):
        if len(crex) == 1:
            singlerow += [i]
        if len(crex) == 2:
            doublerow += [i]
        for c in crex:
            cols[c].add(i)
    for c in range(k):
        w = len(cols[c])
        colsw[w].add(c)

    assert len(colsw[0]) == 0

    def reduce_row(c, orig, new):
        relsex[new].remove(c)
        if len(relsex[new]) == 1:
            singlerow.append(new)
        for i in range(kp):
            rels[new][i] -= rels[orig][i]

    marked = [False] * n
    rn = n
    changed2 = True
    extrarels = 1.1
    while changed2:
        changed = True
        while changed:
            changed = False
            # delete single rows
            tmp = singlerow[:]
            if len(singlerow) > 0:
                singlerow.clear()
                changed = True
            for r in tmp:
                if len(relsex[r]) == 1 and not marked[r]:
                    c = relsex[r][0]
                    if len(cols[c]) > 1:
                        for i in cols[c]:
                            if i != r and not marked[i]:
                                reduce_row(c, r, i)
                        colsw[len(cols[c])].remove(c)
                        cols[c] = set([r])
                        colsw[1].add(c)
            # delete single columns
            torm = []
            for cr in colsw[1]:
                assert len(cols[cr]) == 1
                for r in cols[cr]:
                    if not marked[r]:
                        torm += [r]
                        marked[r] = True
                        rn -= 1
            if len(colsw[1]) > 0:
                changed = True
            for r in torm:
                for c in relsex[r]:
                    cols[c].remove(r)
                    colsw[len(cols[c]) + 1].remove(c)
                    colsw[len(cols[c])].add(c)
        changed2 = False
        for c in range(k):
            assert c in colsw[len(cols[c])], (c, cols[c], colsw[len(cols[c])])
        while rn > extrarels + k - len(colsw[0]) + kp and len(doublerow) > 0:
            r = doublerow.pop()
            if len(relsex[r]) == 2 and not marked[r]:
                changed2 = True
                marked[r] = True
                rn -= 1
                for c in relsex[r]:
                    cols[c].remove(r)
                    colsw[len(cols[c]) + 1].remove(c)
                    colsw[len(cols[c])].add(c)
    # this case handles the case that too many relations remain
    if rn > extrarels * (k - len(colsw[0]) + kp):
        assert len(doublerow) == 0
        assert len(colsw[0]) == k
        rw = []
        for i in range(n):
            if not marked[i]:
                w = len([v for v in rels[i] if v != 0])
                rw += [(w, i)]
        rw.sort()
        for _w, i in reversed(rw):
            rn -= 1
            marked[i] = True
            if rn <= kp * extrarels:
                break
    return lambda i: not marked[i]


def mod_matrix(m, p):
    for i in range(len(m)):
        for j in range(len(m[i])):
            m[i][j] %= p
            if m[i][j] > p // 2:
                m[i][j] -= p


def rels_matrix(rels, relsex, k, rf):
    m = 0
    kmap = [None] * k
    nrels = []
    for i, crex in enumerate(relsex):
        if rf(i):
            for c in crex:
                if kmap[c] is None:
                    kmap[c] = m
                    m += 1
    
    for i, (cr, crex) in enumerate(zip(rels, relsex)):
        if rf(i):
            cur = [0] * m + cr
            if len(crex) == 2:
                u, v = crex
                cur[kmap[u]] = -2
                cur[kmap[v]] = -2
            elif len(crex) == 1:
                v = crex[0]
                cur[kmap[v]] = -2
            nrels += [cur]
    return nrels, kmap


def solve_lanczos(M, q):
    b = [(-v[-1]) % q for v in M]
    M, d = to_sp([v[:-1] for v in M])
    x = block_lanczos(M, d, b, q)
    return x


def to_sp(mat):
    n, m = len(mat), len(mat[0])
    nmat = [[] for _ in range(n)]
    for i in range(n):
        for j, v in enumerate(mat[i]):
            if v != 0:
                nmat[i].append((j, v))
    return nmat, m


def sptr(mat, m):
    nmat = [[] for _ in range(m)]
    n = len(mat)
    for i in range(n):
        for j, v in mat[i]:
            nmat[j].append((i, v))
    return nmat


def spmax(mat):
    mx = 0
    for i in range(len(mat)):
        c1, c2 = 0, 0
        for _, v in mat[i]:
            if v > 0:
                c1 += v
            else:
                c2 -= v
        mx = max(mx, max(c1, c2))
    return mx


def spmul2(mx, mat, vec, q):
    n = len(mat)
    nvec = [0] * n
    xvec = vec[:]
    bt = 63 - mx.bit_length() - 1
    mask = (1 << bt) - 1
    for b in range(0, q.bit_length() + 1, bt):
        cvec = [v & mask for v in xvec]
        for i in range(n):
            a = 0
            for j, v in mat[i]:
                a += cvec[j] * v
                # assert abs(a) <= 1 << 63
            nvec[i] += a << b
        for i in range(len(xvec)):
            xvec[i] >>= bt
    return vmod(nvec, q)


def spmul(mat, vec, q):
    n = len(mat)
    nvec = [0] * n
    for i in range(n):
        a = 0
        for j, v in mat[i]:
            a += vec[j] * v
        nvec[i] = a
        # for some reason this is slower
        # nvec[i] = sum([v * vec[j] for j, v in mat[i]])
    return vmod(nvec, q)


def dot(vec1, vec2, q):
    return sum([v1 * v2 for v1, v2 in zip(vec1, vec2)]) % q


def smul(c, vec, q):
    return [c * v % q for v in vec]


def vadd(vec1, vec2):
    return [v1 + v2 for v1, v2 in zip(vec1, vec2)]


def vsub(vec1, vec2):
    return [v1 - v2 for v1, v2 in zip(vec1, vec2)]


def vmod(vec, q):
    return [c % q for c in vec]


def matmul_choice(q, m):
    # This function chooses the best matrix multiplication routine
    # In my tests spmul2 with integer unpacked works better using pypy3 for
    # bigger numbers, which probably should not be the case
    is_pypy = "PyPy" in sys.version
    if is_pypy:
        mx = spmax(m)
        if q.bit_length() + mx.bit_length() > 70 and mx.bit_length() <= 15:
            return partial(spmul2, mx)
    return spmul


def block_lanczos(mat, d, b, q):
    m, mt = mat, sptr(mat, d)

    # This is a hack to achieve better matrix multiplication performance
    rmul = matmul_choice(q, m)
    rmult = matmul_choice(q, mt)

    v0 = rmult(mt, b, q)

    def mulA(x):
        return rmult(mt, rmul(m, x, q), q)

    x = [0] * d
    denp = None
    v2, v1 = None, v0
    for i in range(1, d + 1):
        Avi = mulA(v1)
        num = dot(Avi, Avi, q)
        den1 = dot(v1, Avi, q)
        if den1 % q == 0:
            # it must finish exactly on the d-th iteration
            return None
        den1i = pow(den1, -1, q)
        ci = smul(num * den1i % q, v1, q)
        numx = dot(v1, v0, q)
        x = vmod(vadd(x, smul(numx * den1i % q, v1, q)), q)
        vi = vsub(Avi, ci)
        if i > 1:
            num = den1
            ci1 = smul(num * denp % q, v2, q)
            vi = vsub(vi, ci1)
        denp = den1i
        v2, v1 = v1, vmod(vi, q)
    if not all(c % q == 0 for c in v1):
        return None
    # assert b == vmod(spmul(mat, x), q)
    return vmod(x, q)


def compute_dlog(dlogs, pf, p):
    return sum([dlogs[q] * e for q, e in pf.items()]) % p


def individual_logs(dlogs, Hlogs, y, g, p, op, qlimit, climit):
    yy = y
    inf = min(op, 2**31)
    bound = isqrt(p) + 1
    for w in range(1, inf):
        yy = yy * g % p
        ab = ratlift(yy, bound, p)
        if ab is not None:
            a, b = ab
            sa = 1 if a > 0 else -1
            a = abs(a)
            iss, af = is_Bsmooth(qlimit, a)
            if iss:
                iss, bf = is_Bsmooth(qlimit, b)
                if iss:
                    ye = (-w + ((p - 1) // 2 if sa == 1 else 0)) % op
                    alog = compute_dlog(dlogs, af, op)
                    blog = compute_dlog(dlogs, bf, op)
                    ye = (ye + alog - blog) % op
                    ye = ye * pow((p - 1) // op, -1, op) % op
                    assert pow(g, ye * (p - 1) // op, p) == y
                    return ye
                

def ratlift(u, bnd, m):
    a1, a2 = m, u
    v1, v2 = 0, 1
    m2 = bnd
    while True:
        if v2 >= m2:
            return None
        if a2 < m2:
            return sign(v2) * a2, abs(v2)
        q = a1 // a2
        a1 -= q * a2
        v1 -= q * v2
        a1, a2, v1, v2 = a2, a1, v2, v1


def sign(a):
    if a > 0:
        return 1
    elif a < 0:
        return -1
    else:
        return 0


def is_Bsmooth(b, n):
    ps = {}
    for p in smooth_primes(b):
        if n > 1 and p**2 > n:
            if n < b:
                ps[n] = 1
                return True, ps
            else:
                return False, None
        e = 0
        while n % p == 0:
            e += 1
            n //= p
        if e > 0:
            ps[p] = e
    return n == 1, ps


def dlp(p, g, h, n, qlimit=None, climit=None):
    gy, y, op = g, h, n
    #assert isprime(p)
    #assert isprime(op)

    opq = op
    while (p - 1) % (opq * n) == 0:
        opq *= n

    # (H + c1) where 0 <= c1 <= climit
    if climit is None:
        climit = default_climit(p)

    # init. factor base
    if qlimit is None:
        qlimit = default_qlimit(p)
    fbq = smooth_primes(qlimit)
    np = len(fbq)
    # mapping between primes and indices in fbq
    ifbq = [None] * qlimit
    for i, q in enumerate(fbq):
        ifbq[q] = i
    # get (real) logs of fbq elements
    fbqlogs = [log(q) for q in fbq]

    # find g such that |<g>| = p-1
    g = primitive_root(p)

    # (H + c1) where 0 <= c1 <= climit
    H = isqrt(p) + 1
    assert H > g
    assert H + climit < p

    # relation involving g = prod(pi^ei) mod p
    rels = []
    relsex = []
    assert g <= qlimit
    if g not in fbq:
        gf = factorint(g)   # adapted
        ig = np
        np += 1
        cr = [0] * np
        for q, e in gf.items():
            assert q < qlimit
            cr[ifbq[q]] = e
        cr[-1] = -1
        rels += [cr]
        relsex += [[]]
    else:
        ig = ifbq[g]

    solved = False
    CF = CongruenceFinder(H, p, opq, fbq, fbqlogs, ifbq, qlimit, np)
    nclimit = None
    first = True
    rounds = 0
    while not solved:
        rounds += 1
        if rounds > 10:
            # solution not converging: increase qlimit -> expand factor base 
            return dlp(p, gy, y, op, qlimit + 50, climit)
        # expand climit
        nclimit = climit if first else max(50, climit // 10)
        
        nrels, nrelsex = CF.get(nclimit)
        rels += nrels
        relsex += nrelsex
        climit += nclimit if not first else 0
        first = False
        if len(nrels) == 0:
            continue
        m = np
        for i in range(climit):
            if CF.infb[i]:
                m += 1
        n = len(rels)
        if n < m + 1:
            continue

        rf = structured_gaussian_elimination(rels, relsex, len(CF.fb))
        mod_matrix(rels, opq)
        Mrels, kmap = rels_matrix(rels, relsex, len(CF.fb), rf)
        if len(Mrels) == 0:
            continue
        n, m = len(Mrels), len(Mrels[0]) - np
        if len(Mrels) < len(Mrels[0]) + 1:
            continue
        ikmap = [None] * (len(Mrels[0]) - np)
        for i, v in enumerate(kmap):
            if v is not None:
                ikmap[v] = i
        for i in range(n):
            Mrels[i][m + ig], Mrels[i][-1] = Mrels[i][-1], Mrels[i][m + ig]

        x = solve_lanczos(Mrels, opq)
        if x is None:
            continue
        
        iinfb = [None] * len(CF.fb)
        for i, v in enumerate(CF.infb):
            if v is not None:
                iinfb[v] = i
        
        gr = pow(g, (p - 1) // opq, p)
        exps = {}
        exps[g] = 1
        Hexps = [None] * climit
        solved = True
        for i in range(m + np - 1):
            if x[i] is not None:
                if i < m:
                    v = iinfb[ikmap[i]]
                    Hexps[v] = x[i]
                elif i - m != ig:
                    exps[fbq[i - m]] = x[i]
                else:
                    exps[fbq[-1]] = x[i]
            elif i >= m:
                solved = False
        
        if solved:
            for i, (rel, relex) in enumerate(zip(rels, relsex)):
                if not rf(i):
                    unknown = [c for c in relex if Hexps[iinfb[c]] is None]
                    if len(unknown) == 1:
                        other = 0
                        for c in relex:
                            if Hexps[iinfb[c]] is not None:
                                other = Hexps[iinfb[c]]
                        v = iinfb[unknown[0]]
                        pf = {fbq[j]: rel[j] for j in range(np) if rel[j] != 0}
                        Hexps[v] = (compute_dlog(exps, pf, opq) * (opq + 1) // 2 - other) % opq
    
    for q, e in exps.items():
        assert pow(gr, e, p) == pow(q, (p - 1) // opq, p)
    Hlogs = 0
    for v, e in enumerate(Hexps):
        if e is not None:
            Hlogs += 1
            assert pow(gr, e, p) == pow(H + v, (p - 1) // opq, p)
    
    # INDIVIDUAL LOG COMPUTATION
    xlogs = []
    for x in [gy, y]:
        xe = individual_logs(exps, Hexps, x, g, p, opq, qlimit, climit)

        if opq != op:
            assert xe % (opq // op) == 0
            xe //= opq // op
            assert pow(g, (p - 1) // op * xe, p) == x

        xlogs += [xe]
    ye = xlogs[1] * pow(xlogs[0], -1, op) % op
    assert pow(gy, ye, p) == y
    return ye
