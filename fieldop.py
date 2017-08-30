import mathTools.field as field
import mathTools.ellipticCurve as ellipticCurve
import mathTools.pairing as pairing
from mathTools.otosEC import OptimAtePairing as e_hat
import mathTools.otosEC as oEC
import gmpy2 as gmpy
from sys import getsizeof
from Crypto.Random.random import randint
# from cryptogroup import CryptoGroup

import unittest
import time

# Setting BN curve parameters
c = gmpy.mpz(1)  # p is 256-bit long
d = gmpy.mpz(1)
b = c ** 4 + d ** 6  # b = c**4+d**6
u = gmpy.mpz(-(2 ** 62 + 2 ** 55 + 1))  # p is 256-bit long


def pr(u):
    return 36 * u ** 4 + 36 * u ** 3 + 24 * u ** 2 + 6 * u + 1


def nr(u):
    return 36 * u ** 4 + 36 * u ** 3 + 18 * u ** 2 + 6 * u + 1


p = pr(u)
n = nr(u)

assert gmpy.is_prime(p)
assert gmpy.is_prime(n)

# t = 6 * u ** 2 + 1

# Fp
Fp = field.Field(p)
fp0 = Fp.zero()
fp1 = Fp.one()

# E[Fp]
C = ellipticCurve.Curve(fp0, b * fp1, Fp)  # Y**2 = X**3+b
PInf = ellipticCurve.ECPoint(infty=True)
EFp = ellipticCurve.ECGroup(Fp, C, PInf)
P = EFp.elem((-d ** 2) * fp1, (c ** 2) * fp1)  # P  is a generator of EFp of order n (n*P = Pinf)
assert n * P == PInf

# E[Fp2]
poly1 = field.polynom(Fp, [fp1, fp0, fp1])  # X**2+1
Fp2 = field.ExtensionField(Fp, poly1, rep='i')  # A**2 = -1
fp2_0 = Fp2.zero()
fp2_1 = Fp2.one()
fp2_ip = field.polynom(Fp, [fp1, fp0])  # 1*A+0
fp2_i = field.ExtensionFieldElem(Fp2, fp2_ip)
xi = (c ** 2) * fp2_1 + (d ** 3) * fp2_i  # c**2+(d**3)*A (4+i)
cxi = (c ** 2) * fp2_1 - (d ** 3) * fp2_i  # c**2-(d**3)*A
C2 = ellipticCurve.Curve(fp2_0, cxi, Fp2)  # Y**2 = X**3+c**2-(d**3)*A The twisted curve
PInf2 = ellipticCurve.ECPoint(infty=True)
EFp2 = ellipticCurve.ECGroup(Fp2, C2, PInf2)

u0 = EFp2.elem((-d) * fp2_i, c * fp2_1)  # EC point (-d*A,c)
h = 2 * p - n
Q = u0 * h  # Q is a generator of G2 of order n
assert n * Q == PInf2

# Fp6
poly3 = field.polynom(Fp2, [fp2_1, fp2_0, fp2_0, -xi])  # X**3-xi
Fp6 = field.ExtensionField(Fp2, poly3)
fp6_0 = Fp6.zero()
fp6_1 = Fp6.one()
fp6_xi = Fp6.elem(xi)  # xi in Fp6

# Fp12
poly6 = field.polynom(Fp6, [fp6_1, fp6_0, -fp6_xi])  # X**2-xi
Fp12 = field.ExtensionField(Fp6, poly6)
fp12_0 = Fp12.zero()
fp12_1 = Fp12.one()
C12 = ellipticCurve.Curve(fp12_0, b * fp12_1, Fp12)  # Y**2 = X**3+b
PInf12 = ellipticCurve.ECPoint(infty=True)
EFp12 = ellipticCurve.ECGroup(Fp12, C12, PInf12)

gamma = oEC.prec_gamma(Fp12, u, c, d)
Qpr = oEC.psi(EFp12, Q)  # Qpr lives in E[Fp12b]
Pair = pairing.Pairing(EFp, EFp12, C, P, Q, n, Qpr, oEC.frobenius, gamma)
gt = e_hat(P, Q, Pair)

r = randint(0, int(n - 1))
gtr = gt ** r

rgp = (randint(0, int(n - 1)) * P, randint(0, int(n - 1)) * P)
rhp = (randint(0, int(n - 1)) * Q, randint(0, int(n - 1)) * Q)


def en(g, h):
    """Evaluates the bilinear operator on pairs of elements of G and H.
    Arguments:
        g, a pair of elements of EFp
        h, a pair of elements of EFp2
    Returns a triple of elements of Fp12
    """
    r0 = e_hat(g[0], h[0], Pair)
    r1 = e_hat(g[0], h[1], Pair) * e_hat(g[1], h[0], Pair)
    r2 = e_hat(g[1], h[1], Pair)
    return r0, r1, r2


def e(g, h):
    """Evaluates the bilinear operator on pairs of elements of G and H.
    Uses Karatsuba's trick
    Arguments:
        g, a pair of elements of EFp
        h, a pair of elements of EFp2
    Returns a triple of elements of Fp12
    """
    r0 = e_hat(g[0], h[0], Pair)
    r2 = e_hat(g[1], h[1], Pair)
    r1 = e_hat(g[0] + g[1], h[0] + h[1], Pair) * Fp12.invert(r0 * r2)

    return r0, r1, r2

def make_EFptable(group, base, max_dl=2 ** 32, max_search=2 ** 12):
    """This function makes a multiplication table to aid in computing discrete
    logarithms to speed up decryption of multiple messages encrypted with the
    same public/private key
    :param group is an EFp group
    :param base is a tuple
    :param max_dl is the highest value that the DL can take
    :param max_search is the maximum number of steps that we agree to make during DL extraction
    """
    # Store the giant steps. Keys are truncated x coordinate of points, values are exponents

    giant_steps = {}
    table_size = max_dl / max_search + 1
    # Size of the giant steps (on the curve)
    giant_step = oEC.mulECP(group, base, max_search)
    # Counter of current value of the exponent
    exponent = max_search
    # Position of that counter on the curve
    running_step = giant_step
    # j performs as many steps as needed.
    # The '+1' handles the case when max_dl is not a square
    for j in xrange(table_size):
        lsb_running_step = int(gmpy.t_mod_2exp(running_step[0], 128))
        giant_steps[lsb_running_step] = exponent
        running_step = oEC.addEFp(group, running_step, giant_step)
        exponent += max_search
        # print("Giant steps: ", giant_steps)
    return giant_steps


def log_group_EFp(Group, a, b, table):
    """Extracts the discrete log of b in base a in Group, which must be an EFp.
    Assumes table contains precomputed values for base a

    :param Group: EFp Group in which a and b lie
    :param a: base element as a tuple
    :param b: tuple from which DL must be extracted
    :param table: precomputed table for a
    :return: x: b == a**x
    """
    if table is None:
        table = make_EFptable(Group, a)

    i = 0
    while not gmpy.t_mod_2exp(b[0], 128) in table:
        b = oEC.addEFp(Group, a, b)
        i += 1
    return table[gmpy.t_mod_2exp(b[0], 128)] - i


class TestFieldOp(unittest.TestCase):
    def setUp(self):
        self.startTime = time.time()

    def test_MulEFp(self):
        _ = r * P
        t = time.time() - self.startTime
        print "%s: %.4f" % ("EFp point multiplication -- generic", t)

    def test_MulEFp_opt(self):
        _ = oEC.mulECP(EFp, oEC.toTupleEFp(P, Jcoord=True), r, sq=False, Jcoord=True)
        t = time.time() - self.startTime
        print "%s: %.4f" % ("EFp point multiplication -- optimized", t)

    def test_MulEFp2(self):
        _ = r * Q
        t = time.time() - self.startTime
        print "%s: %.4f" % ("EFp2 point multiplication -- generic", t)

    def test_MulEFp2_opt(self):
        _ = oEC.mulECP(EFp2, oEC.toTupleEFp2(Q), r, sq=True)
        t = time.time() - self.startTime
        print "%s: %.4f" % ("EFp2 point multiplication -- optimized", t)

    def test_ExpFp12(self):
        _ = gt ** r
        t = time.time() - self.startTime
        print "%s: %.4f" % ("Fp12 point exponentiation -- generic", t)

    def test_ExpFp12_opt(self):
        mul = oEC.tmulFp12
        sqrt = oEC.tsqrtFp12
        sqmu = oEC.squareAndMultiply
        tgt = oEC.toTupleFp12(gt)
        tgtr = sqmu(Fp12, tgt, r, mul, sqrt, gamma)

        t = time.time() - self.startTime
        print "%s: %.4f" % ("Fp12 point exponentiation -- optimized", t)
        self.assertEqual(oEC.toFp12elem(Fp12, tgtr), gt ** r,
                         'Optimised Fp12 exp gives inconsistent results')

    def test_InvFp12(self):
        igtr = Fp12.invert(gtr)

        t = time.time() - self.startTime
        print "%s: %.4f" % ("Fp12 point inversion -- generic", t)
        self.assertEqual(igtr * gtr, Fp12.one(),
                         'Inversion in Fp12 does not work')

    def test_pairing(self):
        _ = e_hat(P, Q, Pair)
        t = time.time() - self.startTime
        print "%s: %.4f" % ("Pairing", t)

    def test_bilinearity(self):
        gta1 = e_hat(r * P, Q, Pair)
        gta2 = e_hat(P, r * Q, Pair)
        gta3 = gt ** r
        self.assertEqual(gta1, gta3, 'Not bilinear wrt G1')
        self.assertEqual(gta2, gta3, 'Not bilinear wrt G2')

    def test_ppairing(self):
        gt1 = en(rgp, rhp)
        t = time.time() - self.startTime
        print "%s: %.4f" % ("pairing on pairs -- generic", t)
        t1 = time.time()
        gt2 = e(rgp, rhp)
        t = time.time() - t1
        print "%s: %.4f" % ("pairing on pairs -- Karatsuba", t)
        self.assertEqual(gt1, gt2, 'pairing on pairs seems inconsistent')

    def test_EFp_DLog(self):
        Ptuple = oEC.toTupleEFp(P)
        table = make_EFptable(EFp, Ptuple, max_dl=2 ** 32, max_search=2 ** 12)
        t = time.time() - self.startTime
        print "%s: %.4f" % ("ECTable computation", t)
        s = getsizeof(table)
        print "%s: %d %s" % ("Size of ECTable", s, "bytes")
        x = 9
        y = oEC.mulECP(EFp, Ptuple, x, sq=False)
        #print y
        t2 = time.time()
        xc = log_group_EFp(EFp, Ptuple, y, table)
        t = time.time() - t2
        print "%s: %.4f" % ("DL Extraction", t)
        self.assertEqual(x, xc, 'DL extraction fails')


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFieldOp)
    unittest.TextTestRunner(verbosity=1).run(suite)
