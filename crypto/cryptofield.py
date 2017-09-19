# -*- coding: utf-8 -*-
from __future__ import print_function
import mathTools.field as field
import mathTools.ellipticCurve as ellipticCurve
import mathTools.pairing as pairing
import mathTools.otosEC as oEC
import gmpy2 as gmpy
from Crypto.Random.random import randint
from Crypto.Random.random import getrandbits
from mathTools.otosEC import OptimAtePairing as e_hat


class CryptoField:

    def __init__(self):
        # Setting BN curve parameters
        c = gmpy.mpz(1)  # p is 256-bit long
        d = gmpy.mpz(1)
        b = c ** 4 + d ** 6  # b = c**4+d**6
        u = gmpy.mpz(-(2 ** 62 + 2 ** 55 + 1))  # p is 256-bit long

        p = CryptoField.pr(u)
        self.n = CryptoField.nr(u)

        assert gmpy.is_prime(p)
        assert gmpy.is_prime(self.n)

        # t = 6 * u ** 2 + 1

        # Fp
        Fp = field.Field(p)
        fp0 = Fp.zero()
        fp1 = Fp.one()

        # E[Fp]
        C = ellipticCurve.Curve(fp0, b * fp1, Fp)  # Y**2 = X**3+b
        PInf = ellipticCurve.ECPoint(infty=True)
        EFp = ellipticCurve.ECGroup(Fp, C, PInf)
        self.P = EFp.elem((-d ** 2) * fp1, (c ** 2) * fp1)  # P  is a generator of EFp of order n (n*P = Pinf)
        assert self.n * self.P == PInf

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
        h = 2 * p - self.n
        self.Q = u0 * h  # Q is a generator of G2 of order n
        assert self.n * self.Q == PInf2

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
        Qpr = oEC.psi(EFp12, self.Q)  # Qpr lives in E[Fp12b]
        self.Pair = pairing.Pairing(EFp, EFp12, C, self.P, self.Q, self.n, Qpr, oEC.frobenius, gamma)
        self.gt = e_hat(self.P, self.Q, self.Pair)
        self.Gt = Fp12
        self.e_hat = e_hat
        self.H = EFp2
        self.G = EFp
        self.p = p
        self.Gamma = oEC.prec_gamma(Fp12, u, c, d)
        self.EFp12=EFp12 
        #r = randint(0, int(n - 1))
        #gtr = gt ** r

        #rgp = (randint(0, int(n - 1)) * P, randint(0, int(n - 1)) * P)
        #rhp = (randint(0, int(n - 1)) * Q, randint(0, int(n - 1)) * Q)

    @staticmethod
    def pr(u):
        return 36 * u ** 4 + 36 * u ** 3 + 24 * u ** 2 + 6 * u + 1

    @staticmethod
    def nr(u):
        return 36 * u ** 4 + 36 * u ** 3 + 18 * u ** 2 + 6 * u + 1

    def e(self, g, h):
        """Evaluates the bilinear operator on pairs of elements of G and H.
        """
        r00 = e_hat(oEC.toEFp(self.G, g[0]), oEC.toEFp2(self.H, h[0]), self.Pair)
        r01 = e_hat(oEC.toEFp(self.G, g[0]), oEC.toEFp2(self.H, h[1]), self.Pair)
        r10 = e_hat(oEC.toEFp(self.G, g[1]), oEC.toEFp2(self.H, h[0]), self.Pair)
        r11 = e_hat(oEC.toEFp(self.G, g[1]), oEC.toEFp2(self.H, h[1]), self.Pair)
        return (r00, r01, r10, r11)
