# Copyright 2017 Ilya Marchenko
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# -*- coding: utf-8 -*-
import time
import numpy as np
import itertools
import mathTools.field as field
import mathTools.ellipticCurve as ellipticCurve
import mathTools.pairing as pairing
import ppat.ppats
import mathTools.otosEC as oEC
import gmpy2 as gmpy
from Crypto.Random.random import randint
from Crypto.Random.random import getrandbits
from mathTools.otosEC import OptimAtePairing as e_hat


class CryptoGroup:
    # maximum size in bits of secret
    MAX_BITS = 100

    def __init__(self):
        # c = gmpy.mpz(2) # p is 160-bit long
        c = gmpy.mpz(1)  # p is 256-bit long
        d = gmpy.mpz(1)
        b = c**4 + d**6  # b = c**4+d**6
        # u = gmpy.mpz(-(2**38 + 2**28 + 1 )) # p is 160-bit long
        u = gmpy.mpz(-(2**62 + 2**55 + 1))  # p is 256-bit long
        #p = 36*u**4 + 36*u**3 + 24*u**2 + 6*u + 1
        p = self.pr(u)
        self.n = self.nr(u)
        #p = 1300829
        #n = 1299721

        #n = 36*u**4 + 36*u**3 + 18*u**2 + 6*u + 1
        # n is 160-bit long with low HW
        t = 6 * u**2 + 1

        ##### Fp #####
        Fp = field.Field(p)
        fp0 = Fp.zero()
        fp1 = Fp.one()

        print (Fp, " ...done")
        ##### E[Fp] #####
        C = ellipticCurve.Curve(fp0, b * fp1, Fp)  # Y**2 = X**3+b
        PInf = ellipticCurve.ECPoint(infty=True)
        EFp = ellipticCurve.ECGroup(Fp, C, PInf)
        # P  is a generetor of EFp of order n (n*P = Pinf)
        self.P = EFp.elem((-d**2) * fp1, (c**2) * fp1)

        ##### Fp2b #####
        poly1 = field.polynom(Fp, [fp1, fp0, fp1])  # X**2+1
        print (poly1)

        Fp2 = field.ExtensionField(Fp, poly1, rep='i')  # A**2 = -1
        print (Fp2, " ...done")
        fp2_0 = Fp2.zero()
        fp2_1 = Fp2.one()
        fp2_ip = field.polynom(Fp, [fp1, fp0])  # 1*A+0
        fp2_i = field.ExtensionFieldElem(Fp2, fp2_ip)
        xi = (c**2) * fp2_1 + (d**3) * fp2_i  # c**2+(d**3)*A (4+i)
        cxi = (c**2) * fp2_1 - (d**3) * fp2_i  # c**2-(d**3)*A
        # ixi = 8*fp2bi-8*fp2b1 # 8*A-8
        #xi = ixi.invert()
        # C2b = EllipticCurve.Curve(fp2b0, 3*ixi,Fp2b) # Y**2 = X**3+3*(8*A-8)
        # Y**2 = X**3+c**2-(d**3)*A The twisted curve
        C2 = ellipticCurve.Curve(fp2_0, cxi, Fp2)
        PInf2 = ellipticCurve.ECPoint(infty=True)
        EFp2 = ellipticCurve.ECGroup(Fp2, C2, PInf2)

        u0 = EFp2.elem((-d) * fp2_i, c * fp2_1)  # EC point (-d*A,c)
        h = 2 * p - self.n
        self.Q = u0 * h  # Q is a generator of G2 of order n

        r = randint(1, int(self.n))
        s = randint(1, int(self.n))
        rP = r * self.P
        sQ = s * self.Q

        ##### Fp6 #####
        poly3 = field.polynom(Fp2, [fp2_1, fp2_0, fp2_0, -xi])  # X**3-xi
        Fp6 = field.ExtensionField(Fp2, poly3)
        fp6_0 = Fp6.zero()
        fp6_1 = Fp6.one()
        fp6_xi = Fp6.elem(xi)  # xi in Fp6

        ##### Fp12 #####
        poly6 = field.polynom(Fp6, [fp6_1, fp6_0, -fp6_xi])  # X**2-xi
        Fp12 = field.ExtensionField(Fp6, poly6)
        print(Fp12, " ...done")
        fp12_0 = Fp12.zero()
        fp12_1 = Fp12.one()
        C12 = ellipticCurve.Curve(fp12_0, b * fp12_1, Fp12)  # Y**2 = X**3+b
        PInf12 = ellipticCurve.ECPoint(infty=True)
        EFp12 = ellipticCurve.ECGroup(Fp12, C12, PInf12)

        Qpr = oEC.psi(EFp12, self.Q)  # Qpr lives in E[Fp12b]
        self.Pair = pairing.Pairing(EFp, EFp12, C, self.P, self.Q, self.n, Qpr,
                                    oEC.frobenius, oEC.prec_gamma(Fp12, u, c, d))
        self.Gt = Fp12
        self.e_hat = e_hat
        self.H = EFp2
        self.G = EFp
        self.p = p
        self.Gamma = oEC.prec_gamma(Fp12, u, c, d)
    def pr(self, u):
        return 36 * u**4 + 36 * u**3 + 24 * u**2 + 6 * u + 1

    def nr(self, u):
        return 36 * u**4 + 36 * u**3 + 18 * u**2 + 6 * u + 1

    # Range for iterating over large numbers (larger than int)
    def range(self, stop): return iter(itertools.count().next, stop)

    ############### Crypto Functions #############
    

    def e(self, g, h):
        """Evaluates the bilinear operator on pairs of elements of G and H.
        Uses Karatsuba's trick.
        Arguments:
            pp, the public parameters describing the groups
            g, a pair of elements of G
            h, a pair of elements of H
        Returns a triple of elements of Gt
        """
        r0 = e_hat(g[0], h[0], self.Pair)
        r2 = e_hat(g[1], h[1], self.Pair)
        r1 = e_hat(g[0] + g[1], h[0] + h[1], self.Pair) * \
            self.Gt.invert(r0 * r2)
        return (r0, r1, r2)

    def KeyGen(self):

        s = getrandbits(self.MAX_BITS)

        P1 = self.P * randint(0, int(self.n))
        G1 = (self.G.neg(P1) * s, P1)  # Description of G1 - (g^{-s},g)

        Q1 = self.Q * randint(0, int(self.n))
        H1 = (self.H.neg(Q1) * s, Q1)  # Description of H1 - (h^{-s},h)

        # g = P*randint(0,int(n)) # Random element of G
        # h = Q*randint(0,int(n))# Random element of H
        g = P1
        h = Q1

        gt = e_hat(g, h, self.Pair)

        pk = {'G': self.G, 'G1': G1, 'H': self.H, 'H1': H1,
              'Gt': self.Gt, 'g': g, 'h': h, 'e': gt}

        def pi_1(g):
            if g[1] == 1:
                output = g[0]
            else:
                output = g[0] + (g[1] * s)  # pi_1((g1,g2)) = g1 * g2^s
            return output

        def pi_2(h):
            if h[1] == 1:
                output = h[0]
            else:
                output = h[0] + (h[1] * s)  # pi_2((h1,h2)) = h1 * h2^s
            return output

        def pi_t(gt):
            if gt['C1'] == 1:
                output = gt[0]
            else:
                # pi_t((gt1,gt2,gt3)) = gt1 * gt2^s * gt3^{s^2}
                output = gt['C0'] * (gt['C1']**s) * (gt['C2']**(s**2))
            return output

        sk = {'pi_1': pi_1, 'pi_2': pi_2, 'pi_t': pi_t}
        return pk, sk

    def Enc_src(self, pk, M):

        a = randint(1, int(self.p))
        b = randint(1, int(self.p))

        G1 = pk['G1']
        H1 = pk['H1']

        g1 = (G1[0] * a, G1[1] * a)
        h1 = (H1[0] * b, H1[1] * b)

        C_0 = (pk['g'] * M + g1[0], g1[1])  # C_0 = (g^M * g^{-as}, g^a)
        C_1 = (pk['h'] * M + h1[0], h1[1])  # C_1 = (h^M * h^{-bs}, h^b)

        C = {'C0': C_0, 'C1': C_1}
        return C
    def Enc_src_neg(self, pk, M):

        a = randint(1, int(self.p))
        b = randint(1, int(self.p))

        G1 = pk['G1']
        H1 = pk['H1']

        g1 = (G1[0] * a, G1[1] * a)
        h1 = (H1[0] * b, H1[1] * b)

        C_0 = ((pk['g'] * M + g1[0]).__neg__(), g1[1].__neg__())  # C_0 = (g^M * g^{-as}, g^a)
        C_1 = ((pk['h'] * M + h1[0]).__neg__(), h1[1].__neg__())  # C_1 = (h^M * h^{-bs}, h^b)

        C = {'C0': C_0, 'C1': C_1}
        return C

    def Enc_tgt(self, pk, M):

        a = randint(1, int(self.p))
        b = randint(1, int(self.p))

        G1 = pk['G1']
        H1 = pk['H1']
        Gt = pk['Gt']

        gt = pk['e']

        g1 = (G1[0] * a, G1[1] * a)
        h1 = (H1[0] * b, H1[1] * b)

        C0 = (gt**M)*e_hat(G1[1],h1[0],self.Pair)*e_hat(g1[0],H1[1],self.Pair)
        C1 = e_hat(G1[1],h1[1],self.Pair)*e_hat(g1[1],H1[1],self.Pair)
        C2 = Gt.one()

        C = {'C0': C0, 'C1': C1, 'C2': C2}
        return C

    def Multiply_src(self, pk, C0, C1):

        Gt = pk['Gt']

        G1 = pk['G1']
        H1 = pk['H1']

        g = pk['g']
        h = pk['h']
        a = randint(1, int(self.p))
        b = randint(1, int(self.p))

        # e(C0,C1):
        #eC = (e_hat(C0[0],C1[0],Pair),e_hat(C0[0],C1[1],Pair)*e_hat(C0[1],C1[0],Pair),e_hat(C0[1],C1[1],Pair))
        eC = self.e(C0, C1)
        # e(g,h1):
        g1 = (e_hat(g, H1[0] * a, self.Pair),
              e_hat(g, H1[1] * a, self.Pair), Gt.one())
        # e(g1,h):
        h1 = (e_hat(G1[0] * b, h, self.Pair),
              e_hat(G1[1] * b, h, self.Pair), Gt.one())

        C = {'C0': (eC[0] * g1[0]) * h1[0], 'C1': eC[1] * (g1[1] * h1[1]), 'C2': eC[2]}

        return C  # C = e(C0,C1) * e(g,h1) * e(g1,h)


    def switch(self, sk, pk, blinded_cipher, blinding_factor_ec, table):
        # simulates encryption switching - needs review

        # Perform a decryption in the Pairing group to recover an integer (blindingfactor + m)
        blinded_plaintext = self.Dec_tgt(sk, pk, blinded_cipher, table)

        # Encrypt blinded integer in EC group
        blinded_cipher_in_ec = self.Enc_src(pk, blinded_plaintext)

        # Negate the blindingfactor in EC
        negated_bf = self.negate_src(blinding_factor_ec)
        #  and add to the newly encrypted value in the EC group
        # thus removing the blinding factor
        switched_cipher = self.Add_src(pk, blinded_cipher_in_ec, negated_bf)

        return switched_cipher

    def blind_pair_cipher(self, pk, C, blinding_factors):
        # accumulate blinding factors
        bf_in_ec = blinding_factors[0]['bfEC']
        bf_in_pair = blinding_factors[0]['bfPair']
        for bf_count in range(1, len(blinding_factors)):
            bf_in_ec = self.Add_src(pk, bf_in_ec, blinding_factors[bf_count]['bfEC'])
            bf_in_pair = self.Add_tgt(pk, bf_in_pair, blinding_factors[bf_count]['bfPair'])

        # Add the blinding factor the ciphertext
        return {'C':self.Add_tgt(pk, C, bf_in_pair), 'bfEC':bf_in_ec}

    def generate_blinding_factor(self, pk):
        MAX_BLINDING_FACTOR = 2**14
        
        # create a blinding factor
        blinding_factor = randint(1, MAX_BLINDING_FACTOR)

        # Encrypt the blinding factor in the source and target
        bf_in_ec = self.Enc_src(pk, blinding_factor)
        bf_in_pair = self.Enc_tgt(pk, blinding_factor)

        # TODO should provide a proof that the plaintext of bf_in_ec == bf_in_pair
        return {'bfEC':bf_in_ec, 'bfPair':bf_in_pair}

    def negate_src(self, C):
        # straightforward EC negation
        return {'C0': (C['C0'][0].__neg__(),C['C0'][1].__neg__()), 'C1': (C['C1'][0].__neg__(),C['C1'][1].__neg__())}

    def Add_src(self, pk, C, C_prime):

        G1 = pk['G1']
        H1 = pk['H1']

        a = randint(1, int(self.n))
        b = randint(1, int(self.n))
        g1 = (G1[0]*a, G1[1]*a)
        h1 = (H1[0]*b, H1[1]*b)

        C_0 = ((C_prime['C0'][0]+g1[0])+C['C0'][0],(C_prime['C0'][1]+g1[1])+C['C0'][1])
        C_1 = ((C_prime['C1'][0]+h1[0])+C['C1'][0],(C_prime['C1'][1]+h1[1])+C['C1'][1])
        
        C_doubleprime = {'C0':C_0,'C1':C_1} # C'' = ( C0*C'0*g1 , C1*C'1*h1 )
        return C_doubleprime

    def Add_tgt(self, pk, C, C_prime):

        Gt = pk['Gt']

        G1 = pk['G1']
        H1 = pk['H1']

        g = pk['g']
        h = pk['h']
        a = randint(1, int(self.p))
        b = randint(1, int(self.p))

        # e(g,h1):
        g1 = (e_hat(g, (H1[0] * a), self.Pair),
              e_hat(g, (H1[1] * a), self.Pair), Gt.one())
        # e(g1,h):
        h1 = (e_hat((G1[0] * b), h, self.Pair),
              e_hat((G1[1] * b), h, self.Pair), Gt.one())

        C0 = C['C0'] * C_prime['C0'] * g1[0] * h1[0]
        C1 = C['C1'] * C_prime['C1'] * g1[1] * h1[1]
        C2 = C['C2'] * C_prime['C2'] * g1[2] * h1[2]

        # C'' = e(C,C') * e(g,h1) * e(g1,h)
        C_doubleprime = {'C0': C0, 'C1': C1, 'C2': C2}
        return C_doubleprime

    # Computes x such that a^x = b over a group of order p using baby step-giant step
    def log_group(self, a, b, Group, table={}):

        #baby_steps = {}

        x = Group.neg(a) * (2**16)
        gamma = b

        if table == {}:
            print('No Table Found')
            table = make_ECtable(Group, a)

        for i in range(2**16):
            if gamma in table:
                return i * (2**16) + table[gamma]
            else:
                gamma = gamma + x

        return "No Match"

    # Computes x such that a^x = b over a field of order p using baby step-giant step
    def log_field(self, a, b, Field, table={}):

        #baby_steps = {}

        x = oEC.toTupleFp12(a.invert()**(2**16))
        gamma = oEC.toTupleFp12(b)

        if table == {}:
            print('No Table Found')
            table = make_Ftable(Field, a)

        for i in range(2**16):
            if gamma in table:
                return (i) * (2**16) + table[gamma]
            else:
                gamma = oEC.tmulFp12(self.Gt, gamma, x, self.Gamma)

        return "No Match"

    def Dec_src(self, sk, pk, C, table={}):

        M = self.log_group(sk['pi_1']((pk['g'], 1)), sk['pi_1']
                      (C['C0']), pk['G'], table)

        return M
        # if sk['pi_2'](C['C1']).__eq__(pk['H'].mul(sk['pi_2']((pk['h'],1)),M)):
        #    return M
        # else:
        #    print "M != M'"
        #    return 3.14159 # Dummy value returned

    def Dec_tgt(self, sk, pk, C, table={}):

        g = pk['g']
        h = pk['h']

        M = self.log_field(e_hat(g, h, self.Pair), sk['pi_t'](C), pk['Gt'], table)
        return M

    def make_ECtable(self, grp, point):
        # This function makes a multiplication table to aid in computing discrete
        # logarithims to speed up decryption of multiple messages encrypted with the
        # same public/private key

        # VT: I edited this to make it terminate at exactly the index it previously did.
        # However, I'd assume that the omission of the j=2**16 + 1 case makes no
        # difference so we should just let the range terminate at 2**16 as before.
        # But I'm reluctant to fiddle in case I mess up something important.

        baby_steps = {}
        pt = point

        for j in xrange(2**16 + 1):
            baby_steps[pt] = j + 1
            pt = pt + point

        return baby_steps

    def make_Ftable(self, field, elt):
        # This function makes a multiplication table to aid in computing discrete
        # logarithims to speed up decryption of multiple messages encrypted with the
        # same public/private key

        baby_steps = {}
        gt = oEC.toTupleFp12(field.one())

        for j in xrange(2**16):
            gt = oEC.tmulFp12(self.Gt, gt, oEC.toTupleFp12(elt), self.Gamma)
            baby_steps[gt] = j + 1

        return baby_steps