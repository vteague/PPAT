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
# import sys
# sys.path.append('C:\Research\STAR Vote\PPAT-master')

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
from mathTools.otosEC import OptimAtePairing as e_hat

range = lambda stop: iter(itertools.count().next, stop)  # Range for iterating over large numbers (larger than int)

# c = gmpy.mpz(2) # p is 160-bit long
c = gmpy.mpz(1)  # p is 256-bit long
d = gmpy.mpz(1)
b = c ** 4 + d ** 6  # b = c**4+d**6
# u = gmpy.mpz(-(2**38 + 2**28 + 1 )) # p is 160-bit long
u = gmpy.mpz(-(2 ** 62 + 2 ** 55 + 1))  # p is 256-bit long


# p = 36*u**4 + 36*u**3 + 24*u**2 + 6*u + 1
def pr(u):
    return 36 * u ** 4 + 36 * u ** 3 + 24 * u ** 2 + 6 * u + 1


def nr(u):
    return 36 * u ** 4 + 36 * u ** 3 + 18 * u ** 2 + 6 * u + 1


p = pr(u)
n = nr(u)

# p = 1300829
# n = 1299721


# n = 36*u**4 + 36*u**3 + 18*u**2 + 6*u + 1
# n is 160-bit long with low HW
t = 6 * u ** 2 + 1

##### Fp #####
Fp = field.Field(p)
fp0 = Fp.zero()
fp1 = Fp.one()

print Fp, " ...done"
##### E[Fp] #####
C = ellipticCurve.Curve(fp0, b * fp1, Fp)  # Y**2 = X**3+b
PInf = ellipticCurve.ECPoint(infty=True)
EFp = ellipticCurve.ECGroup(Fp, C, PInf)
P = EFp.elem((-d ** 2) * fp1, (c ** 2) * fp1)  # P  is a generator of EFp of order n (n*P = Pinf)

##### Fp2b #####
poly1 = field.polynom(Fp, [fp1, fp0, fp1])  # X**2+1
print poly1

Fp2 = field.ExtensionField(Fp, poly1, rep='i')  # A**2 = -1
print Fp2, " ...done"
fp2_0 = Fp2.zero()
fp2_1 = Fp2.one()
fp2_ip = field.polynom(Fp, [fp1, fp0])  # 1*A+0
fp2_i = field.ExtensionFieldElem(Fp2, fp2_ip)
xi = (c ** 2) * fp2_1 + (d ** 3) * fp2_i  # c**2+(d**3)*A (4+i)
cxi = (c ** 2) * fp2_1 - (d ** 3) * fp2_i  # c**2-(d**3)*A
# ixi = 8*fp2bi-8*fp2b1 # 8*A-8
# xi = ixi.invert()
# C2b = EllipticCurve.Curve(fp2b0, 3*ixi,Fp2b) # Y**2 = X**3+3*(8*A-8)
C2 = ellipticCurve.Curve(fp2_0, cxi, Fp2)  # Y**2 = X**3+c**2-(d**3)*A The twisted curve
PInf2 = ellipticCurve.ECPoint(infty=True)
EFp2 = ellipticCurve.ECGroup(Fp2, C2, PInf2)

u0 = EFp2.elem((-d) * fp2_i, c * fp2_1)  # EC point (-d*A,c)
h = 2 * p - n
Q = u0 * h  # Q is a generator of G2 of order n

r = randint(1, int(n))
s = randint(1, int(n))
rP = r * P
sQ = s * Q

##### Fp6 #####
poly3 = field.polynom(Fp2, [fp2_1, fp2_0, fp2_0, -xi])  # X**3-xi
Fp6 = field.ExtensionField(Fp2, poly3)
fp6_0 = Fp6.zero()
fp6_1 = Fp6.one()
fp6_xi = Fp6.elem(xi)  # xi in Fp6

##### Fp12 #####
poly6 = field.polynom(Fp6, [fp6_1, fp6_0, -fp6_xi])  # X**2-xi
Fp12 = field.ExtensionField(Fp6, poly6)
print Fp12, " ...done"
fp12_0 = Fp12.zero()
fp12_1 = Fp12.one()
C12 = ellipticCurve.Curve(fp12_0, b * fp12_1, Fp12)  # Y**2 = X**3+b
PInf12 = ellipticCurve.ECPoint(infty=True)
EFp12 = ellipticCurve.ECGroup(Fp12, C12, PInf12)

Qpr = oEC.psi(EFp12, Q)  # Qpr lives in E[Fp12b]
Pair = pairing.Pairing(EFp, EFp12, C, P, Q, n, Qpr, oEC.frobenius, oEC.prec_gamma(Fp12, u, c, d))

############### PPATS ########################
x1 = randint(1, int(n - 1));
print "x1 is", x1
g1 = x1 * Q
h1td = randint(1, int(n - 1));
print "h1 trapdoor is", h1td
h1 = h1td * P
ppatspp = ppat.ppats.PPATSPublicParameters(P, Q, Pair, 'Ate', optim=True)
print 'public parameters ppatspp created'
ppatspk = ppat.ppats.PPATSPublicKey(ppatspp, g1, h1)
print 'public key ppatspk created'
ppatssk = ppat.ppats.PPATSPrivateKey(ppatspp, ppatspk, x1)
print 'secret key ppatssk created'

############### PPATC ########################
x1c = randint(1, int(n - 1));
print "x1c is", x1c
x2c = randint(1, int(n - 1));
print "x2c is", x2c
g1c = x1c * Q
g2c = x2c * Q
h1tdc = randint(1, int(n - 1));
print "h1tdc trapdoor is", h1tdc
h1c = h1tdc * P
ppatcpp = ppat.ppatc.PPATCPublicParameters(P, Q, Pair, 'Ate', psi=None, optim=True)
print 'public parameters (complex) ppatcpp created'
ppatcpk = ppat.ppatc.PPATCPublicKey(ppatcpp, g1c, g2c, h1c)
print 'public key (complex) ppatcpk created'
ppatcsk = ppat.ppatc.PPATCPrivateKey(ppatcpp, ppatcpk, x1c, x2c)
print 'secret key (complex) ppatcsk created'

############### Crypto Functions #############
Gamma = oEC.prec_gamma(Fp12, u, c, d)


def Setup():
    pp = (p, ppatcpp.h, ppatcpp.g, Fp12, ppatcpp.e)
    return pp


def KeyGen(pp):
    G = pp[1]  # Generator of G
    H = pp[2]  # Generator of H
    s = randint(1, int(p))

    P = EFp.dbleAndAdd(G, randint(1, int(n)))
    G1 = (EFp.dbleAndAdd(EFp.neg(P), s), P)  # Description of G1 - (g^{-s},g)

    Q = EFp2.dbleAndAdd(H, randint(1, int(n)))
    H1 = (EFp2.dbleAndAdd(EFp2.neg(Q), s), Q)  # Description of H1 - (h^{-s},h)

    Gt = pp[3]

    g = EFp.dbleAndAdd(G, randint(1, int(n)))  # Random element of G
    h = EFp2.dbleAndAdd(H, randint(1, int(n)))  # Random element of H

    pk = (G, G1, H, H1, Gt, g, h)

    def pi_1(g):
        if g[1] == 1:
            output = g[0]
        else:
            output = EFp.dadd(g[0], EFp.dbleAndAdd(g[1], s))  # pi_1((g1,g2)) = g1 * g2^s
        return output

    def pi_2(h):
        if h[1] == 1:
            output = h[0]
        else:
            output = EFp2.dadd(h[0], EFp.dbleAndAdd(h[1], s))  # pi_2((h1,h2)) = h1 * h2^s
        return output

    def pi_t(gt):
        if gt[1] == 1:
            output = gt[0]
        else:
            output = Fp12.mul(gt[0], Fp12.mul(Fp12.powop(gt[1], s), Fp12.powop(Fp12.powop(gt[2], s),
                                                                               s)))  # pi_t((gt1,gt2,gt3)) = gt1 * gt2^s * gt3^{s^2}
        return output

    sk = (pi_1, pi_2, pi_t)
    return pk, sk


def Enc_src(pk, M):
    a = randint(1, int(p))
    b = randint(1, int(p))

    g = pk[1]
    h = pk[3]

    g1 = (EFp.dbleAndAdd(g[0], a), EFp.dbleAndAdd(g[1], a))
    h1 = (EFp2.dbleAndAdd(h[0], b), EFp2.dbleAndAdd(h[1], b))

    C_0 = (EFp.dadd(EFp.dbleAndAdd(pk[-2], M), g1[0]), g1[1])  # C_0 = (g^M * g^{-as}, g^a)
    C_1 = (EFp2.dadd(EFp2.dbleAndAdd(pk[-1], M), h1[0]), h1[1])  # C_1 = (h^M * h^{-bs}, h^b)
    return (C_0, C_1)


def Enc_tgt(pk, M):
    # from mathTools.otosEC import OptimAtePairing as e_hat
    # from mathTools.otosEC import OptimTatePairing as e_hat

    a = randint(1, int(p))
    b = randint(1, int(p))

    G1 = pk[1]
    H1 = pk[3]

    g = pk[-2]
    h = pk[-1]

    g1 = (EFp.dbleAndAdd(G1[0], a), EFp.dbleAndAdd(G1[1], a))
    h1 = (EFp2.dbleAndAdd(H1[0], b), EFp2.dbleAndAdd(H1[1], b))

    C0 = Fp12.mul(Fp12.powop(e_hat(g, h, Pair), M), Fp12.mul(e_hat(G1[1], h1[0], Pair), e_hat(g1[0], H1[1], Pair)))
    C1 = Fp12.mul(e_hat(G1[1], h1[1], Pair), e_hat(g1[1], H1[1], Pair))
    C2 = Fp12.one()

    C = C0, C1, C2
    return C


def Multiply_src(pk, C0, C1):
    g = pk[-2]
    h = pk[-1]

    G1 = pk[1]
    H1 = pk[3]

    # e(C0,C1):
    eC = (e_hat(C0[0], C1[0], Pair), Fp12.mul(e_hat(C0[0], C1[1], Pair), e_hat(C0[1], C1[0], Pair)),
          e_hat(C0[1], C1[1], Pair))
    # e(g,h1):
    g1 = (e_hat(g, EFp2.dbleAndAdd(H1[0], randint(1, int(p))), Pair),
          e_hat(g, EFp2.dbleAndAdd(H1[1], randint(1, int(p))), Pair), Fp12.one())
    # e(g1,h):
    h1 = (e_hat(EFp.dbleAndAdd(G1[0], randint(1, int(p))), h, Pair),
          e_hat(EFp.dbleAndAdd(G1[1], randint(1, int(p))), h, Pair), Fp12.one())

    C = (Fp12.mul(Fp12.mul(eC[0], g1[0]), h1[0]), Fp12.mul(eC[1], Fp12.mul(g1[1], h1[1])), eC[2])

    return C  # C = e(C0,C1) * e(g,h1) * e(g1,h)


def Add_src(pk, C, C_prime):
    g1 = (EFp.dbleAndAdd(pk[1][0], randint(1, int(n))), EFp.dbleAndAdd(pk[1][1], randint(1, int(n))))
    h1 = (EFp2.dbleAndAdd(pk[3][0], randint(1, int(n))), EFp2.dbleAndAdd(pk[3][1], randint(1, int(n))))

    C_0 = (EFp.dadd(EFp.dadd(C_prime[0][0], g1[0]), C[0][0]), EFp.dadd(EFp.dadd(C_prime[0][1], g1[1]), C[0][1]))
    C_1 = (EFp2.dadd(EFp2.dadd(C_prime[1][0], h1[0]), C[1][0]), EFp2.dadd(EFp2.dadd(C_prime[1][1], h1[1]), C[1][1]))

    C_doubleprime = (C_0, C_1)  # C'' = ( C0*C'0*g1 , C1*C'1*h1 )

    return C_doubleprime


def Add_tgt(pk, C, C_prime):
    g = pk[-2]
    h = pk[-1]

    G1 = pk[1]
    H1 = pk[3]

    # e(g,h1):
    g1 = (e_hat(g, EFp2.dbleAndAdd(H1[0], randint(1, int(p))), Pair),
          e_hat(g, EFp2.dbleAndAdd(H1[1], randint(1, int(p))), Pair), Fp12.one())
    # e(g1,h):
    h1 = (e_hat(EFp.dbleAndAdd(G1[0], randint(1, int(p))), h, Pair),
          e_hat(EFp.dbleAndAdd(G1[1], randint(1, int(p))), h, Pair), Fp12.one())

    C0 = Fp12.mul(Fp12.mul(Fp12.mul(C[0], C_prime[0]), g1[0]), h1[0])
    C1 = Fp12.mul(Fp12.mul(Fp12.mul(C[1], C_prime[1]), g1[1]), h1[1])
    C2 = Fp12.mul(Fp12.mul(Fp12.mul(C[2], C_prime[2]), g1[2]), h1[2])

    C_doubleprime = C0, C1, C2  # C'' = e(C,C') * e(g,h1) * e(g1,h)
    return C_doubleprime


def log_group(a, b, Group, table={}):  # Computes x such that a^x = b over a group of order p using baby step-giant step

    baby_steps = {}

    x = Group.dbleAndAdd(Group.neg(a), 2 ** 16)
    gamma = b

    if table == {}:
        table = make_ECtable(Group, b) # TODO: looks like there is something wrong here: table is unused, and using b looks odd.  
    else:
        baby_steps = table

    for i in range(2 ** 16):
        if gamma in baby_steps:
            return i * (2 ** 16) + baby_steps[gamma]
        else:
            gamma = Group.dadd(gamma, x)

    return "No Match"


def log_field(a, b, Field, table={}):  # Computes x such that a^x = b over a field of order p using baby step-giant step

    baby_steps = {}

    x = oEC.toTupleFp12(Field.powop(a.invert(), 2 ** 16))
    gamma = oEC.toTupleFp12(b)

    if table == {}:
        print 'No Table Found'
        table = make_Ftable(Field, b)
    else:
        baby_steps = table

    for i in range(2 ** 16):
        if gamma in baby_steps:
            return (i) * (2 ** 16) + baby_steps[gamma]
        else:
            gamma = oEC.tmulFp12(Fp12, gamma, x, Gamma)

    return "No Match"


def Dec_src(sk, pk, C, table={}):
    M = log_group(sk[0]((pk[-2], 1)), sk[0](C[0]), EFp, table)

    if sk[1](C[1]).__eq__(EFp2.dbleAndAdd(sk[1]((pk[-1], 1)), M)):
        return M
    else:
        print "M != M'"
        return 3.14159  # Dummy value returned


def Dec_tgt(sk, pk, C, table={}):
    g = pk[-2]
    h = pk[-1]

    M = log_field(e_hat(g, h, Pair), sk[2](C), Fp12, table)
    return M


def make_ECtable(Group, point):
    # This function makes a multiplication table to aid in computing discrete
    # logarithims to speed up decryption of multiple messages encrypted with the
    # same public/private key

    baby_steps = {}
    pt = point

    for j in range(2 ** 16):
        pt = Group.add(pt, point)
        baby_steps[pt] = j + 2

    return baby_steps


def make_Ftable(Field, elt):
    # This function makes a multiplication table to aid in computing discrete
    # logarithims to speed up decryption of multiple messages encrypted with the
    # same public/private key

    baby_steps = {}
    e = oEC.toTupleFp12(Field.one())

    for j in range(2 ** 16):
        e = oEC.tmulFp12(Fp12, e, oEC.toTupleFp12(elt), Gamma)
        baby_steps[e] = j + 1

    return baby_steps


# def make_random_ballot(c,k,pk):
#    
#    M = np.ones((k,c),dtype = tuple)    
#    
#    chosen_candidates = []
#    for i in range(k):
#        x = np.random.choice([p for p in range(c) if p not in chosen_candidates])
#        M[i][x] = 2
#        chosen_candidates.append(x)
#
#    for j in range(k):
#        for l in range(c):
#            M[j][l] = Enc_src(pk, M[j][l])
#    return M        


############ Testing Crypto Functions ##########
pp = Setup()
print 'Setup Complete'

pk, sk = KeyGen(pp)
print 'KeyGen Complete'

M0 = 123
M1 = 2 ** 30

C0_src = Enc_src(pk, M0)
C1_src = Enc_src(pk, M1)
print 'Enc_src Complete'

C0_tgt = Enc_tgt(pk, M0)
C1_tgt = Enc_tgt(pk, M1)
print 'Enc_tgt Complete'

ECtable = make_ECtable(EFp, pk[-2])
Ftable = make_Ftable(Fp12, e_hat(pk[-2], pk[-1], Pair))
print 'Multiplication Tables Created'

CM = Multiply_src(pk, C0_src[0], C0_src[1])
print 'Multiply_src Complete'

C_doubleprime_src = Add_src(pk, C0_src, C1_src)
print 'Add_src Complete'

C_doubleprime_tgt = Add_tgt(pk, C0_tgt, C1_tgt)
print 'Add_tgt Complete'

t1 = time.clock()
print Dec_src(sk, pk, C0_src, ECtable)
print Dec_src(sk, pk, C1_src, ECtable)
print 'Dec_src Complete'
t2 = time.clock()
print "Source decryption time:", (t2 - t1), "seconds."

t3 = time.clock()
print Dec_tgt(sk, pk, C0_tgt, Ftable)
print Dec_tgt(sk, pk, C1_tgt, Ftable)
print 'Dec_tgt Complete'
t4 = time.clock()
print "Target decryption time:", (t4 - t3), "seconds."
