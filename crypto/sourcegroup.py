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
from __future__ import print_function
import time
import numpy as np
import itertools
import hashlib
import json
import mathTools.field as field
import mathTools.ellipticCurve as ellipticCurve
import mathTools.pairing as pairing
import ppat.ppats
import mathTools.otosEC as oEC
from dltable import DLTable
import gmpy2 as gmpy
from gmpy2 import mpz
from Crypto.Random.random import randint
from Crypto.Random.random import getrandbits
from mathTools.otosEC import OptimAtePairing as e_hat
from math import sqrt
from group import Group

class SourceGroup(Group):
    # maximum size in bits of secret
    MAX_BITS = 100

    def __init__(self, cryptofield, dltable):
        super(SourceGroup, self).__init__(cryptofield, dltable)

    def encrypt(self, public_key, message):
        a = randint(1, int(self.field.p))
        b = randint(1, int(self.field.p))

        G1 = public_key['G1']
        H1 = public_key['H1']

        g1 = (oEC.mulECP(self.field.G, G1[0], a), oEC.mulECP(self.field.G, G1[1], a))
        h1 = (oEC.mulECP(self.field.H, H1[0], b, True), oEC.mulECP(self.field.H, H1[1], b, True))

        C_0 = (oEC.addEFp(self.field.G, oEC.mulECP(self.field.G, public_key['g'], message), g1[0]), g1[1])
        C_1 = (oEC.addEFp2(self.field.H, oEC.mulECP(self.field.H, public_key['h'], message, True), h1[0]), h1[1])

        C = {'C0': C_0, 'C1': C_1}
        return C

    def multiply(self, public_key, cipher_one, cipher_two):
        c3 = e_hat(oEC.toEFp(self.field.G, cipher_one['C0'][0]), oEC.toEFp2(self.field.H, cipher_two['C1'][0]), self.field.Pair)
        c2 = e_hat(oEC.toEFp(self.field.G, cipher_one['C0'][0]), oEC.toEFp2(self.field.H, cipher_two['C1'][1]), self.field.Pair)
        c1 = e_hat(oEC.toEFp(self.field.G, cipher_one['C0'][1]), oEC.toEFp2(self.field.H, cipher_two['C1'][0]), self.field.Pair)
        c0 = e_hat(oEC.toEFp(self.field.G, cipher_one['C0'][1]), oEC.toEFp2(self.field.H, cipher_two['C1'][1]), self.field.Pair)
    
        #c_one = oEC.tmulFp12(self.field.Gt, c_one, G1xH1[1], self.field.Gamma)

        cipher_prime = {'C0': oEC.toTupleFp12(c0), 'C1': oEC.toTupleFp12(c1), 'C2':oEC.toTupleFp12(c2), 'C3':oEC.toTupleFp12(c3)}

        return cipher_prime  # C = e(C0,C1) * e(g,h1) * e(g1,h)

    def negate(self, cipher):
        """
        Negates the specified EC point in the src group. Use for subtraction
        Args:
            C (ECPoint): to negate
        """
        # straightforward EC negation
        # We first convert the tuple to an ECPoint in the appropriate group
        # We then convert back to a tuple
        C0_0 = oEC.toTupleEFp(oEC.toEFp(self.field.G, cipher['C0'][0]).__neg__())
        C0_1 = oEC.toTupleEFp(oEC.toEFp(self.field.G, cipher['C0'][1]).__neg__())
        C1_0 = oEC.toTupleEFp2(oEC.toEFp2(self.field.H, cipher['C1'][0]).__neg__())
        C1_1 = oEC.toTupleEFp2(oEC.toEFp2(self.field.H, cipher['C1'][1]).__neg__())

        return {'C0': (C0_0, C0_1), 'C1': (C1_0, C1_1)}


    def add(self, public_key, cipher_one, ciher_two):

        G1 = public_key['G1']
        H1 = public_key['H1']

        #a = randint(1, int(self.n))
        #b = randint(1, int(self.n))
        #g1 = (G1[0]*a, G1[1]*a)
        #h1 = (H1[0]*b, H1[1]*b)
        c_zero = (oEC.addEFp(self.field.G, oEC.addEFp(self.field.G, ciher_two['C0'][0], G1[0]),
                             cipher_one['C0'][0]),
                  oEC.addEFp(self.field.G, oEC.addEFp(self.field.G, ciher_two['C0'][1], G1[1]),
                             cipher_one['C0'][1]))
        c_one = (oEC.addEFp2(self.field.G, oEC.addEFp2(self.field.G, ciher_two['C1'][0], H1[0]),
                             cipher_one['C1'][0]),
                 oEC.addEFp2(self.field.G, oEC.addEFp2(self.field.G, ciher_two['C1'][1], H1[1]),
                             cipher_one['C1'][1]))

        c_doubleprime = {'C0':c_zero, 'C1':c_one} # C'' = ( C0*C'0*g1 , C1*C'1*h1 )
        return c_doubleprime
    
    def decrypt(self, sk, pk, C):
        g = self.EFpTupleToPoint(pk['g'])
        c0_0 = self.EFpTupleToPoint(C['C0'][0])
        c0_1 = self.EFpTupleToPoint(C['C0'][1])
        return self.dltable.extract(pk['g'], oEC.toTupleEFp(sk['pi_1']((c0_0, c0_1))))
