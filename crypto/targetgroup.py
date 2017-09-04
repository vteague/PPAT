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
import mathTools.otosEC as oEC
from dltable import DLTable
from Crypto.Random.random import randint
from Crypto.Random.random import getrandbits
from mathTools.otosEC import OptimAtePairing as e_hat
from group import Group
class TargetGroup(Group):

    def __init__(self, cryptofield, dltable):
        super(TargetGroup, self).__init__(cryptofield, dltable)

    def encrypt(self, public_key, message):

        a = randint(1, int(self.field.p))
        b = randint(1, int(self.field.p))

        G1 = public_key['G1']
        H1 = public_key['H1']
        Gt = public_key['Gt']

        gt = public_key['e']

        g1 = (oEC.mulECP(self.field.G, G1[0], a), oEC.mulECP(self.field.G, G1[1], a))
        h1 = (oEC.mulECP(self.field.H, H1[0], b, True), oEC.mulECP(self.field.H, H1[1], b, True))

        gt_oec = oEC.toTupleFp12(gt)
        gt_m = oEC.squareAndMultiplyFp12(Gt, gt_oec, message,
                                         oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)

        gt_mXG1h1 = oEC.tmulFp12(Gt,
                                 gt_m,
                                 oEC.toTupleFp12(e_hat(oEC.toEFp(self.field.G, G1[1]),
                                                       oEC.toEFp2(self.field.H, h1[0]),
                                                       self.field.Pair)),
                                 self.field.Gamma)

        C0 = oEC.tmulFp12(Gt, gt_mXG1h1,
                          oEC.toTupleFp12(e_hat(oEC.toEFp(self.field.G, g1[0]),
                                                oEC.toEFp2(self.field.H, H1[1]),
                                                self.field.Pair)),
                          self.field.Gamma)

        C1 = oEC.tmulFp12(Gt,
                          oEC.toTupleFp12(e_hat(oEC.toEFp(self.field.G, G1[1]),
                                                oEC.toEFp2(self.field.H, h1[1]),
                                                self.field.Pair)),
                          oEC.toTupleFp12(e_hat(oEC.toEFp(self.field.G, g1[1]),
                                                oEC.toEFp2(self.field.H, H1[1]),
                                                self.field.Pair)),
                          self.field.Gamma)
        #C1 = e_hat(G1[1],h1[1],self.Pair)*e_hat(g1[1],H1[1],self.Pair)
        C2 = oEC.toTupleFp12(Gt.one())
        C = {'C0': C0, 'C1': C1, 'C2': C2}
        return C

    def add(self, public_key, cipher_one, cipher_two):

        Gt = public_key['Gt']

        G1 = public_key['G1']
        H1 = public_key['H1']

        #g = pk['g']
        #h = pk['h']
        #a = randint(1, int(self.p))
        #b = randint(1, int(self.p))

        # e(g,h1):
        #g1 = (e_hat(g, (H1[0] * a), self.Pair),
        #      e_hat(g, (H1[1] * a), self.Pair), Gt.one())
        # e(g1,h):
        #h1 = (e_hat((G1[0] * b), h, self.Pair),
         #     e_hat((G1[1] * b), h, self.Pair), Gt.one())
        C0 = oEC.tmulFp12(Gt, cipher_one['C0'], cipher_two['C0'], self.field.Gamma)
        C0 = oEC.tmulFp12(Gt, C0, public_key['H1Pair'][0], self.field.Gamma)
        C0 = oEC.tmulFp12(Gt, C0, public_key['G1Pair'][0], self.field.Gamma)

        C1 = oEC.tmulFp12(Gt, cipher_one['C1'], cipher_two['C1'], self.field.Gamma)
        C1 = oEC.tmulFp12(Gt, C1, public_key['H1Pair'][1], self.field.Gamma)
        C1 = oEC.tmulFp12(Gt, C1, public_key['G1Pair'][1], self.field.Gamma)

        C2 = oEC.tmulFp12(Gt, cipher_one['C2'], cipher_two['C2'], self.field.Gamma)
        C2 = oEC.tmulFp12(Gt, C2, public_key['H1Pair'][2], self.field.Gamma)
        C2 = oEC.tmulFp12(Gt, C2, public_key['G1Pair'][2], self.field.Gamma)

        C_doubleprime = {'C0': C0, 'C1': C1, 'C2': C2}
        return C_doubleprime

    def negate(self, cipher):
        return cipher

    def decrypt(self, sk, pk, C):
        s = sk['s']

        c1_oec = oEC.squareAndMultiplyFp12(self.field.Gt, C['C1'], s,
                                           oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)

        c2_oec = oEC.squareAndMultiplyFp12(self.field.Gt, C['C2'], sk['s2'],
                                           oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)

        if C['C1'] == 1:
            #c_oec = oEC.toTupleFp12(C[0])
            c_oec = C[0]
        else:
            c_oec = oEC.tmulFp12(self.field.Gt, C['C0'], c1_oec, self.field.Gamma)
            c_oec = oEC.tmulFp12(self.field.Gt, c_oec, c2_oec, self.field.Gamma)
        return self.dltable.extract_from_full(c_oec)
