"""
# Copyright 2017 Chris Culnane
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
"""
# -*- coding: utf-8 -*-
from __future__ import print_function
import mathTools.otosEC as oEC
from dltable import DLTable
from Crypto.Random.random import randint
from Crypto.Random.random import getrandbits
from mathTools.otosEC import OptimAtePairing as e_hat
from group import Group

class TargetGroup(Group):
    """
    Target group class that provides operations in the target Group. 

    Extends Group
    """
    def __init__(self, cryptofield, dltable):
        super(TargetGroup, self).__init__(cryptofield, dltable)

    def encrypt(self, public_key, message):
        """
        Encrypt the message with the public_key
        """
        u1 = randint(1, int(self.field.p))
        u2 = randint(1, int(self.field.p))
        u3 = randint(1, int(self.field.p))

        Gt = public_key['Gt']

        gt = public_key['e']

        gt_oec = oEC.toTupleFp12(gt)
        gt_m = oEC.squareAndMultiplyFp12(Gt, gt_oec, message,
                                         oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)

        gt_ssabc=oEC.squareAndMultiplyFp12(Gt, public_key['gt_ssprime'], u1 + u2 + u3,
                                           oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)

        c3 = oEC.tmulFp12(Gt, gt_m, gt_ssabc, self.field.Gamma)
        c1 = oEC.squareAndMultiplyFp12(Gt, public_key['gt_sprime'], u3,
                                       oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)
        c2 = oEC.squareAndMultiplyFp12(Gt, public_key['gt_s'], u2,
                                       oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)
        c0 = oEC.squareAndMultiplyFp12(Gt, gt_oec, u1,
                                       oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)

        C = {'C0': c0, 'C1': c1, 'C2': c2, 'C3': c3}
        return C

    def add(self, public_key, cipher_one, cipher_two):
        """
        Add two target cipher texts together
        """
        Gt = public_key['Gt']

        G1 = public_key['G1']
        H1 = public_key['H1']
        
        c0 = oEC.tmulFp12(Gt, cipher_one['C0'], cipher_two['C0'], self.field.Gamma)
        c1 = oEC.tmulFp12(Gt, cipher_one['C1'], cipher_two['C1'], self.field.Gamma)
        c2 = oEC.tmulFp12(Gt, cipher_one['C2'], cipher_two['C2'], self.field.Gamma)
        c3 = oEC.tmulFp12(Gt, cipher_one['C3'], cipher_two['C3'], self.field.Gamma)
        C_doubleprime = {'C0': c0, 'C1': c1, 'C2': c2, 'C3': c3}
        return C_doubleprime

    def negate(self, cipher):
        """
        Not implemented, would negate a target group cipher
        """
        raise NotImplementedError('stills needs to be implemented')

    def decrypt(self, secret_key, public_key, cipher):
        """
        Decrypts a cipher text in the target Group
        """
        c0_oec = oEC.squareAndMultiplyFp12(self.field.Gt, cipher['C0'], secret_key['s']*secret_key['sprime'],
                                           oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)

        c1_oec = oEC.squareAndMultiplyFp12(self.field.Gt, cipher['C1'], secret_key['s'],
                                           oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)

        c2_oec = oEC.squareAndMultiplyFp12(self.field.Gt, cipher['C2'], secret_key['sprime'],
                                           oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)

        if cipher['C1'] == 1:
            #TODO check this condition with new structure
            #c_oec = oEC.toTupleFp12(C[0])
            c_oec = cipher[0]
        else:
            c_oec = oEC.tmulFp12(self.field.Gt, c0_oec, c1_oec, self.field.Gamma)
            c_oec = oEC.tmulFp12(self.field.Gt, c_oec, c2_oec, self.field.Gamma)
            c_oec = oEC.tmulFp12(self.field.Gt, c_oec, cipher['C3'], self.field.Gamma)
        return self.dltable.extract(oEC.toTupleFp12(public_key['e']), c_oec)
