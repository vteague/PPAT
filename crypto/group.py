"""
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
import time
import abc
import numpy as np
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

class Group(object):
    """
    Group superclass containing utility methods and abstract methods approrpiate
    for the source and target group subclasses.
    """
    __metaclass__ = abc.ABCMeta
    MAX_BLINDING_FACTOR = 2**31
    # maximum size in bits of secret
    MAX_BITS = 100
    def __init__(self, cryptofield, dltable):
        self.field = cryptofield
        self.dltable = dltable
        self.gt = None

    def save_public_key(self, public_key, pk_file_path):
        """
        Save the public key to a JSON file

        :param public_key: public_key dictionary
        :param pk_file_path: string path to save public key to
        """
        pk_json = {}
        pk_json['g_0'] = public_key['g'][0].digits(16)
        pk_json['g_1'] = public_key['g'][1].digits(16)
        pk_json['g_2'] = public_key['g'][2]
        pk_json['h_0'] = public_key['h'][0].digits(16)
        pk_json['h_1'] = public_key['h'][1].digits(16)
        pk_json['h_2'] = public_key['h'][2].digits(16)
        pk_json['h_3'] = public_key['h'][3].digits(16)
        pk_json['h_4'] = public_key['h'][4]
        with open(pk_file_path, 'w') as outfile:
            json.dump(pk_json, outfile, indent=True, sort_keys=True)

    def load_public_key(self, pk_file_path):
        """
        load the public key from a JSON file

        :param pk_file_path: string path to load key from

        :return: JSON object contain public key
        """
        with open(pk_file_path) as data_file:
            return json.load(data_file)

    def save_secret_key(self, secret_key, sk_file_path):
        """
        Save the secret key to a JSON file

        :param secret_key: secret_key dictionary
        :param sk_file_path: string path to save secret key to
        """
        sk_json = {}
        sk_json['s'] = secret_key['s']
        sk_json['sprime'] = secret_key['sprime']
        with open(sk_file_path, 'w') as outfile:
            json.dump(sk_json, outfile, indent=True, sort_keys=True)

    def load_secret_key(self, sk_file_path):
        """
        load the secret key from a JSON file

        :param sk_file_path: string path to load key from

        :return: JSON object contain secret key
        """
        with open(sk_file_path) as data_file:
            return json.load(data_file)

    def key_gen(self, key=None):
        """
        Generate key dictionaries. If key is specified read
        from the key dictionary, otherwise generate a new key pair

        :param key: optional key dictionary containing a key pair to use

        :return: (public_key, secret_key) dictionary
        """
        if key != None:
            pub_key = key['pk']
            sec_key = key['sk']
            s = sec_key['s']
            sprime = sec_key['sprime']
            P1 = self.EFpTupleToPoint((mpz(pub_key['g_0'], base=16),
                                       mpz(pub_key['g_1'], base=16),
                                       pub_key['g_2']))
            Q1 = self.EFp2TupleToPoint((mpz(pub_key['h_0'], base=16),
                                        mpz(pub_key['h_1'], base=16),
                                        mpz(pub_key['h_2'], base=16),
                                        mpz(pub_key['h_3'], base=16),
                                        pub_key['h_4']))
            G1 = (self.field.G.neg(P1) * s, P1)  # Description of G1 - (g^{-s},g)
            H1 = (self.field.H.neg(Q1) * sprime, Q1)  # Description of H1 - (h^{-s},h)

            g = P1
            h = Q1
        else:
            s = getrandbits(self.MAX_BITS)
            sprime = getrandbits(self.MAX_BITS)

            P1 = self.field.P * randint(0, int(self.field.n))
            G1 = (self.field.G.neg(P1) * s, P1)  # Description of G1 - (g^{-s},g)

            Q1 = self.field.Q * randint(0, int(self.field.n))
            H1 = (self.field.H.neg(Q1) * sprime, Q1)  # Description of H1 - (h^{-s},h)

            g = P1
            h = Q1

        self.gt = e_hat(g, h, self.field.Pair)
        gtelem = oEC.toTupleFp12(self.gt)
        gt_s = oEC.squareAndMultiplyFp12(self.field.Gt, gtelem, s,
                                         oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)
        gt_sprime = oEC.squareAndMultiplyFp12(self.field.Gt, gtelem, sprime,
                                              oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)
        gt_ssprime = oEC.squareAndMultiplyFp12(self.field.Gt, gtelem, (s*sprime),
                                               oEC.tmulFp12, oEC.tsqrtFp12, self.field.Gamma)
        gt_ssprime = oEC.toTupleFp12(self.field.Gt.invert(oEC.toFp12elem(self.field.Gt,
                                                                         gt_ssprime)))

        # Convert G1 into oEC tuples
        G1_0_elem = oEC.toTupleEFp(self.field.G.elem(G1[0].x, G1[0].y))
        G1_1_elem = oEC.toTupleEFp(self.field.G.elem(G1[1].x, G1[1].y))
        G1Elem = (G1_0_elem, G1_1_elem)
        H1Elem = (oEC.toTupleEFp2(H1[0]), oEC.toTupleEFp2(H1[1]))

        # Convert g to EFp Tuple
        gElem = oEC.toTupleEFp(self.field.G.elem(g.x, g.y))
        # Convert h to EFp2 Tuple
        hElem = oEC.toTupleEFp2(h)

        # Pre-compute and convert G1[0] * H1[0] and G1[1] * H1[1]
        public_key = {'G': self.field.G, 'G1': G1Elem, 'H': self.field.H, 'H1': H1Elem,
                      'Gt': self.field.Gt, 'g': gElem, 'h': hElem, 'e': self.gt,
                      'gt_s': gt_s, 'gt_sprime': gt_sprime, 'gt_ssprime': gt_ssprime}

        secret_key = {'sprime':sprime, 's':s, 'ssprime':s*sprime}

        return public_key, secret_key

    def EFpTupleToPoint(self, elem):
        """
        Converts an Fp tuple to a point in Fp
        """
        return oEC.toEFp(self.field.G, elem)

    def EFp2TupleToPoint(self, elem):
        """
        Converts an Fp2 tuple to a point in Fp2
        """
        return oEC.toEFp2(self.field.H, elem)


    @abc.abstractmethod
    def encrypt(self, public_key, message):
        raise NotImplementedError('users must define encrypt to use this base class')

    @classmethod
    def generate_rand_int(cls):
        return randint(1, cls.MAX_BLINDING_FACTOR)

    @classmethod
    def switch(cls, secret_key, public_key, blinded_cipher, blinding_factor_ec, source_group, target_group):
        """
            simulates encryption switching - needs review
        """
        # Perform a decryption in the Pairing group to recover an integer (blindingfactor + m)
        blinded_plaintext = target_group.decrypt(secret_key, public_key, blinded_cipher)

        # Encrypt blinded integer in EC group
        blinded_cipher_in_ec = source_group.encrypt(public_key, blinded_plaintext)

        # Negate the blindingfactor in EC
        negated_bf = source_group.negate(blinding_factor_ec)
        #  and add to the newly encrypted value in the EC group
        # thus removing the blinding factor
        switched_cipher = source_group.add(public_key, blinded_cipher_in_ec, negated_bf)

        return switched_cipher


    @abc.abstractmethod
    def negate(self, cipher):
        """
        Negates the specified EC point in the src group. Use for subtraction
        Args:
            C (ECPoint): to negate
        """
        raise NotImplementedError('users must define negate to use this base class')

    @abc.abstractmethod
    def add(self, public_key, cipher_one, cipher_two):
        """
        Abstract method to add to cipher texts in the same group together
        """
        raise NotImplementedError('users must define add to use this base class')

    @abc.abstractmethod
    def decrypt(self, secret_key, public_key, cipher_text):
        """
        Abstract method to decrypt a cipher text
        """
        raise NotImplementedError('users must define decrypt to use this base class')

    @classmethod
    def sim_switch(cls, secret_key, public_key, cipher, sourcegrp, targetgrp):
        """
        Simulates performing a switch with 3 blinding factors representing 3 parties
        """
        blinding_factor = cls.generate_rand_int()

        cipher_bf_source = sourcegrp.encrypt(public_key, blinding_factor)
        cipher_bf_target = targetgrp.encrypt(public_key, blinding_factor)

        blinded_cipher = targetgrp.add(public_key, cipher, cipher_bf_target)

        switched_cipher = cls.switch(secret_key, public_key,
                                     blinded_cipher,
                                     cipher_bf_source, sourcegrp, targetgrp)
        return switched_cipher
