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

# -*- coding: utf-8 -*-
from __future__ import print_function
import time
import abc
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

class Group(object):
    __metaclass__ = abc.ABCMeta
    def __init__(self, cryptofield, dltable):
        self.field = cryptofield
        self.dltable = dltable

    def save_public_key(self, pk, pk_file_path):
        pk_json = {}
        pk_json['g_0'] = pk['g'][0].digits(16)
        pk_json['g_1'] = pk['g'][1].digits(16)
        pk_json['g_2'] = pk['g'][2]
        pk_json['h_0'] = pk['h'][0].digits(16)
        pk_json['h_1'] = pk['h'][1].digits(16)
        pk_json['h_2'] = pk['h'][2].digits(16)
        pk_json['h_3'] = pk['h'][3].digits(16)
        pk_json['h_4'] = pk['h'][4]
        with open(pk_file_path, 'w') as outfile:
            json.dump(pk_json, outfile, indent=True, sort_keys=True)

    def load_public_key(self, pk_file_path):
        with open(pk_file_path) as data_file:
            return json.load(data_file)

    def save_secret_key(self, sk, sk_file_path):
        sk_json = {}
        sk_json['s'] = sk['s']
        with open(sk_file_path, 'w') as outfile:
            json.dump(sk_json, outfile, indent=True, sort_keys=True)

    def load_secret_key(self, sk_file_path):
        with open(sk_file_path) as data_file:
            return json.load(data_file)

    def key_gen(self, key=None):
        if key != None:
            pubKey = key['pk']
            secKey = key['sk']
            s = secKey['s']
            P1 = self.EFpTupleToPoint((mpz(pubKey['g_0'], base=16),
                                       mpz(pubKey['g_1'], base=16),
                                       pubKey['g_2']))
            Q1 = self.EFp2TupleToPoint((mpz(pubKey['h_0'], base=16),
                                        mpz(pubKey['h_1'], base=16),
                                        mpz(pubKey['h_2'], base=16),
                                        mpz(pubKey['h_3'], base=16),
                                        pubKey['h_4']))
            G1 = (self.field.G.neg(P1) * s, P1)  # Description of G1 - (g^{-s},g)
            H1 = (self.field.H.neg(Q1) * s, Q1)  # Description of H1 - (h^{-s},h)
            g = P1
            h = Q1
        else:
            s = getrandbits(self.MAX_BITS)

            P1 = self.P * randint(0, int(self.n))
            G1 = (self.G.neg(P1) * s, P1)  # Description of G1 - (g^{-s},g)

            Q1 = self.Q * randint(0, int(self.n))
            H1 = (self.H.neg(Q1) * s, Q1)  # Description of H1 - (h^{-s},h)

            # g = P*randint(0,int(n)) # Random element of G
            # h = Q*randint(0,int(n))# Random element of H
            g = P1
            h = Q1

        self.gt = e_hat(g, h, self.field.Pair)
        G1Pair = (oEC.toTupleFp12(e_hat(G1[0], h, self.field.Pair)),
                  oEC.toTupleFp12(e_hat(G1[1], h, self.field.Pair)), oEC.toTupleFp12(self.field.Gt.one()))
        H1Pair = (oEC.toTupleFp12(e_hat(g, H1[0], self.field.Pair)),
                  oEC.toTupleFp12(e_hat(g, H1[1], self.field.Pair)), oEC.toTupleFp12(self.field.Gt.one()))

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
        G1xH1_0 = oEC.tmulFp12(self.field.Gt, G1Pair[0], H1Pair[0], self.field.Gamma)
        G1xH1_1 = oEC.tmulFp12(self.field.Gt, G1Pair[1], H1Pair[1], self.field.Gamma)
        pk = {'G': self.field.G, 'G1': G1Elem, 'H': self.field.H, 'H1': H1Elem,
              'Gt': self.field.Gt, 'g': gElem, 'h': hElem, 'e': self.gt,
              'G1xH1':(G1xH1_0, G1xH1_1), 'G1Pair':G1Pair, 'H1Pair':H1Pair}

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

        sk = {'pi_1': pi_1, 'pi_2': pi_2, 'pi_t': pi_t, 's':s, 's2':s**2}
        return pk, sk
    def EFpTupleToPoint(self, elem):
        return oEC.toEFp(self.field.G, elem)

    def EFp2TupleToPoint(self, elem):
        return oEC.toEFp2(self.field.H, elem)


    @abc.abstractmethod
    def encrypt(self, public_key, message):
        raise NotImplementedError('users must define encrypt to use this base class')

    @classmethod
    def generate_blinding_factor(cls, public_key, source_group, target_group):
        """
        Generates a blinding factor within the fixed range

        TODO:
            Decide MAX_BLINDING_FACTOR range
            Provide proof of equality
        """
        MAX_BLINDING_FACTOR = 2**9

        # create a blinding factor
        blinding_factor = randint(1, MAX_BLINDING_FACTOR)

        # Encrypt the blinding factor in the source and target
        bf_in_ec = source_group.encrypt(public_key, blinding_factor)
        bf_in_pair = target_group.encrypt(public_key, blinding_factor)

        # TODO should provide a proof that the plaintext of bf_in_ec == bf_in_pair
        return {'bfEC':bf_in_ec, 'bfPair':bf_in_pair}
    
    @classmethod
    def switch(cls, sk, pk, blinded_cipher, blinding_factor_ec, source_group, target_group):
        """
            simulates encryption switching - needs review
        """
        start = time.time()
        # Perform a decryption in the Pairing group to recover an integer (blindingfactor + m)
        blinded_plaintext = self.Dec_tgt(sk, pk, blinded_cipher, table)
        end = time.time()
        print("Dec:", end-start)

        # Encrypt blinded integer in EC group
        start = time.time()
        blinded_cipher_in_ec = self.Enc_src(pk, blinded_plaintext)
        end = time.time()
        print("Enc:", end-start)

        start = time.time()
        # Negate the blindingfactor in EC
        negated_bf = self.negate_src(blinding_factor_ec)
        #  and add to the newly encrypted value in the EC group
        # thus removing the blinding factor
        switched_cipher = self.Add_src(pk, blinded_cipher_in_ec, negated_bf)
        end = time.time()
        print("Neg:", end-start)

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
        raise NotImplementedError('users must define add to use this base class')

    @abc.abstractmethod
    def decrypt(self, sk, pk, C):
        raise NotImplementedError('users must define decrypt to use this base class')

    