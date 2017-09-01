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

class CryptoSwitch:

    # maximum size in bits of secret
    MAX_BITS = 100

    def __init__(self, cryptofield, dltable):
        self.field = cryptofield
        self.dltable = dltable

    def save_public_key(self, public_key, pk_file_path):
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
        with open(pk_file_path) as data_file:
            return json.load(data_file)

    def save_secret_key(self, secret_key, sk_file_path):
        sk_json = {}
        sk_json['s'] = secret_key['s']
        with open(sk_file_path, 'w') as outfile:
            json.dump(sk_json, outfile, indent=True, sort_keys=True)

    def load_secret_key(self, sk_file_path):
        with open(sk_file_path) as data_file:
            return json.load(data_file)

    def key_gen(self, key=None):
        if key != None:
            public_key = key['pk']
            secret_key = key['sk']['s']
            P1 = self.EFpTupleToPoint((mpz(public_key['g_0'], base=16),
                                       mpz(public_key['g_1'], base=16),
                                       public_key['g_2']))
            Q1 = self.EFp2TupleToPoint((mpz(public_key['h_0'], base=16),
                                        mpz(public_key['h_1'], base=16),
                                        mpz(public_key['h_2'], base=16),
                                        mpz(public_key['h_3'], base=16),
                                        public_key['h_4']))
            G1 = (self.G.neg(P1) * secret_key, P1)  # Description of G1 - (g^{-s},g)
            H1 = (self.H.neg(Q1) * secret_key, Q1)  # Description of H1 - (h^{-s},h)
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

        self.gt = e_hat(g, h, self.Pair)
        G1Pair = (oEC.toTupleFp12(e_hat(G1[0], h, self.Pair)),
                  oEC.toTupleFp12(e_hat(G1[1], h, self.Pair)), oEC.toTupleFp12(self.Gt.one()))
        H1Pair = (oEC.toTupleFp12(e_hat(g, H1[0], self.Pair)),
                  oEC.toTupleFp12(e_hat(g, H1[1], self.Pair)), oEC.toTupleFp12(self.Gt.one()))

        # Convert G1 into oEC tuples
        G1_0_elem = oEC.toTupleEFp(self.G.elem(G1[0].x, G1[0].y))
        G1_1_elem = oEC.toTupleEFp(self.G.elem(G1[1].x, G1[1].y))
        G1Elem = (G1_0_elem, G1_1_elem)
        H1Elem = (oEC.toTupleEFp2(H1[0]), oEC.toTupleEFp2(H1[1]))

        # Convert g to EFp Tuple
        gElem = oEC.toTupleEFp(self.G.elem(g.x, g.y))
        # Convert h to EFp2 Tuple
        hElem = oEC.toTupleEFp2(h)

        # Pre-compute and convert G1[0] * H1[0] and G1[1] * H1[1]
        G1xH1_0 = oEC.tmulFp12(self.Gt, G1Pair[0], H1Pair[0], self.Gamma)
        G1xH1_1 = oEC.tmulFp12(self.Gt, G1Pair[1], H1Pair[1], self.Gamma)
        pk = {'G': self.G, 'G1': G1Elem, 'H': self.H, 'H1': H1Elem,
              'Gt': self.Gt, 'g': gElem, 'h': hElem, 'e': self.gt,
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

    def Enc_src(self, pk, M):
        a = randint(1, int(self.p))
        b = randint(1, int(self.p))

        G1 = pk['G1']
        H1 = pk['H1']

        g1 = (oEC.mulECP(self.G, G1[0], a), oEC.mulECP(self.G, G1[1], a))
        h1 = (oEC.mulECP(self.H, H1[0], b, True), oEC.mulECP(self.H, H1[1], b, True))

        C_0 = (oEC.addEFp(self.G, oEC.mulECP(self.G, pk['g'], M), g1[0]), g1[1])
        C_1 = (oEC.addEFp2(self.H, oEC.mulECP(self.H, pk['h'], M, True), h1[0]), h1[1])

        C = {'C0': C_0, 'C1': C_1}
        return C

    def Enc_tgt(self, pk, M):

        a = randint(1, int(self.p))
        b = randint(1, int(self.p))

        G1 = pk['G1']
        H1 = pk['H1']
        Gt = pk['Gt']

        gt = pk['e']

        g1 = (oEC.mulECP(self.G, G1[0], a), oEC.mulECP(self.G, G1[1], a))
        h1 = (oEC.mulECP(self.H, H1[0], b, True), oEC.mulECP(self.H, H1[1], b, True))

        gt_oec = oEC.toTupleFp12(gt)
        gt_m = oEC.squareAndMultiplyFp12(self.Gt, gt_oec, M,
                                         oEC.tmulFp12, oEC.tsqrtFp12, self.Gamma)

        gt_mXG1h1 = oEC.tmulFp12(self.Gt,
                                 gt_m,
                                 oEC.toTupleFp12(e_hat(oEC.toEFp(self.G, G1[1]),
                                                       oEC.toEFp2(self.H, h1[0]), self.Pair)),
                                 self.Gamma)

        C0 = oEC.tmulFp12(self.Gt, gt_mXG1h1,
                          oEC.toTupleFp12(e_hat(oEC.toEFp(self.G, g1[0]),
                                                oEC.toEFp2(self.H, H1[1]),
                                                self.Pair)),
                          self.Gamma)

        C1 = oEC.tmulFp12(self.Gt,
                          oEC.toTupleFp12(e_hat(oEC.toEFp(self.G, G1[1]),
                                                oEC.toEFp2(self.H, h1[1]),
                                                self.Pair)),
                          oEC.toTupleFp12(e_hat(oEC.toEFp(self.G, g1[1]),
                                                oEC.toEFp2(self.H, H1[1]),
                                                self.Pair)),
                          self.Gamma)
        #C1 = e_hat(G1[1],h1[1],self.Pair)*e_hat(g1[1],H1[1],self.Pair)
        C2 = oEC.toTupleFp12(Gt.one())
        C = {'C0': C0, 'C1': C1, 'C2': C2}
        return C

    def Multiply_src(self, pk, C0, C1):

        Gt = pk['Gt']

        G1xH1 = pk['G1xH1']

        eC = self.e(C0, C1)

        ec0 = oEC.toTupleFp12(eC[0])
        ec1 = oEC.toTupleFp12(eC[1])
        ec2 = oEC.toTupleFp12(eC[2])

        c0 = oEC.tmulFp12(self.Gt, ec0, G1xH1[0], self.Gamma)
        c1 = oEC.tmulFp12(self.Gt, ec1, G1xH1[1], self.Gamma)

        C = {'C0': c0, 'C1': c1, 'C2': ec2}

        return C  # C = e(C0,C1) * e(g,h1) * e(g1,h)


    def sim_switch(self, sk, pk, cipher):
        """
        Simulates performing a switch with 3 blinding factors representing 3 parties
        """
        start = time.time()
        blinding_factors = []
        #for bfcount in range(0, 3):
        startone = time.time()
        blinding_factors.append(self.generate_blinding_factor(pk))
        endone = time.time()
        print("genblind", endone-startone)

        startone = time.time()
        blinded_cipher = self.blind_pair_cipher(pk, cipher, blinding_factors)
        endone = time.time()
        print("blind", endone-startone)

        switched_cipher = self.switch(sk, pk,
                                      blinded_cipher['C'],
                                      blinded_cipher['bfEC'], self.Ftable)
        end = time.time()
        print("switch:", end-start)
        return switched_cipher

    def switch(self, sk, pk, blinded_cipher, blinding_factor_ec, table):
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

    def blind_pair_cipher(self, pk, C, blinding_factors):
        """
        Blind a cipher in the target (Pairing based) group. This is just a matter
        of summing the encrypted blinding factors in both the EC and Pairing based groups
        then adding the pairing based sum to the cipher.

        Returns:
            Dictionary containing:
                C: Blinded Cipher (Pairing group)
                bfEC: Cipher of sum of encrypted blinding factors (EC Group)
        """
        # accumulate blinding factors
        bf_in_ec = blinding_factors[0]['bfEC']
        bf_in_pair = blinding_factors[0]['bfPair']
        for bf_count in range(1, len(blinding_factors)):
            bf_in_ec = self.Add_src(pk, bf_in_ec, blinding_factors[bf_count]['bfEC'])
            bf_in_pair = self.Add_tgt(pk, bf_in_pair, blinding_factors[bf_count]['bfPair'])

        # Add the blinding factor the ciphertext
        return {'C':self.Add_tgt(pk, C, bf_in_pair), 'bfEC':bf_in_ec}

    def EFpTupleToPoint(self, elem):
        return oEC.toEFp(self.G, elem)
    
    def EFp2TupleToPoint(self, elem):
        return oEC.toEFp2(self.H, elem)
    
    def generate_blinding_factor(self, pk):
        """
        Generates a blinding factor within the fixed range

        TODO:
            Decide MAX_BLINDING_FACTOR range
            Provide proof of equality
        """
        MAX_BLINDING_FACTOR = 2**19

        # create a blinding factor
        blinding_factor = randint(1, MAX_BLINDING_FACTOR)

        # Encrypt the blinding factor in the source and target
        bf_in_ec = self.Enc_src(pk, blinding_factor)
        bf_in_pair = self.Enc_tgt(pk, blinding_factor)

        # TODO should provide a proof that the plaintext of bf_in_ec == bf_in_pair
        return {'bfEC':bf_in_ec, 'bfPair':bf_in_pair}

    def negate_src(self, C):
        """
        Negates the specified EC point in the src group. Use for subtraction
        Args:
            C (ECPoint): to negate
        """
        # straightforward EC negation
        # We first convert the tuple to an ECPoint in the appropriate group
        # We then convert back to a tuple
        C0_0 = oEC.toTupleEFp(oEC.toEFp(self.G, C['C0'][0]).__neg__())
        C0_1 = oEC.toTupleEFp(oEC.toEFp(self.G, C['C0'][1]).__neg__())
        C1_0 = oEC.toTupleEFp2(oEC.toEFp2(self.H, C['C1'][0]).__neg__())
        C1_1 = oEC.toTupleEFp2(oEC.toEFp2(self.H, C['C1'][1]).__neg__())

        return {'C0': (C0_0, C0_1), 'C1': (C1_0, C1_1)}


    def Add_src(self, pk, C, C_prime):


        G1 = pk['G1']
        H1 = pk['H1']

        #a = randint(1, int(self.n))
        #b = randint(1, int(self.n))
        #g1 = (G1[0]*a, G1[1]*a)
        #h1 = (H1[0]*b, H1[1]*b)
        C_0 = (oEC.addEFp(self.G, oEC.addEFp(self.G, C_prime['C0'][0], G1[0]), C['C0'][0]),
               oEC.addEFp(self.G, oEC.addEFp(self.G, C_prime['C0'][1], G1[1]), C['C0'][1]))
        C_1 = (oEC.addEFp2(self.G, oEC.addEFp2(self.G, C_prime['C1'][0], H1[0]), C['C1'][0]),
               oEC.addEFp2(self.G, oEC.addEFp2(self.G, C_prime['C1'][1], H1[1]), C['C1'][1]))

        C_doubleprime = {'C0':C_0, 'C1':C_1} # C'' = ( C0*C'0*g1 , C1*C'1*h1 )
        return C_doubleprime

    def Add_tgt(self, pk, C, C_prime):

        Gt = pk['Gt']

        G1 = pk['G1']
        H1 = pk['H1']

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
        C0 = oEC.tmulFp12(self.Gt, C['C0'], C_prime['C0'], self.Gamma)
        C0 = oEC.tmulFp12(self.Gt, C0, pk['H1Pair'][0], self.Gamma)
        C0 = oEC.tmulFp12(self.Gt, C0, pk['G1Pair'][0], self.Gamma)

        C1 = oEC.tmulFp12(self.Gt, C['C1'], C_prime['C1'], self.Gamma)
        C1 = oEC.tmulFp12(self.Gt, C1, pk['H1Pair'][1], self.Gamma)
        C1 = oEC.tmulFp12(self.Gt, C1, pk['G1Pair'][1], self.Gamma)

        C2 = oEC.tmulFp12(self.Gt, C['C2'], C_prime['C2'], self.Gamma)
        C2 = oEC.tmulFp12(self.Gt, C2, pk['H1Pair'][2], self.Gamma)
        C2 = oEC.tmulFp12(self.Gt, C2, pk['G1Pair'][2], self.Gamma)

        C_doubleprime = {'C0': C0, 'C1': C1, 'C2': C2}
        return C_doubleprime

    # Computes x such that a^x = b over a group of order p using baby step-giant step
    def log_groupnew(self, a, b, Group, table={}):
        table = self.ECtable
        #baby_steps = {}

        x = Group.neg(self.EFpTupleToPoint(a)) * (2**16)
        x = oEC.toTupleEFp(Group.elem(x.x, x.y))
        gamma = b

        if table == {}:
            print('No Table Found')
            table = make_ECtable(Group, a)

        for i in range(2**16):
            if gamma in table:
                return i * (2**16) + table[gamma]
            else:
                gamma = oEC.addEFp(Group, gamma, x)
                #gamma = gamma + x

        return "No Match"

    def log_group(self, a, b, Group):
        """Extracts the discrete log of b in base a in Group.
        Assumes that self.ECtable contains precomputed values for base a

        :param a: base element as a tuple
        :param b: element from which DL must be extracted as a tuple
        :param Group: Group in which a and b lie
        :param table: precomputed table for a
        :return: x: b == a**x
        """
        table = self.ECtable
        i = 0
        while not b in table:
            b = oEC.addEFp(Group, a, b)
            i += 1
            if (i % 2**13 == 0):
                print("i =", i)
        return table[b] - i

    # Computes x such that a^x = b over a field of order p using baby step-giant step
    def log_field(self, a, b, Field, table={}):
        table = self.Ftable
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
        g = self.EFpTupleToPoint(pk['g'])
        c0_0=self.EFpTupleToPoint(C['C0'][0])
        c0_1=self.EFpTupleToPoint(C['C0'][1])
        
        M = self.log_group(sk['pi_1']((g, 1)), sk['pi_1']
                           ((c0_0,c0_1)), pk['G'], table)
        
        return M
        """
        s = sk['s']
        a = pk['g']
        partialDec = oEC.addEFp(pk['G'], C['C0'][0], oEC.mulECP(pk['G'], C['C0'][1], s))
        #M = self.log_group(pk['g'], partialDec, pk['G'], table)
        #self.EFpTupleToPoint(partialDec).__neg__()
        x = self.EFpTupleToPoint(pk['g']).__neg__()
        #x = pk['G'].neg(self.EFpTupleToPoint(pk['g']))
        partialDec = oEC.addEFp(pk['G'], partialDec, oEC.toTupleEFp(pk['G'].elem(x.x, x.y)))
        print("partiaal", partialDec)
        M = self.log_full_group(partialDec)

        #M = self.log_group(sk['pi_1']((pk['g'], 1)), sk['pi_1']
        #                   (C['C0']), pk['G'], table)

        return M
        # if sk['pi_2'](C['C1']).__eq__(pk['H'].mul(sk['pi_2']((pk['h'],1)),M)):
        #    return M
        # else:
        #    print "M != M'"
        #    return 3.14159 # Dummy value returned
"""
    def Dec_tgt(self, sk, pk, C, table={}):

        s = sk['s']

        c1_oec = oEC.squareAndMultiplyFp12(self.Gt, C['C1'], s,
                                           oEC.tmulFp12, oEC.tsqrtFp12, self.Gamma)

        c2_oec = oEC.squareAndMultiplyFp12(self.Gt, C['C2'], sk['s2'],
                                           oEC.tmulFp12, oEC.tsqrtFp12, self.Gamma)

        if C['C1'] == 1:
            #c_oec = oEC.toTupleFp12(C[0])
            c_oec = C[0]
        else:
            c_oec = oEC.tmulFp12(self.Gt, C['C0'], c1_oec, self.Gamma)
            c_oec = oEC.tmulFp12(self.Gt, c_oec, c2_oec, self.Gamma)
        M = self.log_full_field(c_oec)
        return M

    def make_ECtable(self, grp, point, max_dl=2**32, max_search=2**12):
        """This function makes a multiplication table to aid in computing discrete
        logarithms to speed up decryption of multiple messages encrypted with the
        same public/private key
        Inputs:
        - grp is an elliptic curve group on Fp
        - point is a tuple
        - max_dl is the highest value that the DL can take
        """

        # Store the giant steps. Keys are points, values are exponents
        giant_steps = {}
        table_size = max_dl / max_search + 1
        # Size of the giant steps (on the curve)
        giant_step = oEC.mulECP(grp, point, max_search)
        # Counter of current value of the exponent
        exponent = max_search
        # Position of that counter on the curve
        running_step = giant_step
        # j performs as many steps as needed.
        # The '+1' handles the case when max_dl is not a square
        for j in xrange(table_size):
            giant_steps[running_step] = exponent
            running_step = oEC.addEFp(grp, running_step, giant_step)
            exponent += max_search
            # print("Giant steps: ", giant_steps)

        self.ECtable = giant_steps
        return self.ECtable

    def make_Ftable(self, field, elt):
        # This function makes a multiplication table to aid in computing discrete
        # logarithims to speed up decryption of multiple messages encrypted with the
        # same public/private key

        baby_steps = {}
        gt = oEC.toTupleFp12(field.one())

        for j in xrange(2**16):
            gt = oEC.tmulFp12(self.Gt, gt, oEC.toTupleFp12(elt), self.Gamma)
            baby_steps[gt] = j + 1

        self.Ftable = baby_steps
        return self.Ftable

    def make_full_Ftable(self, field, elt, dltable = None):
        # This function makes a multiplication table to aid in computing discrete
        # logarithims to speed up decryption of multiple messages encrypted with the
        # same public/private key

        baby_steps = {}
        gt = oEC.toTupleFp12(field.one())

        for j in range((2**20)+1):
            gt = oEC.tmulFp12(self.Gt, gt, oEC.toTupleFp12(elt), self.Gamma)
            if not dltable == None:
                dltable.add_row(gt, j+1)
            else:
                baby_steps[gt] = j + 1

        self.fieldtable_full = baby_steps
        return self.fieldtable_full

    def make_full_ECtable(self, grp, point):
        # This function makes a multiplication table to aid in computing discrete
        # logarithims to speed up decryption of multiple messages encrypted with the
        # same public/private key

        # VT: I edited this to make it terminate at exactly the index it previously did.
        # However, I'd assume that the omission of the j=2**16 + 1 case makes no
        # difference so we should just let the range terminate at 2**16 as before.
        # But I'm reluctant to fiddle in case I mess up something important.

        baby_steps = {}
        pt = point
        for j in range(2**5 + 1):
            baby_steps[pt] = j + 1
            pt = oEC.addEFp(self.G, pt, point)

        self.ectable_full = baby_steps
        return self.ectable_full
    def log_full_field(self, dec):
        if self.use_dltable:
            ret = self.dltable.lookup(dec)
            if not ret == None:
                return int(ret)
            else:
                return "No Match"
        else:
            table = self.fieldtable_full
            if dec in table:
                return table[dec]
            return "No Match"
    def log_full_group(self, dec):
        table = self.ectable_full
        if dec in table:
            return table[dec]
        return "No Match"
