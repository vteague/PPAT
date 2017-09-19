"""
# Copyright 2017 Olivier Pereira & Chris Culnane
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
from __future__ import print_function
import mathTools.otosEC as oEC
import gmpy2 as gmpy
from dltable import DLTable

class EFp12Table(DLTable):
    """
    Class to build and lookup a baby step giant step Fp12 Discrete Log table
    """

    def __init__(self, table, group):
        super(EFp12Table, self).__init__(table, group)
        self.max_extract_steps = 0

    def build(self, base, max_dl=2 ** 32, max_search=2 ** 12):
        """This function makes a multiplication table to aid in computing discrete
        logarithms to speed up decryption of multiple messages encrypted with the
        same public/private key
        :param base is a tuple
        :param max_dl is the highest value that the DL can take
        :param max_search is the maximum number of steps that we agree to make during DL extraction
        """
        self.max_extract_steps = max_search * 2
        giant_step = oEC.squareAndMultiplyFp12(self.group.Gt, oEC.toTupleFp12(base), max_search,
                                               oEC.tmulFp12, oEC.tsqrtFp12, self.group.Gamma)
        # Store the giant steps. Keys are truncated x coordinate of points, values are exponents
        table_size = max_dl / max_search + 1
        # Size of the giant steps (on the curve)
        # Counter of current value of the exponent
        exponent = max_search
        # Position of that counter on the curve
        running_step = giant_step
        # j performs as many steps as needed.
        # The '+1' handles the case when max_dl is not a square
        for j in xrange(table_size):
            lsb_running_step = int(gmpy.t_mod_2exp(running_step[0], 128))
            self.table.add_row(lsb_running_step, exponent)
            running_step = oEC.tmulFp12(self.group.Gt, running_step, giant_step, self.group.Gamma)
            exponent += max_search


    def extract(self, a, b):
        """Extracts the discrete log of b in base a in Group, which must be an EFp.
        Assumes table contains precomputed values for base a

        :param a: base element as a tuple
        :param b: tuple from which DL must be extracted

        :return: x: b == a**x
        """
        i = 0
        lsb_running_step = int(gmpy.t_mod_2exp(b[0], 128))
        lookup = self.table.lookup(lsb_running_step)
        while lookup is None and i <= self.max_extract_steps:
            b = oEC.tmulFp12(self.group.Gt, a, b, self.group.Gamma)
            i += 1
            lsb_running_step = int(gmpy.t_mod_2exp(b[0], 128))
            lookup = self.table.lookup(lsb_running_step)
        if lookup is None:
            return lookup
        else:
            adl = oEC.squareAndMultiplyFp12(self.group.Gt, a, lookup,
                                            oEC.tmulFp12, oEC.tsqrtFp12, self.group.Gamma)
            assert adl == b
            return lookup - i
