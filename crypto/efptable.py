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

class EFpTable(DLTable):
    """
    Class to build and lookup a baby step giant step Fp Discrete Log table
    """
    def __init__(self, table, group):
        super(EFpTable, self).__init__(table, group)
        self.max_extract_steps = 0

    def build(self, base, max_dl=2 ** 32, max_search=2 ** 12):
        """This function makes a multiplication table to aid in computing discrete
        logarithms to speed up decryption of multiple messages encrypted with the
        same public/private key
        :param base is a tuple
        :param max_dl is the highest value that the DL can take
        :param max_search is the maximum number of steps that we agree to make during DL extraction
        """
        # Store the giant steps. Keys are truncated x coordinate of points, values are exponents
        self.max_extract_steps = max_search * 2

        table_size = max_dl / max_search + 1
        # Size of the giant steps (on the curve)
        giant_step = oEC.mulECP(self.group, base, max_search)
        # Counter of current value of the exponent
        exponent = max_search
        # Position of that counter on the curve
        running_step = giant_step
        # j performs as many steps as needed.
        # The '+1' handles the case when max_dl is not a square
        for j in xrange(table_size):
            lsb_running_step = int(gmpy.t_mod_2exp(running_step[0], 128))
            self.table.add_row(lsb_running_step, exponent)
            running_step = oEC.addEFp(self.group, running_step, giant_step)
            exponent += max_search



    def extract(self, a, b):
        """Extracts the discrete log of b in base a in Group, which must be an EFp.
        Assumes table contains precomputed values for base a

        :param a: base element as a tuple
        :param b: tuple from which DL must be extracted

        :return: x: b == a**x
        """
        i = 0
        lookup = self.table.lookup(int(gmpy.t_mod_2exp(b[0], 128)))
        while lookup is None and i <= self.max_extract_steps:
            b = oEC.addEFp(self.group, a, b)
            i += 1
            lookup = self.table.lookup(int(gmpy.t_mod_2exp(b[0], 128)))
        if lookup is None:
            return lookup
        else:
            adl = oEC.mulECP(self.group, a, lookup, sq=False)
            assert adl == b
            return lookup - i
