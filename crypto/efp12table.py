from __future__ import print_function
from sys import getsizeof
import mathTools.otosEC as oEC
import gmpy2 as gmpy
from dltable import DLTable

class EFp12Table(DLTable):

    def __init__(self, table, group):
        super(EFp12Table, self).__init__(table, group)
        self.giant_step_size=0;

    def build(self, base, max_dl=2 ** 32, max_search=2 ** 12):
        """This function makes a multiplication table to aid in computing discrete
        logarithms to speed up decryption of multiple messages encrypted with the
        same public/private key
        :param base is a tuple
        :param max_dl is the highest value that the DL can take
        :param max_search is the maximum number of steps that we agree to make during DL extraction
        """
        self.giant_step_size = max_search * 2
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
        while lookup is None and i <= self.giant_step_size:
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
