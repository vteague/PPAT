from __future__ import print_function
from sys import getsizeof
import mathTools.otosEC as oEC
import gmpy2 as gmpy
from dltable import DLTable

class EFp12Table(DLTable):

    def __init__(self, table, group):
        super(EFp12Table, self).__init__(table, group)

    def make_full_Ftable(self, elt):
        # This function makes a multiplication table to aid in computing discrete
        # logarithims to speed up decryption of multiple messages encrypted with the
        # same public/private key

        gt = oEC.toTupleFp12(self.group.Gt.one())
        self.table.add_row(gt, 0)
        for j in range((2**10)+1):
            gt = oEC.tmulFp12(self.group.Gt, gt, oEC.toTupleFp12(elt), self.group.Gamma)
            self.table.add_row(gt, j+1)
            
    def extract_from_full(self,elt):
        return self.table.lookup(elt)
        


    def build(self, base, max_dl=2 ** 32, max_search=2 ** 12):
        """This function makes a multiplication table to aid in computing discrete
        logarithms to speed up decryption of multiple messages encrypted with the
        same public/private key
        :param base is a tuple
        :param max_dl is the highest value that the DL can take
        :param max_search is the maximum number of steps that we agree to make during DL extraction
        """

        giant_step =base * max_search
        # Store the giant steps. Keys are truncated x coordinate of points, values are exponents
        print(giant_step)
        giant_steps = {}
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
            #giant_steps[lsb_running_step] = exponent
            running_step = oEC.addEFp(self.group, running_step, giant_step)
            exponent += max_search
            # print("Giant steps: ", giant_steps)
        #return giant_steps


    def extract(self, a, b):
        """Extracts the discrete log of b in base a in Group, which must be an EFp.
        Assumes table contains precomputed values for base a

        :param a: base element as a tuple
        :param b: tuple from which DL must be extracted

        :return: x: b == a**x
        """
        i = 0
        lookup = self.table.lookup(gmpy.t_mod_2exp(b[0], 128))
        while lookup is None:
            b = oEC.addEFp(self.group, a, b)
            i += 1
            lookup = self.table.lookup(gmpy.t_mod_2exp(b[0], 128))
        return lookup - i
