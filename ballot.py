import sys
from copy import deepcopy
from copy import copy
from cryptogroup import CryptoGroup
class Ballot:
    def __init__(self):
        self.preferences = []
        self.encprefs = [];

    def add_pref_row(self, row):
        prefrow = [int(e) if e.isdigit() else e for e in row.rstrip('\n').split(' ')]
        self.preferences.append(prefrow)
        return len(prefrow)

    def __str__(self):
        return self.preferences.__str__()

    def encrypt_prefs(self, group, pubkey):
        for preflist in self.preferences:
            encprefrow = []
            for pref in preflist:
                encprefrow.append(group.Enc_src(pubkey, pref))
            self.encprefs.append(encprefrow)

    def add_to_tally(self, group, pubkey, tallies):

        row_sums = []
        # sum rows
        for encprefrow in self.encprefs:
            row_sum = encprefrow[0]
            for pref_counter in range(1, len(encprefrow)):
                row_sum = group.Add_src(pubkey, row_sum, encprefrow[pref_counter])
            row_sums.append(row_sum)
        print(len(row_sums))
        # encryption of one
        enc_running_total = row_sums[0]
        for row_counter in range(1, len(row_sums)):
            row_sums[row_counter] = group.Add_src(pubkey, row_sums[row_counter], group.negate_src(enc_running_total))
            enc_running_total = group.Add_src(pubkey, row_sums[row_counter], enc_running_total)

        # sum the columns apply row sum (1 or 0) to select only current preference
        firstrow = self.encprefs[0]
        column_tallies = []
        for pref_counter in range(0, len(firstrow)):
            print(firstrow[pref_counter])
            print(row_sums[0])
            column_tallies[pref_counter] = group.Multiply_src(pubkey, firstrow[pref_counter], row_sums[0])
        for row_counter in range(1, len(row_sums)):
            for pref_counter in range(0, len(firstrow)):
                mult_pref = group.Multiply_src(pubkey, self.encprefs[row_counter][pref_counter], row_sums[row_counter])
                column_tallies[pref_counter] = group.Add_tgt(pubkey, column_tallies[pref_counter], mult_pref)
        if len(tallies) == 0:
            for col_counter in range(0, len(column_tallies)):
                tallies[0] = column_tallies[col_counter]
        else:
            for col_counter in range(0, len(column_tallies)):
                group.Add_tgt(pubkey, tallies[0], column_tallies[col_counter])
