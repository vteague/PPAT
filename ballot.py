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
from __future__ import print_function
import time
class Ballot:
    """
    Ballot Class that contains the preferences for a single ballot.

    Also used to hold the encrypted preferences and to
    add those preferences into the irv count."""
    def __init__(self):
        self.preferences = []
        self.encprefs = []

    def add_pref_row(self, row):
        """
        Adds a string from the vote file of the form 0 1 0 to this ballot. Converts the
        string into a list of integer preferences

        Args:
            row (string): string contain preference
        """
        prefrow = [int(e) if e.isdigit() else e for e in row.rstrip('\n').split(' ')]
        self.preferences.append(prefrow)
        return len(prefrow)

    def __str__(self):
        return self.preferences.__str__()

    def encrypt_prefs(self, group, pubkey):
        """
        Encrypts all the preferences (0 and 1's) using the group and public key"""
        for preflist in self.preferences:
            encprefrow = []
            for pref in preflist:
                encprefrow.append(group.Enc_src(pubkey, pref))
            self.encprefs.append(encprefrow)

    def add_to_tally(self, group, pubkey, secretkey, tallies, eliminated):
        """
        Adds the prefernecs into the irv tally. This is where the heart of the
        irv count is performed. Sums each row homomorphically, calculates row
        multiplier (0 or 1) to determine current preference, multiplies each
        possible preference value in each column with the respective multiplier,
        before finally summing the columns"""
        row_sums = []
        # sum rows
        for encprefrow in self.encprefs:
            row_sum = None
            for pref_counter in range(0, len(encprefrow)):
                if pref_counter not in eliminated:
                    if row_sum is None:
                        row_sum = encprefrow[pref_counter]
                    else:
                        row_sum = group.Add_src(pubkey, row_sum, encprefrow[pref_counter])
            row_sums.append(row_sum)
        
        # encryption of one and running total set to 0
        enc_one = group.Enc_src(pubkey, 1)
        running_total = group.Enc_src(pubkey, 0)
        enc_previous_row = enc_one

        #for each row calculate the multiplier by multipling it by 1-sum(previous multipliers)
        for row_counter in range(0, len(row_sums)):
            row_adjust = group.Multiply_src(pubkey, enc_previous_row['C0'], row_sums[row_counter]['C1'])
            row_sums[row_counter] = group.sim_switch(secretkey, pubkey, row_adjust)
            running_total = group.Add_src(pubkey, running_total, row_sums[row_counter])
            enc_previous_row = group.Add_src(pubkey, enc_one, group.negate_src(running_total))
        
        # apply row multiplier (1 or 0) to select only current preference
        firstrow = self.encprefs[0]
        column_tallies = []
        for pref_counter in range(0, len(firstrow)):
            #Skip if we have eliminated this index - appending None
            if pref_counter not in eliminated:
                column_tallies.append(group.Multiply_src(pubkey, firstrow[pref_counter]['C0'], row_sums[0]['C1']))
            else:
                column_tallies.append(None)

        # sum the tallies in each column (should be a list with a 1 in the column of current pref)
        for row_counter in range(1, len(row_sums)):
            for pref_counter in range(0, len(firstrow)):
                #ignore if we have eliminated
                if pref_counter not in eliminated:
                    mult_pref = group.Multiply_src(pubkey, self.encprefs[row_counter][pref_counter]['C0'], row_sums[row_counter]['C1'])
                    column_tallies[pref_counter] = group.Add_tgt(pubkey, column_tallies[pref_counter], mult_pref)
        # check if this is the first ballot to add to the overall tallies
        tallylength = len(tallies)
        if tallylength == 0:
            #if it is just append out column tallies
            for col_counter in range(0, len(column_tallies)):
                tallies.append(column_tallies[col_counter])
        else:
            # otherwise homomorphically sum with overall tally for column
            for col_counter in range(0, len(column_tallies)):
                # skip column if eliminated
                if col_counter not in eliminated:
                    tallies[col_counter] = group.Add_tgt(pubkey, tallies[col_counter], column_tallies[col_counter])

