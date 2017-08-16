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
from ballot import Ballot

class Ballots:
    """Wrapper for ballot objects
    Can be used to load ballots from a file and to then
    perform an irv count

    Args:
        ballot_file_path (string): path to ballot file"""
    def __init__(self, ballot_file_path):
        self.file_path = ballot_file_path
        self.ballots = []
        self.candcount = 0
        self.eliminated = []

    @property
    def ballotcount(self):
        """
        Gets the number of ballots
        """
        return len(self.ballots)

    def eliminate(self, cand_idx):
        """Eliminate a candidate from the count. This does not
        remove the candidate column, it just treats it a None

        Args:
            cand_idx (int): candidate index to eliminate"""
        self.eliminated.append(cand_idx)

    def load(self):
        """Loads ballots from the file specified during construction.

        The ballot file should be preference string of the form
        0 0 1, separated by spaces, with each preference on a new line.
        For example, a three candidate, three preference entry would look like:

        0 0 1\n
        1 0 0\n
        0 1 0

        Each ballot is sepasrated by a blank line
        """
        ballotfile = open(self.file_path, "r")
        currentballot = Ballot()
        row_added = False

        for line in ballotfile:
            if line == "\n":
                self.ballots.append(currentballot)
                currentballot = Ballot()
                row_added = False
                # new ballot
            else:
                prefcount = currentballot.add_pref_row(line)
                if prefcount > self.candcount:
                    self.candcount = prefcount
                row_added = True
        if row_added:
            self.ballots.append(currentballot)
        print(len(self.ballots), 'ballots added with', self.candcount, 'candidates.')

    def encrypt_prefs(self, group, pubkey):
        """Encrypt the preferences in all the ballots

        Args:
            group (CryptoGroup): crypto group to use for encryption
            pubkey: Public key object
        """
        for ballot in self.ballots:
            ballot.encrypt_prefs(group, pubkey)

    def run_count(self, group, pubkey, secretkey):
        """Run the irv count on the ballots and return the tallies
        Args:
            group (CryptoGroup): crypto group to use for encryption
            pubkey: Public key object
            secretkey: Secret key object
        Returns:
            tallies list tally for each candidate
        """
        tallies = []
        counter = 1
        for ballot in self.ballots:
            print("Adding ballot:", counter)
            ballot.add_to_tally(group, pubkey, secretkey, tallies, self.eliminated)
            print("Finished adding ballot:", counter)
            counter = counter + 1
        return tallies
       