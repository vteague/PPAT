import sys
from ballot import Ballot


class Ballots:
    def __init__(self, ballotFilePath):
        self.file_path = ballotFilePath
        self.ballots = []
        self.candcount =0 

    def load(self):
        file = open(self.file_path, "r")
        currentballot = Ballot()
        row_added = False

        for line in file:
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
        print(len(self.ballots), ' ballots added with ', self.candcount, ' candidates.')

    def encrypt_prefs(self, group, pubkey):
        for ballot in self.ballots:
            ballot.encrypt_prefs(group, pubkey)

    def run_count(self, group, pubkey):
        tallies = []
        
        for ballot in self.ballots:
            ballot.add_to_tally(group, pubkey, tallies)
        return tallies
       