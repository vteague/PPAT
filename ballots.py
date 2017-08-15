from ballot import Ballot

class Ballots:
    def __init__(self, ballotFilePath):
        self.file_path = ballotFilePath
        self.ballots = []
        self.candcount = 0
        self.eliminated = []
    def eliminate(self, cand_idx):
        self.eliminated.append(cand_idx)

    def load(self):
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
        print(len(self.ballots), ' ballots added with ', self.candcount, ' candidates.')

    def encrypt_prefs(self, group, pubkey):
        for ballot in self.ballots:
            ballot.encrypt_prefs(group, pubkey)

    def run_count(self, group, pubkey, sk):
        tallies = []
        counter = 1
        for ballot in self.ballots:
            print("Adding ballot: ", counter)
            ballot.add_to_tally(group, pubkey, sk, tallies, self.eliminated)
            print("Finished adding ballot: ", counter)
            counter = counter + 1
        return tallies
       