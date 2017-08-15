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

    def add_to_tally(self, group, pubkey, sk, tallies, eliminated):

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

        # encryption of one
        enc_one = group.Enc_src(pubkey, 1)
        running_total = group.Enc_src(pubkey, 0)
        enc_previous_row = enc_one
        for row_counter in range(0, len(row_sums)):
            row_adjust = group.Multiply_src(pubkey, enc_previous_row['C0'], row_sums[row_counter]['C1'])
            row_sums[row_counter] = group.sim_switch(sk, pubkey, row_adjust)
            running_total = group.Add_src(pubkey, running_total, row_sums[row_counter])
            enc_previous_row = group.Add_src(pubkey, enc_one, group.negate_src(running_total))
            
        # sum the columns apply row sum (1 or 0) to select only current preference
        firstrow = self.encprefs[0]
        column_tallies = []
        for pref_counter in range(0, len(firstrow)):
            if pref_counter not in eliminated:
                column_tallies.append(group.Multiply_src(pubkey, firstrow[pref_counter]['C0'], row_sums[0]['C1']))
            else:
                column_tallies.append(None)
        for row_counter in range(1, len(row_sums)):
            for pref_counter in range(0, len(firstrow)):
                if pref_counter not in eliminated:
                    mult_pref = group.Multiply_src(pubkey, self.encprefs[row_counter][pref_counter]['C0'], row_sums[row_counter]['C1'])
                    column_tallies[pref_counter] = group.Add_tgt(pubkey, column_tallies[pref_counter], mult_pref)
        tallylength = len(tallies)
        if tallylength == 0:
            for col_counter in range(0, len(column_tallies)):
                tallies.append(column_tallies[col_counter])
        else:
            for col_counter in range(0, len(column_tallies)):
                if col_counter not in eliminated:
                    tallies[col_counter] = group.Add_tgt(pubkey, tallies[col_counter], column_tallies[col_counter])
