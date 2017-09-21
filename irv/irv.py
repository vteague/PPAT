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

class IRV:
    """ IRV Counter
    for use with encryption switching """

    def __init__(self, sourcegrp, targetgrp, pubkey, secretkey):
        self.sourcegrp = sourcegrp
        self.targetgrp = targetgrp
        self.pubkey = pubkey
        self.secretkey = secretkey

    def perform_count(self, ballots):
        """start an IRV count on the specified ballots"""
        for round_count in range(0, ballots.candcount):
            print("Starting tally round", round_count)
            round_result = self.perform_round(ballots, round_count)
            if round_result['stop']:
                print("Result Found")
                return
            if round_result['removeidx'] != None:
                ballots.eliminate(round_result['removeidx'])

    def perform_round(self, ballots, irvround):
        """Perform a single round of IRV counting"""
        tallies = ballots.run_threaded_count(self.sourcegrp, self.targetgrp, self.pubkey, self.secretkey, irvround)
        dec_tally = []
        min_indexes = []
        min_idx = ballots.ballotcount
        candidateindex = 0
        total_remaining_votes = 0
        max_total = 0
        for tally in tallies:
            if tally is not None:
                decrypted_total = self.targetgrp.decrypt(self.secretkey, self.pubkey, tally)
                total_remaining_votes = total_remaining_votes + decrypted_total
                if decrypted_total > max_total:
                    max_total = decrypted_total
                dec_tally.append(decrypted_total)
                if decrypted_total < min_idx:
                    min_idx = decrypted_total
                    del min_indexes[:]
                    min_indexes.append(candidateindex)
            else:
                dec_tally.append(None)
            candidateindex = candidateindex + 1
        print("Finished tally round")
        print("Tally:", dec_tally.__str__())
        print("Candidates with mins:", min_indexes.__str__())
        remove_cand_idx = -1
        if len(min_indexes) > 1:
            print("Selecting candidate to remove at random")
            # TODO make this random
            remove_cand_idx = min_indexes[0]
        elif len(min_indexes) == 0:
            remove_cand_idx = None
        else:
            remove_cand_idx = min_indexes[0]
        if max_total > total_remaining_votes / 2:
            stop = True
        else:
            stop = False
        return {'removeidx':remove_cand_idx, 'tally':dec_tally, 'stop':stop}
