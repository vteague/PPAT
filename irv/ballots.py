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
import multiprocessing
from multiprocessing import Queue
from ballot import Ballot
from ballotworker import BallotWorker
from encryptionworker import EncryptionWorker
"""def pcount(ballot,sourcegrp,targetgrp,pubkey,secretkey,eliminated):
    return ballot.add_to_tally(sourcegrp, targetgrp, pubkey, secretkey, [], eliminated)
"""
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
        ballot_count=0
        total_count=0
        for ballot in self.ballots:
            ballot.encrypt_prefs(group, pubkey)
            ballot_count = ballot_count + 1
            total_count = total_count + 1
            if ballot_count >= 100:
                ballot_count = 0
                print("Encrypted:", total_count, "ballots")
    
    def encrypt_prefs_threaded(self, group, pubkey):
        """Encrypt the preferences in all the ballots

        Args:
            group (CryptoGroup): crypto group to use for encryption
            pubkey: Public key object
        """
        core_count = multiprocessing.cpu_count()
        thread_share = len(self.ballots) / core_count
        ballot_queue = Queue()
        workers = []
        worker_ids = []
        for x in range(core_count):
            start_idx = x * thread_share
            if x == (core_count - 1):
                end_idx = len(self.ballots)
            else:
                end_idx = start_idx + thread_share
            worker = EncryptionWorker(self.ballots, start_idx, end_idx, ballot_queue, group, pubkey, x)
            worker_ids.append(x)
            worker.start()
            workers.append(worker)
        waiting = True
        self.ballots = []
        while waiting:
            entry = ballot_queue.get()
            if entry in worker_ids:
                worker_ids.remove(entry)
            else:
                self.ballots.extend(entry)
            if len(worker_ids) == 0:
                waiting = False
        for worker in workers:
            worker.join()
        #self.ballots = ballot_queue.get()
        #while not ballot_queue.empty():
        #    self.ballots.extend(ballot_queue.get())

    def run_threaded_count(self, sourcegrp, targetgrp, pubkey, secretkey, round_count):
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
        core_count = multiprocessing.cpu_count()
        thread_share = len(self.ballots) / core_count
        tally_queue = Queue()
        workers = []
        worker_ids = []
        for x in range(core_count):
            start_idx = x * thread_share
            if x == (core_count - 1):
                end_idx = len(self.ballots)
            else:
                end_idx = start_idx + thread_share
            worker = BallotWorker(self.ballots, start_idx, end_idx, tally_queue, round_count, self.eliminated, sourcegrp, targetgrp, pubkey, secretkey, x)
            worker.start()
            worker_ids.append(x)
            workers.append(worker)
        waiting = True
        while waiting:
            entry = tally_queue.get()
            if entry in worker_ids:
                worker_ids.remove(entry)
            else:
                if len(tallies) == 0:
                    tallies = entry
                else:
                    for col_counter in range(0, len(entry)):
                        if tallies[col_counter] is not None and entry[col_counter] is not None:
                            tallies[col_counter] = targetgrp.add(pubkey, tallies[col_counter], entry[col_counter])
            if len(worker_ids) == 0:
                waiting = False
        for worker in workers:
            worker.join()
        return tallies
    def run_count(self, sourcegrp, targetgrp, pubkey, secretkey, round_count):
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
        """pool = ThreadPool()
        chunksplit = max(1,len(self.ballots)/8)
        print("ChunkSize:", chunksplit)
        results = pool.map(partial(pcount, sourcegrp=sourcegrp,targetgrp=targetgrp, pubkey=pubkey,secretkey=secretkey, eliminated=self.eliminated), self.ballots,chunksize=chunksplit)
        for column_tallies in results:
            tallylength = len(tallies)
            if tallylength == 0:
                #if it is just append out column tallies
                for col_counter in range(0, len(column_tallies)):
                    tallies.append(column_tallies[col_counter])
            else:
                # otherwise homomorphically sum with overall tally for column
                for col_counter in range(0, len(column_tallies)):
                    # skip column if eliminated
                    if col_counter not in self.eliminated:
                        tallies[col_counter] = targetgrp.add(pubkey,
                                                             tallies[col_counter],
                                                             column_tallies[col_counter])
        """
        for ballot in self.ballots:
            print("Adding ballot:", counter)
            ballot.add_to_tally(sourcegrp, targetgrp, pubkey, secretkey, tallies, self.eliminated, round_count)
            print("Finished adding ballot:", counter)
            counter = counter + 1
        return tallies

    