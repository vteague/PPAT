from __future__ import print_function
from multiprocessing import Queue
from threading import Thread
from multiprocessing import Process
from Crypto.Random import atfork

class BallotWorker(Process):
    def __init__(self, ballots, start_idx, end_idx, tally_queue, round_count, eliminated, sourcegrp, targetgrp, pubkey, secretkey, thread_id):
        super(BallotWorker, self).__init__()
        self.ballots = ballots
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.tally_queue = tally_queue
        self.round_count = round_count
        self.eliminated = eliminated
        self.targetgrp = targetgrp
        self.sourcegrp = sourcegrp
        self.pubkey = pubkey
        self.secretkey = secretkey
        self.thread_id = thread_id

    def run(self):
        atfork()
        tallies = []
        counter = 1
        for i in range(self.start_idx, self.end_idx):
            #print("Adding ballot: Thread", self.thread_id, " :", counter)
            self.ballots[i].add_to_tally(self.sourcegrp, self.targetgrp, self.pubkey, self.secretkey, tallies, self.eliminated, self.round_count)
            #print("Finished adding ballot: Thread", self.thread_id, " :", counter)
            counter = counter + 1
        self.tally_queue.put(tallies)
        self.tally_queue.put(self.thread_id)
