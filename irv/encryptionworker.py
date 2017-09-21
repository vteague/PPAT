from __future__ import print_function
from multiprocessing import Queue
from threading import Thread
from multiprocessing import Process
from Crypto.Random import atfork

class EncryptionWorker(Process):
    def __init__(self, ballots, start_idx, end_idx, ballot_queue, sourcegrp, pubkey, thread_id):
        super(EncryptionWorker, self).__init__()
        self.ballots = ballots
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.ballot_queue = ballot_queue
        self.sourcegrp = sourcegrp
        self.pubkey = pubkey
        self.thread_id = thread_id

    def run(self):
        atfork()
        counter = 0
        total = 0
        ballotrets = []
        for i in range(self.start_idx, self.end_idx):
            self.ballots[i].encrypt_prefs(self.sourcegrp, self.pubkey)
            ballotrets.append(self.ballots[i])
            #print("Finished adding ballot: Thread", self.thread_id, " :", counter)
            counter = counter + 1
            total = total + 1
            if counter >= 100:
                print("Thread", self.thread_id, "encrypted", total)
                counter = 0
        self.ballot_queue.put(ballotrets)
        self.ballot_queue.put(self.thread_id)

