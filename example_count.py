from __future__ import print_function
import sys
import os
import time
from irv.ballots import Ballots
from irv.irv import IRV
from crypto.cryptofield import CryptoField
from crypto.memtable import MemTable
from crypto.efptable import EFpTable
from crypto.efp12table import EFp12Table
from crypto.sourcegroup import SourceGroup
from crypto.targetgroup import TargetGroup

ballots = Ballots('./data/auburn.txt')
ballots.load()

field = CryptoField()
memtable = MemTable()
table = EFpTable(memtable, field.G)
sourcegrp = SourceGroup(field, table)
key = {}
key['pk'] = sourcegrp.load_public_key('./data/pubkey2.json')
key['sk'] = sourcegrp.load_secret_key('./data/secretkey2.json')
pk, sk = sourcegrp.key_gen(key=key)
print("Starting to build source table")
start=time.time()
table.build(pk['g'], 2 ** 32, 2 ** 12)
end=time.time()
print("Finished building source table, took:", end-start)
memtable = MemTable()
table = EFp12Table(memtable, field)
targetgrp = TargetGroup(field, table)
print("Starting to build target table")
table.build(pk['e'], 2 ** 32, 2 ** 12)
print("Finished building target table, took:", end-start)

print("Starting ballot encryption")
start=time.time()
ballots.encrypt_prefs_threaded(sourcegrp,pk)
end=time.time()
print("Finished ballot encryption, took:", end-start)

irv = IRV(sourcegrp, targetgrp, pk, sk)
start=time.time()
irv.perform_count(ballots)
end=time.time()
print("Count:",end-start)

