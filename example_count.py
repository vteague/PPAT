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

ballots = Ballots('./data/samplevotes')
ballots.load()

field = CryptoField()
memtable = MemTable()
table = EFpTable(memtable, field.G)
sourcegrp = SourceGroup(field, table)
key = {}
key['pk'] = sourcegrp.load_public_key('./data/pubkey.json')
key['sk'] = sourcegrp.load_secret_key('./data/secretkey.json')
pk, sk = sourcegrp.key_gen(key=key)
table.build(pk['g'], 2 ** 32, 2 ** 12)

memtable = MemTable()
table = EFp12Table(memtable, field)
targetgrp = TargetGroup(field, table)
table.make_full_Ftable(pk['e'])

ballots.encrypt_prefs(sourcegrp,pk)

irv = IRV(sourcegrp, targetgrp, pk, sk)
start=time.time()
irv.perform_count(ballots)
end=time.time()
print("Count:",end-start)

