import sys
import os
import time
from ballots import Ballots
from irv import IRV
from dltable import DLTable
from cryptogroup import CryptoGroup

ballots = Ballots('samplevotes')
ballots.load()
group = CryptoGroup()
#pk,sk = group.KeyGen()
#group.save_public_key(pk,'pubkey.json')
#group.saveSecretKey(sk,'secretkey.json')
key = {}
key['pk']=group.loadPublicKey('pubkey.json')
key['sk']=group.loadSecretKey('secretkey.json')
pk,sk = group.KeyGen(key)

if os.path.exists('sortedtable.tbl'):
    group.load_dltable('sortedtable.tbl', 72)
else:
    print("Creating DL Field Table")
    # Create a DLTable with the appropraite line length and file name
    dltable = DLTable(group, 'unsortedtable.tbl', 72)

    # Open the table for writing - i.e. we are creating a new table
    dltable.open(for_writing=True)

    print("Starting to build table")
    group.make_full_Ftable(group.Gt, pk['e'], dltable)

    # Close the table to ensure buffers are flushed
    dltable.close()
    print("Finished building table")

    # Call the sort method, with an output file name - performs Unix Sort
    print("Starting sort")
    dltable.sort("sortedtable.tbl")
    print("Finished sort")

    # Create a new DLTable pointing to the sorted table
    group.load_dltable('sortedtable.tbl', 72)

#group.make_full_Ftable(group.Gt, pk['e'],"mytable.tbl")
#ECtable = group.make_ECtable(group.G, pk['g'])
#Ftable = group.make_full_Ftable(group.Gt, pk['e'])

ballots.encrypt_prefs(group,pk)

irv = IRV(group,pk,sk)
start=time.time()
irv.perform_count(ballots)
end=time.time()
print("Count:",end-start)
#tallies = ballots.run_count(group,pk,sk)

#tallycounter=1
#for tally in tallies:
#    print "candidate", tallycounter, group.Dec_tgt(sk, pk, tally)
#    tallycounter = tallycounter + 1
