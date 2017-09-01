import sys
import os
import time
from irv.ballots import Ballots
from irv.irv import IRV
from crypto.dltable import DLTable
from crypto.cryptogroup import CryptoGroup

ballots = Ballots('./data/samplevotes')
ballots.load()
group = CryptoGroup()
#pk,sk = group.KeyGen()
#group.save_public_key(pk,'./data/pubkey.json')
#group.saveSecretKey(sk,'./data/secretkey.json')
key = {}
key['pk']=group.loadPublicKey('./data/pubkey.json')
key['sk']=group.loadSecretKey('./data/secretkey.json')
pk,sk = group.KeyGen(key)

if os.path.exists('./data/sortedtable.tbl'):
    group.load_dltable('./data/sortedtable.tbl', 72)
else:
    print("Creating DL Field Table")
    # Create a DLTable with the appropraite line length and file name
    dltable = DLTable(group, './data/unsortedtable.tbl', 72)

    # Open the table for writing - i.e. we are creating a new table
    dltable.open(for_writing=True)

    print("Starting to build table")
    group.make_full_Ftable(group.Gt, pk['e'], dltable)

    # Close the table to ensure buffers are flushed
    dltable.close()
    print("Finished building table")

    # Call the sort method, with an output file name - performs Unix Sort
    print("Starting sort")
    dltable.sort("./data/sortedtable.tbl")
    print("Finished sort")

    # Create a new DLTable pointing to the sorted table
    group.load_dltable('./data/sortedtable.tbl', 72)

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
