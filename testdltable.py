import sys
from ballots import Ballots
from dltable import DLTable
from cryptogroup import CryptoGroup
import mathTools.otosEC as oEC


def createval(idx, group, elt):
    """
    Utlity method for creating test values
    """
    gt = oEC.toTupleFp12(group.Gt.one())
    for j in range(idx):
        gt = oEC.tmulFp12(group.Gt, gt, elt, group.Gamma)
    return gt

# Calculate line length based on how big an integer value we wish to store
linelength = DLTable.calculate_line_length(20)
print("LineLength:", linelength)
# Create group and keys
group = CryptoGroup()
key = {}
key['pk']=group.loadPublicKey('pubkey.json')
key['sk']=group.loadSecretKey('secretkey.json')
pk,sk = group.KeyGen(key)

print("Created group")

# Create a DLTable with the appropraite line length and file name
dltable = DLTable(group, 'mytesttable.tbl', linelength)

# Open the table for writing - i.e. we are creating a new table
dltable.open(for_writing=True)

print("Starting to build table")
# Iterate through values as before
elt = oEC.toTupleFp12(pk['e'])
gt = oEC.toTupleFp12(group.Gt.one())
for j in range((2**20)+1):
    gt = oEC.tmulFp12(group.Gt, gt, elt, group.Gamma)
    # add the value to the table
    dltable.add_row(gt, j+1)

# Close the table to ensure buffers are flushed
dltable.close()
print("Finished building table")

# Call the sort method, with an output file name - performs Unix Sort
print("Starting sort")
dltable.sort("sortedtest.tbl")
print("Finished sort")

# Create a new DLTable pointing to the sorted table
dltable = DLTable(group, 'sortedtest.tbl', linelength)

# Open the table, by default this is read only
dltable.open()

print("Starting test")
# Perform test
found = 0
notfound = 0
for j in range((2**10)+1):
    lookup = createval(j+1, group, elt)
    ret = dltable.lookup(lookup)
    if ret != None:
        ret = int(ret)
        if ret != j+1:
            print("Not found", j+1, dltable._create_hash_string(lookup))
            notfound = notfound + 1
        else:
            found = found + 1
    else:
        notfound = notfound + 1
        print("Not found", j+1, dltable._create_hash_string(lookup))

print("Finished test")
print("found:", found)
print("notfound:", notfound)
print("Avg Lookup Time:", dltable.get_lookup_time())
dltable.close()
