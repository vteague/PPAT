import sys
import time
from ballots import Ballots
from irv import IRV
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

group.make_offline_Ftable(group.Gt, pk['e'],"mytable.tbl")
#ECtable = group.make_ECtable(group.G, pk['g'])
Ftable = group.make_full_Ftable(group.Gt, pk['e'])

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
