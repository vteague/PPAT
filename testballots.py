import sys
from ballots import Ballots
from irv import IRV
from cryptogroup import CryptoGroup

ballots = Ballots('samplevotes')
ballots.load()
group = CryptoGroup()
pk,sk = group.KeyGen()

ECtable = group.make_ECtable(group.G, pk['g'])
Ftable = group.make_full_Ftable(group.Gt, pk['e'])

ballots.encrypt_prefs(group,pk)

irv = IRV(group,pk,sk)
irv.perform_count(ballots)

#tallies = ballots.run_count(group,pk,sk)

#tallycounter=1
#for tally in tallies:
#    print "candidate", tallycounter, group.Dec_tgt(sk, pk, tally)
#    tallycounter = tallycounter + 1
