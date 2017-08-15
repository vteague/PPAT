import sys
from ballots import Ballots
from cryptogroup import CryptoGroup

ballots = Ballots('samplevotes')
ballots.load()
group = CryptoGroup()
pk,sk = group.KeyGen()

ECtable = group.make_ECtable(group.G,pk['g'])
Ftable = group.make_Ftable(group.Gt,pk['e'])

ballots.encrypt_prefs(group,pk)
tallies = ballots.run_count(group,pk)