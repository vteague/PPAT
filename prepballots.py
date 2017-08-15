import sys
from ballots import Ballots
from cryptogroup import CryptoGroup

ballots = Ballots('samplevotes')
ballots.load()
group = CryptoGroup()
pk,sk = group.KeyGen()

ECtable = group.make_ECtable(group.G,pk['g'])
Ftable = group.make_Ftable(group.Gt,pk['e'])
print('Multiplication Tables Created')

ballots.encrypt_prefs(group, pk)