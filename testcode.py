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

M0 = 0
M1 = 5
C0_src = group.Enc_src(pk,M0)
C1_src = group.Enc_src(pk,M1)

#C0_tgt = group.Enc_tgt(pk,M0)
#C1_tgt = group.Enc_tgt(pk,M1)
print('Enc_src Complete')

M0xM1= group.Multiply_src(pk,C0_src['C0'],C1_src['C1'])
M0plusM1=group.Add_src(pk,C0_src,C1_src)
#M0plusM1Tgt=group.Add_tgt(pk,C0_tgt,C1_tgt)

print group.Dec_src(sk,pk,C0_src,ECtable)
print group.Dec_src(sk,pk,C1_src,ECtable)
print group.Dec_src(sk,pk,M0plusM1,ECtable)
#print group.Dec_tgt(sk, pk, C0_tgt,Ftable)
#print group.Dec_tgt(sk, pk, C1_tgt,Ftable)
print group.Dec_tgt(sk, pk, M0xM1,Ftable)
print group.Dec_tgt(sk, pk, M0plusM1Tgt,Ftable)

#simulate multiple blinding factors
blinding_factors = []
for bfcount in range(0,3):
    blinding_factors.append(group.generate_blinding_factor(pk))

blinded_cipher = group.blind_pair_cipher(pk,M0xM1,blinding_factors)
switched_cipher = group.switch(sk,pk,blinded_cipher['C'],blinded_cipher['bfEC'],Ftable)

print "Switched:", group.Dec_src(sk,pk,switched_cipher,ECtable)

