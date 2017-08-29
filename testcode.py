import sys
from ballots import Ballots
from cryptogroup import CryptoGroup

ballots = Ballots('samplevotes')
ballots.load()
group = CryptoGroup()
pk,sk = group.KeyGen()

ECtable = group.make_ECtable(group.G,group.EFpTupleToPoint(pk['g']))
Ftable = group.make_Ftable(group.Gt,pk['e'])
Ftable = group.make_full_Ftable(group.Gt,pk['e'])
print('Multiplication Tables Created')

M0 = 1
M1 = 1
C0_src = group.Enc_src(pk,M0)
C1_src = group.Enc_src(pk,M1)

#C0_tgt = group.Enc_tgt(pk,M0)
C1_tgt = group.Enc_tgt(pk,M1)
print('Enc_src Complete')

M0xM1= group.Multiply_src(pk,C0_src['C0'],C1_src['C1'])
M0plusM1=group.Add_src(pk,C0_src,C1_src)
M0plusM1Tgt=group.Add_tgt(pk,C1_tgt,C1_tgt)
print "test:", group.Dec_tgt(sk,pk,M0plusM1Tgt)

M0plusM1Tgt=group.Add_tgt(pk,M0xM1,C1_tgt)
neg = group.negate_src(C1_src)
test = group.Add_src(pk, C0_src, neg)
test2 = group.Add_src(pk, test, C1_src)
test3 = group.Multiply_src(pk,test['C0'],C1_src['C1'])
test4 = group.sim_switch(sk,pk,test3)
print "test4:", group.Dec_src(sk,pk,test4)
test5 = group.Add_src(pk, test4, C1_src)
print "test5:", group.Dec_src(sk,pk,test5)


print group.Dec_src(sk,pk,test2)

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

