import sys
from cryptogroup import CryptoGroup

group = CryptoGroup()
pk,sk = group.KeyGen()

ECtable = group.make_ECtable(group.G,group.EFpTupleToPoint(pk['g']))
#ECtable = group.make_ECtable(group.G,pk['g'])
Ftable = group.make_Ftable(group.Gt,pk['e'])
group.make_full_ECtable(group.G,pk['g'])
group.make_full_Ftable(group.Gt,pk['e'])
print('Multiplication Tables Created')

M0 = 2
M1 = 3
C0_src = group.Enc_src(pk,M0)
C1_src = group.Enc_src(pk,M1)
print('Enc_src Complete')

C0_tgt = group.Enc_tgt(pk,M0)
C1_tgt = group.Enc_tgt(pk,M1)
print('Enc_tgt Complete')

M0plusM1=group.Add_src(pk,C0_src,C1_src)
M0plusM1Tgt=group.Add_tgt(pk,C0_tgt,C1_tgt)
print('Cipher text addition Complete')

M0xM1 = group.Multiply_src(pk,C0_src['C0'],C1_src['C1'])
print('Cipher text multiplication Complete')

M0xM1plusM1 = group.Add_tgt(pk,M0xM1,C1_tgt)
print('Cipher text multiplication then add Complete')


print('')
print('Starting Decryption Tests')
print('=========================')
if group.Dec_src(sk,pk,C0_src,ECtable)==M0:
    print('Decrypt_src_C0:PASS')
else:
    print('Decrypt_src_C0:FAIL')

if group.Dec_src(sk,pk,C1_src,ECtable)==M1:
    print('Decrypt_src_C1:PASS')
else:
    print('Decrypt_src_C1:FAIL')

if group.Dec_tgt(sk,pk,C0_tgt,Ftable)==M0:
    print('Decrypt_tgt_C0:PASS')
else:
    print('Decrypt_tgt_C0:FAIL')

if group.Dec_tgt(sk,pk,C1_tgt,Ftable)==M1:
    print('Decrypt_tgt_C1:PASS')
else:
    print('Decrypt_tgt_C1:FAIL')

if group.Dec_tgt(sk,pk,M0plusM1Tgt,Ftable)==M0+M1:
    print('Decrypt_tgt_C0+C1:PASS')
else:
    print('Decrypt_tgt_C0+C1:FAIL')

if group.Dec_src(sk,pk,M0plusM1,ECtable)==M0+M1:
    print('Decrypt_src_C0+C1:PASS')
else:
    print('Decrypt_src_C0+C1:FAIL')

if group.Dec_tgt(sk,pk,M0xM1,Ftable)==(M0*M1):
    print('Decrypt_tgt_C0*C1:PASS')
else:
    print('Decrypt_tgt_C0*C1:FAIL')

if group.Dec_tgt(sk,pk,M0xM1plusM1,Ftable)==(M0*M1)+M1:
    print('Decrypt_tgt_C0*C1+C1:PASS')
else:
    print('Decrypt_tgt_C0*C1+C1:FAIL')

switched = group.sim_switch(sk,pk,M0xM1)

if group.Dec_src(sk,pk,switched,ECtable)==M0*M1:
    print('Decrypt_src_C0*C1Switch:PASS')
else:
    print('Decrypt_src_C0*C1Switch:FAIL')


