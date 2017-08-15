#Copyright 2017 Ilya Marchenko
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

# -*- coding: utf-8 -*-
import sys
sys.path.append('C:\Research\STAR Vote\PPAT-master')

import time
import numpy as np
import itertools
import mathTools.field as field
import mathTools.ellipticCurve as ellipticCurve
import mathTools.pairing as pairing
import ppat.ppats
import mathTools.otosEC as oEC
import gmpy2 as gmpy
from Crypto.Random.random import randint
from Crypto.Random.random import getrandbits
from mathTools.otosEC import OptimAtePairing as e_hat    

# maximum size in bits of secret
MAX_BITS = 100

range = lambda stop: iter(itertools.count().next, stop) # Range for iterating over large numbers (larger than int)

#c = gmpy.mpz(2) # p is 160-bit long
c = gmpy.mpz(1) # p is 256-bit long
d = gmpy.mpz(1)
b = c**4+d**6 # b = c**4+d**6
#u = gmpy.mpz(-(2**38 + 2**28 + 1 )) # p is 160-bit long
u = gmpy.mpz(-(2**62 + 2**55 + 1 )) # p is 256-bit long
#p = 36*u**4 + 36*u**3 + 24*u**2 + 6*u + 1
def pr(u):
    return 36*u**4 + 36*u**3 + 24*u**2 + 6*u + 1
def nr(u):
    return 36*u**4 + 36*u**3 + 18*u**2 + 6*u + 1
p = pr(u)
n = nr(u)

#p = 1300829
#n = 1299721


#n = 36*u**4 + 36*u**3 + 18*u**2 + 6*u + 1
#n is 160-bit long with low HW
t = 6*u**2 + 1

##### Fp #####
Fp = field.Field(p)
fp0 = Fp.zero()
fp1 = Fp.one()

print Fp, " ...done"
##### E[Fp] #####
C = ellipticCurve.Curve(fp0,b*fp1,Fp) # Y**2 = X**3+b
PInf = ellipticCurve.ECPoint(infty = True)
EFp = ellipticCurve.ECGroup(Fp,C,PInf)
P = EFp.elem((-d**2)*fp1,(c**2)*fp1)  # P  is a generetor of EFp of order n (n*P = Pinf)


##### Fp2b #####
poly1 = field.polynom(Fp,[fp1,fp0,fp1]) # X**2+1
print poly1

Fp2 = field.ExtensionField(Fp,poly1,rep='i') # A**2 = -1
print Fp2, " ...done"
fp2_0 = Fp2.zero()
fp2_1 = Fp2.one()
fp2_ip = field.polynom(Fp,[fp1,fp0]) # 1*A+0
fp2_i = field.ExtensionFieldElem(Fp2,fp2_ip)
xi = (c**2)*fp2_1+(d**3)*fp2_i # c**2+(d**3)*A (4+i)
cxi = (c**2)*fp2_1-(d**3)*fp2_i # c**2-(d**3)*A
#ixi = 8*fp2bi-8*fp2b1 # 8*A-8
#xi = ixi.invert()
#C2b = EllipticCurve.Curve(fp2b0, 3*ixi,Fp2b) # Y**2 = X**3+3*(8*A-8)
C2 = ellipticCurve.Curve(fp2_0, cxi,Fp2) # Y**2 = X**3+c**2-(d**3)*A The twisted curve
PInf2 = ellipticCurve.ECPoint(infty = True)
EFp2 = ellipticCurve.ECGroup(Fp2,C2,PInf2)

u0 = EFp2.elem((-d)*fp2_i,c*fp2_1) #EC point (-d*A,c)
h = 2*p-n
Q = u0*h # Q is a generator of G2 of order n

r= randint(1,int(n))
s= randint(1,int(n))
rP = r*P
sQ = s*Q


##### Fp6 #####
poly3 = field.polynom(Fp2,[fp2_1,fp2_0,fp2_0,-xi]) #X**3-xi
Fp6 = field.ExtensionField(Fp2,poly3)
fp6_0 = Fp6.zero()
fp6_1 = Fp6.one()
fp6_xi = Fp6.elem(xi) # xi in Fp6

##### Fp12 #####
poly6 = field.polynom(Fp6,[fp6_1,fp6_0,-fp6_xi]) # X**2-xi
Fp12 = field.ExtensionField(Fp6,poly6)
print Fp12, " ...done"
fp12_0 = Fp12.zero()
fp12_1 = Fp12.one()
C12 = ellipticCurve.Curve(fp12_0,b*fp12_1,Fp12) # Y**2 = X**3+b
PInf12 = ellipticCurve.ECPoint(infty = True)
EFp12 = ellipticCurve.ECGroup(Fp12,C12,PInf12)

Qpr = oEC.psi(EFp12,Q) # Qpr lives in E[Fp12b]
Pair = pairing.Pairing(EFp,EFp12,C,P,Q,n,Qpr,oEC.frobenius,oEC.prec_gamma(Fp12,u,c,d))

############### Crypto Functions #############
Gamma = oEC.prec_gamma(Fp12,u,c,d)

def e(g, h):
    """Evaluates the bilinear operator on pairs of elements of G and H.
    Uses Karatsuba's trick.
    Arguments:
        pp, the public parameters describing the groups
        g, a pair of elements of G
        h, a pair of elements of H
    Returns a triple of elements of Gt
    """
    r0 = e_hat(g[0], h[0], Pair)
    r2 = e_hat(g[1], h[1], Pair)
    r1 = e_hat(g[0] + g[1], h[0] + h[1], Pair) * pp['Gt'].invert(r0 * r2)
    return (r0, r1, r2)

def Setup():
    
    pp={'p': p, 'G': EFp, 'H': EFp2,'Gt': Fp12,'e': e_hat}
    return pp
    
def KeyGen(pp):
    
    G = pp['G'] # Representation of G
    H = pp['H'] # Representation of H    
    s = getrandbits(MAX_BITS)
    
    P1 = P*randint(0,int(n))     
    G1 = (G.neg(P1)*s,P1) # Description of G1 - (g^{-s},g)
    
    Q1 = Q*randint(0,int(n))   
    H1 = (H.neg(Q1)*s,Q1) # Description of H1 - (h^{-s},h)
    
    Gt=pp['Gt']   
    
    #g = P*randint(0,int(n)) # Random element of G
    #h = Q*randint(0,int(n))# Random element of H
    g = P1
    h = Q1

    gt = e_hat(g,h,Pair)
    
    pk = {'G':G,'G1':G1,'H':H,'H1':H1,'Gt':Gt,'g':g,'h':h, 'e':gt}

    def pi_1(g):
        if g[1]==1:
            output = g[0]
        else:
            output = g[0]+(g[1]*s) # pi_1((g1,g2)) = g1 * g2^s
        return output
    
    def pi_2(h):
        if h[1]==1:
            output = h[0]
        else:
            output = h[0]+(h[1]*s) # pi_2((h1,h2)) = h1 * h2^s
        return output 

    def pi_t(gt):
        if gt['C1'] == 1:
            output = gt[0]
        else:
            output = gt['C0']*(gt['C1']**s)*(gt['C2']**(s**2)) # pi_t((gt1,gt2,gt3)) = gt1 * gt2^s * gt3^{s^2}
        return output 

    sk={'pi_1':pi_1, 'pi_2':pi_2,'pi_t':pi_t}
    return pk, sk
    
def Enc_src(pk, M):
    
    a = randint(1,int(p))
    b = randint(1,int(p))
        
    G1 = pk['G1']
    H1 = pk['H1']
    
    g1 = (G1[0]*a,G1[1]*a)
    h1 = (H1[0]*b,H1[1]*b)

    C_0 = (pk['g']*M + g1[0],g1[1]) # C_0 = (g^M * g^{-as}, g^a)
    C_1 = (pk['h']*M + h1[0],h1[1]) # C_1 = (h^M * h^{-bs}, h^b)
    
    C = {'C0':C_0,'C1':C_1}    
    return C

def Enc_tgt(pk,M):
        
    a = randint(1,int(p))
    b = randint(1,int(p))
    
    G1 = pk['G1']
    H1 = pk['H1']
    Gt = pk['Gt']
    
    gt = pk['e']
    
    g1 = (G1[0]*a,G1[1]*a)
    h1 = (H1[0]*b,H1[1]*b)
    
    
    C0 = (gt**M)*e_hat(G1[1],h1[0],Pair)*e_hat(g1[0],H1[1],Pair)
    C1 = e_hat(G1[1],h1[1],Pair)*e_hat(g1[1],H1[1],Pair)
    C2 = Gt.one()    
        
    C = {'C0':C0,'C1':C1,'C2':C2}
    return C
    
def Multiply_src(pk,C0,C1): 

    Gt = pk['Gt']
    
    G1 = pk['G1']
    H1 = pk['H1']
    
    g = pk['g']
    h = pk['h']  
    a = randint(1,int(p))
    b = randint(1,int(p))
    
    #e(C0,C1):
    #eC = (e_hat(C0[0],C1[0],Pair),e_hat(C0[0],C1[1],Pair)*e_hat(C0[1],C1[0],Pair),e_hat(C0[1],C1[1],Pair))
    eC = e(C0, C1)
    #e(g,h1):
    g1 = (e_hat(g,H1[0]*a,Pair),e_hat(g,H1[1]*a,Pair),Gt.one())
    #e(g1,h):
    h1 = (e_hat(G1[0]*b,h,Pair),e_hat(G1[1]*b,h,Pair),Gt.one())
    
    C = {'C0':(eC[0]*g1[0])*h1[0],'C1':eC[1]*(g1[1]*h1[1]),'C2':eC[2]}
    
    return C # C = e(C0,C1) * e(g,h1) * e(g1,h)

def Add_src(pk,C,C_prime):
    
    G1 = pk['G1']
    H1 = pk['H1']
    
    a = randint(1, int(n))
    b = randint(1, int(n))
    g1 = (G1[0]*a,G1[1]*a)
    h1 = (H1[0]*b,H1[1]*b)

    C_0 = ((C_prime['C0'][0]+g1[0])+C['C0'][0],(C_prime['C0'][1]+g1[1])+C['C0'][1])
    C_1 = ((C_prime['C1'][0]+h1[0])+C['C1'][0],(C_prime['C1'][1]+h1[1])+C['C1'][1])
    
    C_doubleprime = {'C0':C_0,'C1':C_1} # C'' = ( C0*C'0*g1 , C1*C'1*h1 )
        
    return C_doubleprime
    
def Add_tgt(pk,C,C_prime):
    
    Gt = pk['Gt']
    
    G1 = pk['G1']
    H1 = pk['H1']
    
    g = pk['g']
    h = pk['h']    
    a = randint(1,int(p))
    b = randint(1,int(p))
    
    #e(g,h1):
    g1 = (e_hat(g,(H1[0]*a),Pair),e_hat(g,(H1[1]*a),Pair),Gt.one())
    #e(g1,h):
    h1 = (e_hat((G1[0]*b),h,Pair),e_hat((G1[1]*b),h,Pair),Gt.one())
    
    C0 = C['C0']*C_prime['C0']*g1[0]*h1[0]
    C1 = C['C1']*C_prime['C1']*g1[1]*h1[1]
    C2 = C['C2']*C_prime['C2']*g1[2]*h1[2]
    
    C_doubleprime = {'C0':C0,'C1':C1,'C2':C2} # C'' = e(C,C') * e(g,h1) * e(g1,h)
    return C_doubleprime
    
def log_group(a,b,Group, table = {}): # Computes x such that a^x = b over a group of order p using baby step-giant step

    #baby_steps = {}
    
    x = Group.neg(a)*(2**16)    
    gamma = b
    
    if table == {}:
        print 'No Table Found'        
        table = make_ECtable(Group,a)

    for i in range(2**16):
        if gamma in table:
            return i*(2**16)+table[gamma]
        else: gamma = gamma+x

    return "No Match"

def log_field(a,b,Field,table = {}): # Computes x such that a^x = b over a field of order p using baby step-giant step

    #baby_steps = {}
    
    x = oEC.toTupleFp12(a.invert()**(2**16))
    gamma = oEC.toTupleFp12(b)

    if table == {}:
        print 'No Table Found'
        table = make_Ftable(Field,a)
        
    for i in range(2**16):
        if gamma in table:
            return (i)*(2**16)+table[gamma]
        else: gamma = oEC.tmulFp12(Fp12,gamma,x,Gamma)

    return "No Match"

def Dec_src(sk,pk,C,table = {}):
        
    M = log_group(sk['pi_1']((pk['g'],1)),sk['pi_1'](C['C0']),pk['G'],table)

    return M
    #if sk['pi_2'](C['C1']).__eq__(pk['H'].mul(sk['pi_2']((pk['h'],1)),M)):
    #    return M
    #else:
    #    print "M != M'"
    #    return 3.14159 # Dummy value returned    
    
def Dec_tgt(sk,pk,C,table = {}):
        
    g = pk['g']
    h = pk['h']        
        
    M = log_field(e_hat(g,h,Pair),sk['pi_t'](C),pk['Gt'],table) 
    return M

def make_ECtable(Group,point):
    #This function makes a multiplication table to aid in computing discrete 
    #logarithims to speed up decryption of multiple messages encrypted with the
    #same public/private key
 
#VT: I edited this to make it terminate at exactly the index it previously did. 
# However, I'd assume that the omission of the j=2**16 + 1 case makes no 
# difference so we should just let the range terminate at 2**16 as before.
# But I'm reluctant to fiddle in case I mess up something important.

    baby_steps = {}    
    pt = point    
    
    for j in xrange(2**16+1):
        baby_steps[pt] = j+1
        pt = pt+point        

       

    return baby_steps
    
def make_Ftable(Field,elt):
    #This function makes a multiplication table to aid in computing discrete 
    #logarithims to speed up decryption of multiple messages encrypted with the
    #same public/private key
        
    baby_steps = {}
    gt = oEC.toTupleFp12(Field.one())

    for j in xrange(2**16):
        gt = oEC.tmulFp12(Fp12,gt,oEC.toTupleFp12(elt),Gamma)
        baby_steps[gt]=j+1

    return baby_steps


#def make_random_ballot(c,k,pk):
#    
#    M = np.ones((k,c),dtype = tuple)    
#    
#    chosen_candidates = []
#    for i in range(k):
#        x = np.random.choice([p for p in range(c) if p not in chosen_candidates])
#        M[i][x] = 2
#        chosen_candidates.append(x)
#
#    for j in range(k):
#        for l in range(c):
#            M[j][l] = Enc_src(pk, M[j][l])
#    return M        


############ Testing Crypto Functions ##########
pp = Setup()
print 'Setup Complete'

pk,sk = KeyGen(pp)
print 'KeyGen Complete'

#M0 = 123
#M1 = 2**30
ECtable = make_ECtable(EFp,pk['g'])
Ftable = make_Ftable(Fp12,pk['e'])
print 'Multiplication Tables Created'

M0 = 0
M1 = 1

C0_src = Enc_src(pk,M0)
C1_src = Enc_src(pk,M1)
print 'Enc_src Complete'

M0xM1= Multiply_src(pk,C0_src['C0'],C1_src['C1'])



# sample candidates
# 0 0 1
# 0 1 0
# 1 0 0

# 0 0 1
# 1 0 0
# 0 1 0

# 0 1 0
# 0 0 1
# 1 0 0

# 1 0 0
# 0 1 0
# 0 0 1

# 0 1 0
# 1 0 0
# 0 0 1