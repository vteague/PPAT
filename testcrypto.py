import time
from crypto.cryptofield import CryptoField
from crypto.memtable import MemTable
from crypto.efptable import EFpTable
from crypto.efp12table import EFp12Table
from crypto.targetgroup import TargetGroup
from crypto.sourcegroup import SourceGroup
from crypto.group import Group
import mathTools.otosEC as oEC
import unittest





class TestCrypto(unittest.TestCase):
    def setUp(self):
        self.startTime = time.time()
        self.sourcegrp = None
        self.targetgrp = None
        self.public_key = None
        self.secret_key = None
    def init_source_group(self):
        if self.sourcegrp == None:
            field = CryptoField()
            memtable = MemTable()
            table = EFpTable(memtable, field.G)
            self.sourcegrp = SourceGroup(field, table)
            key = {}
            key['pk'] = self.sourcegrp.load_public_key('./data/pubkey2.json')
            key['sk'] = self.sourcegrp.load_secret_key('./data/secretkey2.json')
            pk, sk = self.sourcegrp.key_gen(key=key)
            self.public_key = pk;
            self.secret_key = sk;
            table.build(pk['g'], 2 ** 16, 2 ** 6)

    def init_target_group(self):
        if self.targetgrp == None:
            field = CryptoField()
            memtable = MemTable()
            table = EFp12Table(memtable, field)
            self.targetgrp = TargetGroup(field, table)
            key = {}
            key['pk'] = self.targetgrp.load_public_key('./data/pubkey2.json')
            key['sk'] = self.targetgrp.load_secret_key('./data/secretkey2.json')
            pk, sk = self.targetgrp.key_gen(key=key)
            self.public_key = pk;
            self.secret_key = sk;
            table.build(pk['e'], 2 ** 16, 2 ** 6)
    """    
    def test_EFp_DLog(self):
        field = CryptoField()
        memtable = MemTable()
        table = EFpTable(memtable,field.G)
        Ptuple = oEC.toTupleEFp(field.P)
        table.build(Ptuple, 2 ** 10, 2 ** 6)
        t = time.time() - self.startTime
        print "%s: %.4f" % ("ECTable computation", t)
        #s = getsizeof(table)
        #print "%s: %d %s" % ("Size of ECTable", s, "bytes")
        x = 9
        y = oEC.mulECP(field.G, Ptuple, x, sq=False)
        #print y
        t2 = time.time()
        xc = table.extract(Ptuple, y)
        t = time.time() - t2
        print "%s: %.4f" % ("DL Extraction", t)
        self.assertEqual(x, xc, 'DL extraction fails')
    
    def test_EFp12_DLog(self):
        field = CryptoField()

        memtable = MemTable()
        table = EFp12Table(memtable,field)
        table.make_full_Ftable(field.gt)
        t = time.time() - self.startTime
        print "%s: %.4f" % ("ECTable computation", t)
        
    """
    
    def test_target_encdec(self):
        self.init_target_group()
        for j in range(0,1):
            testmsg = 6
            start = time.time()
            ciphertest = self.targetgrp.encrypt(self.public_key,testmsg)
            t = time.time() - start
            print "%s: %.4f" % ("Target Encrypt", t)
            start = time.time()
            dec = self.targetgrp.decrypt(self.secret_key, self.public_key, ciphertest)
            t = time.time() - start
            print "%s: %.4f" % ("Target Decrypt", t)
            
            self.assertEqual(dec, testmsg, 'Decrytion target fails')

    def test_target_addition(self):
        self.init_target_group()
        testmsgone = 2
        testmsgtwo = 6
        cipherone = self.targetgrp.encrypt(self.public_key,testmsgone)
        ciphertwo = self.targetgrp.encrypt(self.public_key,testmsgtwo)
        start = time.time()
        cipheroneplustwo = self.targetgrp.add(self.public_key,cipherone,ciphertwo)
        t = time.time() - start
        print "%s: %.4f" % ("Target Addition", t)
            
        dec = self.targetgrp.decrypt(self.secret_key, self.public_key, cipheroneplustwo)
        self.assertEqual(dec, testmsgone + testmsgtwo, 'addition target fails')
        
    def test_source_encdec(self):
        self.init_source_group()
        for j in range(1,2):
            testmsg = j
            start = time.time()
            ciphertest = self.sourcegrp.encrypt(self.public_key, testmsg)
            t = time.time() - start
            print "%s: %.4f" % ("Source Encrypt", t)
            start = time.time()
            dec = self.sourcegrp.decrypt(self.secret_key, self.public_key, ciphertest)
            t = time.time() - start
            print "%s: %.4f" % ("Source Encrypt", t)
            
            self.assertEqual(dec, testmsg, 'Decrytion source fails')

    def test_source_addition(self):
        self.init_source_group()
        testmsgone = 2
        testmsgtwo = 6
        cipherone = self.sourcegrp.encrypt(self.public_key,testmsgone)
        ciphertwo = self.sourcegrp.encrypt(self.public_key,testmsgtwo)
        start = time.time()
        cipheroneplustwo = self.sourcegrp.add(self.public_key,cipherone,ciphertwo)
        t = time.time() - start
        print "%s: %.4f" % ("Source Addition", t)
            
        dec = self.sourcegrp.decrypt(self.secret_key, self.public_key, cipheroneplustwo)
        self.assertEqual(dec, testmsgone + testmsgtwo, 'addition source fails')
    
    def test_source_multiply(self):
        self.init_source_group()
        self.init_target_group()
        testmsgone = 10
        testmsgtwo = 2
        cipherone = self.sourcegrp.encrypt(self.public_key,testmsgone)
        ciphertwo = self.sourcegrp.encrypt(self.public_key,testmsgtwo)
        start = time.time()
        cipheroneplustwo = self.sourcegrp.multiply(self.public_key,cipherone,ciphertwo)
        t = time.time() - start
        print "%s: %.4f" % ("Source Multiply", t)
        
        dec = self.targetgrp.decrypt(self.secret_key, self.public_key, cipheroneplustwo)
        self.assertEqual(dec, testmsgone * testmsgtwo, 'multiply source fails')

    def test_source_multiply_then_add(self):
        self.init_source_group()
        self.init_target_group()
        testmsgone = 2
        testmsgtwo = 6
        testmsgthree = 3
        cipherone = self.sourcegrp.encrypt(self.public_key,testmsgone)
        ciphertwo = self.sourcegrp.encrypt(self.public_key,testmsgtwo)
        cipherthree = self.targetgrp.encrypt(self.public_key,testmsgthree)    
        cipheronetimestwo = self.sourcegrp.multiply(self.public_key,cipherone,ciphertwo)
        ciphertest = self.targetgrp.add(self.public_key,cipheronetimestwo,cipherthree)
        dec = self.targetgrp.decrypt(self.secret_key, self.public_key, ciphertest)
        self.assertEqual(dec, testmsgone * testmsgtwo + testmsgthree, 'multiply source then add target fails')

    def test_switch(self):
        self.init_source_group()
        self.init_target_group()
        testmsgone = 2
        cipherone = self.targetgrp.encrypt(self.public_key, testmsgone)
        start = time.time()
        blinding_factor = Group.generate_rand_int()
        target_bf = self.targetgrp.encrypt(self.public_key,blinding_factor)
        source_bf = self.sourcegrp.encrypt(self.public_key,blinding_factor)
        blinded_cipher = self.targetgrp.add(self.public_key, cipherone, target_bf)
        switched_cipher = Group.switch(self.secret_key, self.public_key,
                                       blinded_cipher,
                                       source_bf, self.sourcegrp, self.targetgrp)
        dec = self.sourcegrp.decrypt(self.secret_key, self.public_key, switched_cipher)
        t = time.time() - start
        print "%s: %.4f" % ("Switch", t)
        
        test_blinded_cipher = self.targetgrp.decrypt(self.secret_key, self.public_key, blinded_cipher)
        self.assertEqual(test_blinded_cipher, testmsgone + blinding_factor, 'blinding failed')

        self.assertEqual(dec, testmsgone, 'multiply source then add target fails')


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCrypto)
    unittest.TextTestRunner(verbosity=1).run(suite)