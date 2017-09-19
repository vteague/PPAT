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
    def init_target_group(self):
        field = CryptoField()
        memtable = MemTable()
        table = EFp12Table(memtable, field)
        self.targetgrp = TargetGroup(field, table)
        key = {}
        key['pk'] = self.targetgrp.load_public_key('./data/pubkey.json')
        key['sk'] = self.targetgrp.load_secret_key('./data/secretkey.json')
        pk, sk = self.targetgrp.key_gen(key=key)
        self.public_key = pk;
        self.secret_key = sk;
        table.build(pk['e'], 2 ** 4, 2)
        x = 0
        gt = pk['e']
        Gt = pk['Gt']
        gt_oec = oEC.toTupleFp12(gt)
        y = oEC.squareAndMultiplyFp12(Gt, gt_oec, x,
                                         oEC.tmulFp12, oEC.tsqrtFp12, field.Gamma)
        xc = table.extract(gt_oec, y)
        self.assertEqual(xc, x, 'Table lookup failed')
        x = (2 ** 4) + 3
        gt = pk['e']
        Gt = pk['Gt']
        gt_oec = oEC.toTupleFp12(gt)
        y = oEC.squareAndMultiplyFp12(Gt, gt_oec, x,
                                         oEC.tmulFp12, oEC.tsqrtFp12, field.Gamma)
        xc = table.extract(gt_oec, y)
        self.assertEqual(xc, None, 'Table lookup failed')
        
    
    def test_target_encdec(self):
        self.init_target_group()
        testmsg = 9
        """start = time.time()
        ciphertest = self.targetgrp.encrypt(self.public_key,testmsg)
        t = time.time() - start
        print "%s: %.4f" % ("Target Encrypt", t)
        start = time.time()
        dec = self.targetgrp.decrypt(self.secret_key, self.public_key, ciphertest)
        t = time.time() - start
        print "%s: %.4f" % ("Target Decrypt", t)
        
        self.assertEqual(dec, testmsg, 'Decrytion target fails')"""

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCrypto)
    unittest.TextTestRunner(verbosity=1).run(suite)