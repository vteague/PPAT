import time
from crypto.cryptofield import CryptoField
from crypto.memtable import MemTable
from crypto.efptable import EFpTable
from crypto.efp12table import EFp12Table
from crypto.targetgroup import TargetGroup
import mathTools.otosEC as oEC
import unittest





class TestCrypto(unittest.TestCase):
    def setUp(self):
        self.startTime = time.time()

    def test_EFp_DLog(self):
        field = CryptoField()
        memtable = MemTable()
        table = EFpTable(memtable,field.G)
        Ptuple = oEC.toTupleEFp(field.P)
        table.build(Ptuple, 2 ** 32, 2 ** 12)
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
        

    def test_target_encrypt(self):
        field = CryptoField()
        memtable = MemTable()
        table = EFp12Table(memtable, field)
        group = TargetGroup(field, table)
        key = {}
        key['pk']=group.load_public_key('./data/pubkey.json')
        key['sk']=group.load_secret_key('./data/secretkey.json')
        pk, sk = group.key_gen(key = key)
        table.make_full_Ftable(pk['e'])
        testmsg = 4
        ciphertest = group.encrypt(pk,testmsg)

        dec = group.decrypt(sk, pk, ciphertest)
        self.assertEqual(dec, testmsg, 'Decrytion target fails')
        print(dec)
       
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCrypto)
    unittest.TextTestRunner(verbosity=1).run(suite)