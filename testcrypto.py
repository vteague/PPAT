import time
from crypto.cryptofield import CryptoField
from crypto.memtable import MemTable
from crypto.efptable import EFpTable
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


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCrypto)
    unittest.TextTestRunner(verbosity=1).run(suite)