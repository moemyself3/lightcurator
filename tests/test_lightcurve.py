import unittest
from lightcurator import lightcurve as lc

class TestLightcuratorMethods(unittest.TestCase):

    def test_matchcat(self):
        cc = lc.matchcat(23)
        self.assertEqual(cc, 1)

if __name__ == '__main__':
    unittest.main()
