import unittest
import rmsd

class Test_RMSD(unittest.TestCase):
    def test_simple(self):
        set1 = [ [0.,0.,0.], [1.,1.,1.] ]
        set2 = [ [1.,1.,1.], [0.,0.,0.] ]
        rmsd.rmsd(set1, set2)
        set1_aligned = rmsd.alignPointsOnPoints(set1, set2)
        self.assertAlmostEqual( rmsd.rmsd(set2, set1_aligned), 0.)

    def test_harder(self):
        set1 = [ [0.,0.,0.], [1.,1.,1.] ]
        set2 = [ [1.,1.,1.], [2.,2.,2.] ]
        rmsd.rmsd(set1, set2)
        set1_aligned = rmsd.alignPointsOnPoints(set1, set2)
        self.assertAlmostEqual( rmsd.rmsd(set2, set1_aligned), 0.)

    def test_even_harder(self):
        set1 = [ [0.,0.,0.], [1.,1.,1.], [-1.,-1.,-1.] ]
        set2 = [ [2.,2.,2.], [3.,3.,3.], [4.,4.,4.] ]
        rmsd.rmsd(set1, set2)
        set1_aligned = rmsd.alignPointsOnPoints(set1, set2)
        self.assertAlmostEqual( rmsd.rmsd(set2, set1_aligned), 0.)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_RMSD)
    unittest.TextTestRunner(verbosity=4).run(suite)
