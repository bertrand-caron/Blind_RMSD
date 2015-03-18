import unittest
import sys
import rmsd
import logging
logging.basicConfig(level=logging.DEBUG, format='    [%(levelname)s] - %(message)s')

class Test_RMSD(unittest.TestCase):
    pass

def test_alignment_generator(points1, points2):
    def test(self):
        sys.stderr.write("\n")
        logging.info("RMSD before alignment: {0:.4f}".format(rmsd.rmsd(points1, points2)))
        points1_aligned = rmsd.alignPointsOnPoints(points1, points2)
        logging.info("RMSD after alignment: {0:.4f}".format(rmsd.rmsd(points1_aligned, points2)))
        self.assertAlmostEqual( rmsd.rmsd(points1_aligned, points2), 0.)
    return test


if __name__ == "__main__":
    batch_tests = (("simple", 
                    [ [0.,0.,0.], [1.,1.,1.] ], 
                    [ [1.,1.,1.], [0.,0.,0.] ]),
                   ("harder", 
                    [ [0.,0.,0.], [1.,1.,1.] ], 
                    [ [1.,1.,1.], [2.,2.,2.] ]),
                   ("even_harder", 
                    [ [0.,0.,0.], [1.,1.,1.], [-1.,-1.,-1.] ], 
                    [ [2.,2.,2.], [3.,3.,3.], [4.,4.,4.] ]),
                   )
    for name, points1, points2 in batch_tests:
        # Create the tests dynamically
        # Source : http://stackoverflow.com/questions/32899/how-to-generate-dynamic-parametrized-unit-tests-in-python
        test = test_alignment_generator(points1, points2)
        setattr(Test_RMSD, "test_" + name, test)
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_RMSD)
    unittest.TextTestRunner(verbosity=4).run(suite)
