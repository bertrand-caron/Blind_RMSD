import unittest
import sys
import os
import rmsd
import logging
import math
from Vector import rotaxis2m, Vector, rotmat
from numpy import pi, array, dot, sqrt
from random import uniform
import pmx
from copy import deepcopy
logging.basicConfig(level=logging.DEBUG, format='    [%(levelname)s] - %(message)s')
import yaml
sys.path.append("../")
import urllib2
from os.path import exists, dirname, isdir

numerical_tolerance = 1e-5

def test_alignment_generator(points1, points2, expected_rmsd):
    def test(self):
        sys.stderr.write("\n")
        logging.info("RMSD before alignment: {0:.4f}".format(rmsd.rmsd(points1, points2)))
        points1_aligned = rmsd.alignPointsOnPoints(points1, points2)
        logging.info("RMSD after alignment: {0:.4f}".format(rmsd.rmsd(points1_aligned, points2)))
        logging.info("Maximum Tolerated RMSD: {0:.4f}".format(expected_rmsd))
        self.assertLessEqual( rmsd.rmsd(points1_aligned, points2), expected_rmsd)
    return test

FILE_TEMPLATE = "testing/{molecule_name}/{molecule_name}{version}.{extension}"

DOWNLOAD_TEMPLATES = { 'pdb': 'http://compbio.biosci.uq.edu.au/atb/download.py?molid={molid}&outputType=top&dbfile=pdb_fromuser',
                       'yml': 'http://compbio.biosci.uq.edu.au/atb-new/api/molecules/mol_data.py?molid={molid}',
                     }

def download_molecule_files(molecule_name, molids):
    for version in [1,2]:
        molid = molids[version - 1]
        for extension in ['pdb', 'yml']:
            file_name = FILE_TEMPLATE.format(molecule_name=molecule_name, extension=extension, version=version)
            if not exists( dirname(file_name)): os.mkdir( dirname(file_name) )
            if not exists(file_name):
                with open(file_name, 'w') as fh:
                    response = urllib2.urlopen(DOWNLOAD_TEMPLATES[extension])
                    fh.write( response.read() )

def molecule_test_alignment_generator(test_datum, expected_rmsd):
    def test(self):
        molecule_name = test_datum['molecule_name']

        download_molecule_files(molecule_name, [test_datum['id1'], test_datum['id2']])
        point_list1 = []
        point_list2 = []
        m1 = pmx.Model('testing/{molecule_name}1.pdb'.format(molecule_name=molecule_name))
        point_list1 = [ atom.x[:] for atom in m1.atoms]
        m2 = pmx.Model('testing/{molecule_name}2.pdb'.format(molecule_name=molecule_name))
        point_list2 = [ atom.x[:] for atom in m2.atoms]

        sys.stderr.write("\n")
        logging.info("RMSD before alignment: {0:.4f}".format(rmsd.rmsd(point_list1, point_list2)))

        aligned_point_list1 = rmsd.alignPointsOnPoints(point_list1, point_list2, silent=False, use_AD=False, flavour_list1=[atom.symbol for atom in m1.atoms], flavour_list2=[atom.symbol for atom in m2.atoms], show_graph=False, rmsd_tolerance=expected_rmsd)
        for i, atom in enumerate(m1.atoms):
            atom.x = aligned_point_list1[i]
        m1.write('testing/{molecule_name}1_aligned.pdb'.format(molecule_name=molecule_name))

        logging.info("RMSD after alignment: {0:.4f}".format(rmsd.rmsd(aligned_point_list1, point_list2)))
        logging.info("Maximum Tolerated RMSD: {0:.4f}".format(expected_rmsd))
        self.assertLessEqual( rmsd.rmsd(aligned_point_list1, point_list2), expected_rmsd)
    return test

TEST_MOLECULES = ['methanol', 'ethanol', 'benzene', 'sulfuric_acid']

class Test_RMSD(unittest.TestCase):
    def run(self, result=None):
        if result.failures or result.errors:
          #print "Aborting due to first failed test ..."
            pass
        else:
            super(Test_RMSD, self).run(result)

    def assertLessEqual(self, a, b, msg=None):
        if not a <= b + numerical_tolerance:
            self.fail('%s not less than or equal to %s within %f' % (repr(a), repr(b), numerical_tolerance))

if __name__ == "__main__":
    with open('test_data.yml') as fh:
        test_data = yaml.load(fh.read())
        print "Test data is:\n{0}\n".format(yaml.dump(test_data))
        for test_datum in test_data:
            test = molecule_test_alignment_generator(test_datum, 0.15)
            setattr(Test_RMSD, "test_" + test_datum['molecule_name'], test)
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_RMSD)
    unittest.TextTestRunner(verbosity=4).run(suite)
