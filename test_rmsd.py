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

SHOW_GRAPH = True

# Differentiate -1's
def split_equivalence_group(eq_list):
    accu = 0
    split_eq_list = []
    for eq in eq_list:
        if eq != -1: split_eq_list.append(eq)
        else:
            split_eq_list.append(eq-accu)
            accu += 1
    return split_eq_list

def download_molecule_files(molecule_name, molids):
    for version in [1,2]:
        molid = molids[version - 1]
        for extension in ['pdb', 'yml']:
            file_name = FILE_TEMPLATE.format(molecule_name=molecule_name, extension=extension, version=version)
            if not exists( dirname(file_name)): os.mkdir( dirname(file_name) )
            if not exists(file_name):
                with open(file_name, 'w') as fh:
                    download_url = DOWNLOAD_TEMPLATES[extension].format(molid=molid)
                    print "Downloading: {0}".format(download_url)
                    response = urllib2.urlopen(download_url)
                    fh.write( response.read() )

def molecule_test_alignment_generator(test_datum, expected_rmsd):
    def test(self):
        molecule_name = test_datum['molecule_name']

        download_molecule_files(molecule_name, [test_datum['id1'], test_datum['id2']])
        point_list1 = []
        point_list2 = []
        m1 = pmx.Model(FILE_TEMPLATE.format(molecule_name=molecule_name, version=1, extension='pdb'))
        point_list1 = [ atom.x[:] for atom in m1.atoms]
        m2 = pmx.Model(FILE_TEMPLATE.format(molecule_name=molecule_name, version=2, extension='pdb'))
        point_list2 = [ atom.x[:] for atom in m2.atoms]

        with open(FILE_TEMPLATE.format(molecule_name=molecule_name, version=1, extension='yml')) as fh: data1 = yaml.load(fh.read())
        with open(FILE_TEMPLATE.format(molecule_name=molecule_name, version=2, extension='yml')) as fh: data2 = yaml.load(fh.read())
        flavour_list1 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data1['atoms'].items()])
        flavour_list2 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data2['atoms'].items()])
        element_list1 = [ atom['type'] for index, atom in data1['atoms'].items()]
        element_list2 = [ atom['type'] for index, atom in data2['atoms'].items()]

        sys.stderr.write("\n")
        logging.info("RMSD before alignment: {0:.4f}".format(rmsd.rmsd(point_list1, point_list2)))

        aligned_point_list1 = rmsd.alignPointsOnPoints(point_list1, point_list2, silent=False, use_AD=False, element_list1=element_list1, element_list2=element_list2, flavour_list1=flavour_list1, flavour_list2=flavour_list2, show_graph=SHOW_GRAPH, rmsd_tolerance=expected_rmsd)
        for i, atom in enumerate(m1.atoms):
            atom.x = aligned_point_list1[i]
        m1.write('testing/{molecule_name}1_aligned.pdb'.format(molecule_name=molecule_name))

        logging.info("RMSD after alignment: {0:.4f}".format(rmsd.rmsd(aligned_point_list1, point_list2)))
        logging.info("Maximum Tolerated RMSD: {0:.4f}".format(expected_rmsd))
        self.assertLessEqual( rmsd.rmsd(aligned_point_list1, point_list2), expected_rmsd)
    return test

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
