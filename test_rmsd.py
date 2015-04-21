import unittest
import sys
import os
import logging
import pmx
logging.basicConfig(level=logging.DEBUG, format='    [%(levelname)s] - %(message)s')
import yaml
#sys.path.append("../")
import urllib2
from os.path import exists, dirname

import align
from scoring import rmsd, ad

numerical_tolerance = 1e-5
scoring_function = rmsd if True else ad

def test_alignment_generator(points1, points2, expected_rmsd):
    def test(self):
        sys.stderr.write("\n")
        logging.info("Score before alignment: {0:.4f}".format(scoring_function(points1, points2)))
        points1_aligned = align.pointsOnPoints(points1, points2)
        logging.info("Score after alignment: {0:.4f}".format(scoring_function(points1_aligned, points2)))
        logging.info("Maximum Tolerated Score: {0:.4f}".format(expected_rmsd))
        self.assertLessEqual( scoring_function(points1_aligned, points2), expected_rmsd)
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
        m1 = pmx.Model(FILE_TEMPLATE.format(molecule_name=molecule_name, version=1, extension='pdb'))
        point_list1 = [ atom.x[:] for atom in m1.atoms]
        m2 = pmx.Model(FILE_TEMPLATE.format(molecule_name=molecule_name, version=2, extension='pdb'))
        point_list2 = [ atom.x[:] for atom in m2.atoms]

        point_lists = [point_list1, point_list2]

        with open(FILE_TEMPLATE.format(molecule_name=molecule_name, version=1, extension='yml')) as fh: data1 = yaml.load(fh.read())
        with open(FILE_TEMPLATE.format(molecule_name=molecule_name, version=2, extension='yml')) as fh: data2 = yaml.load(fh.read())
        flavour_list1 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data1['atoms'].items()])
        flavour_list2 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data2['atoms'].items()])
        element_list1 = [ atom['type'] for index, atom in data1['atoms'].items()]
        element_list2 = [ atom['type'] for index, atom in data2['atoms'].items()]

        flavour_lists, element_lists = [flavour_list1, flavour_list2], [element_list1, element_list2]

        sys.stderr.write("\n")
        logging.info("Score before alignment: {0:.4f}".format(scoring_function(point_list1, point_list2)))

        aligned_point_list1 = align.pointsOnPoints(point_lists, silent=False, use_AD=False, element_lists=element_lists, flavour_lists=flavour_lists, show_graph=SHOW_GRAPH, score_tolerance=expected_rmsd)
        for i, atom in enumerate(m1.atoms):
            atom.x = aligned_point_list1[i]
        m1.write(FILE_TEMPLATE.format(molecule_name=molecule_name, version="1_aligned", extension='pdb'))

        logging.info("Score after alignment: {0:.4f}".format(scoring_function(aligned_point_list1, point_list2)))
        logging.info("Maximum Tolerated Score: {0:.4f}".format(expected_rmsd))
        self.assertLessEqual( scoring_function(aligned_point_list1, point_list2), expected_rmsd)
    return test

class Test_RMSD(unittest.TestCase):
    def run(self, result=None):
        if result.failures or result.errors:
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
