import argparse
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
from copy import deepcopy
import shutil
import numpy
numpy.set_printoptions(precision=3, linewidth=300)

numerical_tolerance = 1e-5
scoring_function = rmsd

FILE_TEMPLATE = "testing/{molecule_name}/{molecule_name}{version}.{extension}"

HOST = 'http://compbio.biosci.uq.edu.au/atb-new'
#HOST = 'http://scmb-atbweb.biosci.uq.edu.au/atb'
API_TOKEN = 'E1A54AB5008F1E772EBC3A51BAEE98BF'

SEARCH_API_URL = '{HOST}/api/current/molecules/search.py?InChI={inchi}&api_token={API_TOKEN}'

DOWNLOAD_TEMPLATES = { 'pdb': '{HOST}/download.py?molid={molid}&outputType=top&dbfile=pdb_allatom_optimised',
                       'yml': '{HOST}/api/current/molecules/generate_mol_data.py?molid={molid}&api_token={API_TOKEN}',
                     }

SHOW_GRAPH = False

SCHEDULED_FOR_DELETION_MOLECULES_FILE = 'testing/{molecule_name}/delete_indexes.ids'
DELETION_THRESHOLD = 1E-2

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

def IDs_for_InChI(inchi):
    url = SEARCH_API_URL.format(HOST=HOST, inchi=urllib2.quote(inchi), API_TOKEN=API_TOKEN)
    print "Getting molIDs for InChI {inchi} at: {url}".format(inchi=inchi, url=url)
    response = urllib2.urlopen(url)
    data = yaml.load(response.read())
    sorted_molecules = sorted(data['molecules'], key=lambda molecule: (not molecule['has_TI'], molecule['molid']))
    molids = map(lambda molecule: molecule['molid'], sorted_molecules)
    print 'Results: {0}'.format(molids)
    return molids

def download_molecule_files(molecule_name, inchi):
    molids =  IDs_for_InChI(inchi)
    for version, molid  in enumerate(molids):
        try:
            for extension in ['pdb', 'yml']:
                file_name = FILE_TEMPLATE.format(molecule_name=molecule_name, extension=extension, version=version)
                if not exists( dirname(file_name)): os.mkdir( dirname(file_name) )
                if not exists(file_name):
                    with open(file_name, 'w') as fh:
                        download_url = DOWNLOAD_TEMPLATES[extension].format(molid=molid, HOST=HOST, API_TOKEN=API_TOKEN)
                        print "Downloading: {0}".format(download_url)
                        response = urllib2.urlopen(download_url)
                        fh.write( response.read() )
        except Exception, e:
            directory = dirname(FILE_TEMPLATE.format(molecule_name=molecule_name, version='', extension=''))
            if exists(directory): shutil.rmtree(directory)
            raise e
    return molids

def molecule_test_alignment_generator(test_datum):
    def test(self):
        molecule_name = test_datum['molecule_name']
        expected_rmsd = test_datum['expected_rmsd']

        molids = download_molecule_files(molecule_name, test_datum['InChI'])
        molid_to_index = dict(zip(molids, range(len(molids))))
        version1, version2 = map(lambda an_id:molid_to_index[test_datum[an_id]], ['id1', 'id2'])
        file1, file2 = map(lambda version: FILE_TEMPLATE.format(molecule_name=molecule_name, version=version, extension='{extension}'), [version1, version2])

        if len(molids) == 1:
            print "Error: Can't run as there are only one molecule in the test suite."
            return
        m1 = pmx.Model(file1.format(extension='pdb'))
        point_list1 = [ atom.x[:] for atom in m1.atoms]
        m2 = pmx.Model(file2.format(extension='pdb'))
        point_list2 = [ atom.x[:] for atom in m2.atoms]

        point_lists = [point_list1, point_list2]

        with open(file1.format(extension='yml')) as fh: data1 = yaml.load(fh.read())
        with open(file2.format(extension='yml')) as fh: data2 = yaml.load(fh.read())
        flavour_list1 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data1['atoms'].items()])
        flavour_list2 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data2['atoms'].items()])
        element_list1 = [ atom['type'] for index, atom in data1['atoms'].items()]
        element_list2 = [ atom['type'] for index, atom in data2['atoms'].items()]

        flavour_lists, element_lists = [flavour_list1, flavour_list2], [element_list1, element_list2]
        
        def bond_matrix(data):
            def bond_line(data):
                return [ 0 if index not in atom['conn'] else 1 for index in map(lambda x:x+1, range(0, len(data['atoms'])))]
            return [ bond_line(data) for _, atom in sorted(data['atoms'].items())]
               
        bonds = map(bond_matrix, [data1, data2])

        sys.stderr.write("\n")
        #logging.info("Score before alignment: {0:.4f}".format(scoring_function(point_list1, point_list2)))

        aligned_point_list1, best_score = align.pointsOnPoints(deepcopy(point_lists), silent=False, use_AD=False, element_lists=element_lists, flavour_lists=flavour_lists, show_graph=SHOW_GRAPH, score_tolerance=expected_rmsd, bonds=bonds )
        #for i, atom in enumerate(m1.atoms):
        #    atom.x = aligned_point_list1[i]
        #m1.write(FILE_TEMPLATE.format(molecule_name=molecule_name, version="1_aligned_on_0", extension='pdb'))

        logging.info("Score after alignment: {0:.4f}".format(best_score))
        logging.info("Maximum Tolerated Score: {0:.4f}".format(expected_rmsd))
        logging.info("To debug these results, run 'pymol {0} {1}'".format( *map(lambda file: file.format(extension='pdb'), [FILE_TEMPLATE.format(molecule_name=molecule_name, version='{0}_aligned_on_{1}'.format(version2, version1), extension='{extension}'), file1]) ))
        self.assertLessEqual( best_score, expected_rmsd)
    return test

def get_distance_matrix(test_datum, overwrite_results=True):
    molecule_name = test_datum['molecule_name']
    expected_rmsd = test_datum['expected_rmsd']

    molids = download_molecule_files(molecule_name, test_datum['InChI'])
    mol_number = len(molids)
    index_to_molid = dict(zip(range(mol_number), molids))

    matrix = numpy.zeros((mol_number, mol_number))
    matrix[:] = numpy.NAN
    matrix_log_file = FILE_TEMPLATE.format(molecule_name=molecule_name, version='', extension='log')

    to_delete_molids = []

    molids_file = FILE_TEMPLATE.format(molecule_name=molecule_name, version='', extension='ids')
    for version1 in range(1, mol_number):
        m1 = pmx.Model(FILE_TEMPLATE.format(molecule_name=molecule_name, version=version1, extension='pdb'))
        point_list1 = [ atom.x[:] for atom in m1.atoms]
        for version2 in range(version1):
            aligned_pdb_file = FILE_TEMPLATE.format(molecule_name=molecule_name, version="{0}_aligned_on_{1}".format(version1, version2), extension='pdb')
            if exists(aligned_pdb_file) and not overwrite_results: continue

            m2 = pmx.Model(FILE_TEMPLATE.format(molecule_name=molecule_name, version=version2, extension='pdb'))
            point_list2 = [ atom.x[:] for atom in m2.atoms]

            point_lists = [point_list1, point_list2]

            with open(FILE_TEMPLATE.format(molecule_name=molecule_name, version=0, extension='yml')) as fh: data1 = yaml.load(fh.read())
            with open(FILE_TEMPLATE.format(molecule_name=molecule_name, version=1, extension='yml')) as fh: data2 = yaml.load(fh.read())
            flavour_list1 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data1['atoms'].items()])
            flavour_list2 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in data2['atoms'].items()])
            element_list1 = [ atom['type'] for index, atom in data1['atoms'].items()]
            element_list2 = [ atom['type'] for index, atom in data2['atoms'].items()]

            flavour_lists, element_lists = [flavour_list1, flavour_list2], [element_list1, element_list2]

            aligned_point_list1, best_score = align.pointsOnPoints(deepcopy(point_lists), silent=True, use_AD=False, element_lists=element_lists, flavour_lists=flavour_lists, show_graph=SHOW_GRAPH, score_tolerance=expected_rmsd)
            matrix[version1, version2] = best_score
            if best_score <= DELETION_THRESHOLD and index_to_molid[version1] not in to_delete_molids:
                to_delete_molids.append(index_to_molid[version1])

            for i, atom in enumerate(m1.atoms):
                atom.x = aligned_point_list1[i]
            m1.write(aligned_pdb_file)
        #break
    if not exists(matrix_log_file): numpy.savetxt(matrix_log_file, matrix, fmt='%4.3f')
    with open(molids_file, 'w') as fh: fh.write("\n".join([ "{0}: {1}".format(i, molid) for i, molid in enumerate(molids)]))
    print 'Debug these results by running: "pymol {0} {1}"'.format(FILE_TEMPLATE.format(molecule_name=molecule_name, version='0', extension='pdb'), FILE_TEMPLATE.format(molecule_name=molecule_name, version='*_aligned_on_0', extension='pdb'))
    print matrix

    # Write the list of molids to delete in a file
    deletion_file = SCHEDULED_FOR_DELETION_MOLECULES_FILE.format(molecule_name=molecule_name)
    with open(deletion_file, 'w') as fj:
        for molid in to_delete_molids:
            fj.write(' wget "{HOST}/api/current/molecules/delete_duplicate.py?molid={molid}&confirm=true"\n'.format(HOST=HOST, molid=molid))
    print "Could delete following molids: {0}".format(to_delete_molids)
    print 'To do so, run: "chmod +x {deletion_file} && ./{deletion_file}"'.format(deletion_file=deletion_file)
    print '\n\n'

class Test_RMSD(unittest.TestCase):
    def run(self, result=None):
        if result.failures or result.errors:
            pass
            #super(Test_RMSD, self).run(result)
        else:
            super(Test_RMSD, self).run(result)

    def assertLessEqual(self, a, b, msg=None):
        if not a <= b + numerical_tolerance:
            self.fail('%s not less than or equal to %s within %f' % (repr(a), repr(b), numerical_tolerance))

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('--only', nargs='*', help="Only run test cases matching these molecule names")
    parser.add_argument('--debug', help="Be overly verbose", action='store_true')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_command_line()
    print args.only
    with open('test_data.yml') as fh:
        test_data = yaml.load(fh.read())
        print "Test data is:\n{0}\n".format(yaml.dump(test_data))
        for test_datum in test_data:
            if args.only and test_datum['molecule_name'] not in args.only: continue
            if 'id1' in test_datum and 'id2' in test_datum:
                test = molecule_test_alignment_generator(test_datum)
                setattr(Test_RMSD, "test_" + test_datum['molecule_name'], test)
            get_distance_matrix(test_datum, silent=not args.debug)
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_RMSD)
    unittest.TextTestRunner(verbosity=4).run(suite)
