import argparse
import unittest
import sys
import os
import logging
import pmx
logging.basicConfig(level=logging.DEBUG, format='    [%(levelname)s] - %(message)s')
from yaml import load, dump
import urllib.request, urllib.error, urllib.parse
from os.path import exists, dirname, join
from copy import deepcopy
import shutil
import numpy
numpy.set_printoptions(precision=3, linewidth=300)

from Blind_RMSD.align import pointsOnPoints
from Blind_RMSD.helpers.scoring import rmsd, ad, INFINITE_RMSD
from API_client.api import API
from Blind_RMSD.helpers.moldata import group_by, split_equivalence_group, point_list, flavour_list, element_list, pdb_str
from Blind_RMSD.pdb import pdb_data_for, align_pdb_on_pdb
from Blind_RMSD.helpers.exceptions import Topology_Error

numerical_tolerance = 1e-5
scoring_function = rmsd

FILE_TEMPLATE = "testing/{molecule_name}/{version}.{extension}"

API_TOKEN = 'E1A54AB5008F1E772EBC3A51BAEE98BF'
ATB_HOST = 'https://atb.uq.edu.au'
api = API(api_token=API_TOKEN, debug=True, host=ATB_HOST, timeout=600)

SHOW_GRAPH = False

SCHEDULED_FOR_DELETION_MOLECULES_FILE = 'testing/{molecule_name}/delete_indexes.ids'
DELETION_THRESHOLD = 2E-1
TINY_RMSD_SHOULD_DELETE = 2E-1

TEST_DATA_FILE = 'test_data_2.yml'

ERROR_LOG_FILE = 'log.err'
ERROR_LOG = open(ERROR_LOG_FILE, 'w')

UNITED = True

faulty_inchis = []

def download_molecule_files(molecule_name, inchi):
    def sorted_mols_for_InChI(inchi):
        matches = api.Molecules.search(InChI=inchi)
        sorted_molecules = sorted(matches, key=lambda m: (not m.has_TI, m.molid))
        print('Results: {0} (Inchi was: "{1}")'.format([m.molid for m in sorted_molecules], inchi))
        return sorted_molecules

    print('Testing molecule: {0}'.format(molecule_name))
    molecules =  sorted_mols_for_InChI(inchi)

    for version, molecule  in enumerate(molecules):
        try:
            for extension in ['pdb_aa', 'yml']:
                file_name = FILE_TEMPLATE.format(molecule_name=molecule_name, extension=extension, version=version)
                if not exists( dirname(file_name)): os.mkdir( dirname(file_name) )
                # This was a disaster waiting to happen, don't assume that the mapping molid -> temporary index is permanent (which is it not, since we are deleting molecules !)
                if not exists(file_name) or True: molecule.download_file(fnme=file_name, atb_format=extension)
        except Exception as e:
            directory = dirname(FILE_TEMPLATE.format(molecule_name=molecule_name, version='', extension=''))
            if exists(directory): shutil.rmtree(directory)
            raise e
    return molecules

def molecule_test_alignment_generator(test_datum, verbosity=0):
    def test(self):
        molecule_name = test_datum['molecule_name']
        expected_rmsd = test_datum['expected_rmsd']

        molids = download_molecule_files(molecule_name, test_datum['InChI'])
        molid_to_index = dict(list(zip(molids, list(range(len(molids))))))
        version1, version2 = [molid_to_index[test_datum[an_id]] for an_id in ['id1', 'id2']]
        file1, file2 = [FILE_TEMPLATE.format(molecule_name=molecule_name, version=version, extension='{extension}') for version in [version1, version2]]

        if len(molids) == 1:
            print("Error: Can't run as there are only one molecule in the test suite.")
            return
        raise Exception("This is broken")
        m1 = pmx.Model(file1.format(extension='pdb'))
        point_list1 = [ atom.x[:] for atom in m1.atoms]
        m2 = pmx.Model(file2.format(extension='pdb'))
        point_list2 = [ atom.x[:] for atom in m2.atoms]

        point_lists = [point_list1, point_list2]

        with open(file1.format(extension='yml')) as fh: data1 = load(fh.read())
        with open(file2.format(extension='yml')) as fh: data2 = load(fh.read())
        flavour_list1 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in list(data1['atoms'].items())])
        flavour_list2 = split_equivalence_group([ atom['equivalenceGroup'] for index, atom in list(data2['atoms'].items())])
        element_list1 = [ atom['type'] for index, atom in list(data1['atoms'].items())]
        element_list2 = [ atom['type'] for index, atom in list(data2['atoms'].items())]

        flavour_lists, element_lists = [flavour_list1, flavour_list2], [element_list1, element_list2]

        def bond_matrix(data):
            def bond_line(data):
                return [ 0 if index not in atom['conn'] else 1 for index in [x+1 for x in range(0, len(data['atoms']))]]
            return [ bond_line(data) for _, atom in sorted(data['atoms'].items())]

        sys.stderr.write("\n")
        #logging.info("Score before alignment: {0:.4f}".format(scoring_function(point_list1, point_list2)))

        aligned_point_list1, best_score, extra_points, final_permutation = pointsOnPoints(
            deepcopy(point_lists),
            use_AD=False,
            element_lists=element_lists,
            flavour_lists=flavour_lists,
            show_graph=SHOW_GRAPH,
            score_tolerance=expected_rmsd,
            verbosity=verbosity,
        )
        #for i, atom in enumerate(m1.atoms):
        #    atom.x = aligned_point_list1[i]
        #m1.write(FILE_TEMPLATE.format(molecule_name=molecule_name, version="1_aligned_on_0", extension='pdb'))

        logging.info("Score after alignment: {0:.4f}".format(best_score))
        logging.info("Maximum Tolerated Score: {0:.4f}".format(expected_rmsd))
        logging.info("To debug these results, run 'pymol {0} {1}'".format( *[file.format(extension='pdb') for file in [FILE_TEMPLATE.format(molecule_name=molecule_name, version='{0}_aligned_on_{1}'.format(version2, version1), extension='{extension}'), file1]] ))
        self.assertLessEqual( best_score, expected_rmsd)
    return test

def get_distance_matrix(test_datum, debug=False, no_delete=False, max_matrix_size=None, verbosity=0):
    OVERWRITE_RESULTS = True
    ONLY_DO_ONE_ROW = False
    NEXT_TEST_STR = '\n\n'

    molecule_name = test_datum['molecule_name']
    expected_rmsd = DELETION_THRESHOLD

    molecules = download_molecule_files(molecule_name, test_datum['InChI'])

    if len(molecules) == 1:
        print("Found only 1 molecule matching this InChI string. Good job !" + NEXT_TEST_STR)
        return
    elif len(molecules) == 0:
        print("All instance of this molecule are private. Aborting this particular test.")
        return

    if max_matrix_size:
        molecules = molecules[0:min(max_matrix_size, len(molecules))]

    mol_number = len(molecules)

    matrix = numpy.zeros((mol_number, mol_number))
    matrix[:] = numpy.NAN
    matrix_log_file = FILE_TEMPLATE.format(molecule_name=molecule_name, version='', extension='log')

    to_delete_molecules, to_delete_NOW_molecules = [], []
    pymol_files = []

    for i, mol1 in enumerate(molecules):
        with open(FILE_TEMPLATE.format(molecule_name=molecule_name, version=i, extension='pdb_aa')) as fh:
            data1 = pdb_data_for(fh.read())

        for j, mol2 in enumerate(molecules):
            if j >= i:
                continue

            aligned_pdb_file = FILE_TEMPLATE.format(molecule_name=molecule_name, version="{0}_aligned_on_{1}".format(i, j), extension='pdb')
            if exists(aligned_pdb_file) and not OVERWRITE_RESULTS:
                continue

            with open(FILE_TEMPLATE.format(molecule_name=molecule_name, version=j, extension='pdb_aa')) as fh:
                data2 = pdb_data_for(fh.read())

            try:
                aligned_pdb_str, alignment_score, alignment_results = align_pdb_on_pdb(
                    reference_pdb_data=data1,
                    other_pdb_data=data2,
                    soft_fail=False,
                    verbosity=verbosity,
                )
                assert alignment_score is not INFINITE_RMSD
            except Topology_Error:
                print('WARNING: Faulty inchi: {0}'.format(mol1.inchi))
                faulty_inchis.append(mol1.inchi)
                continue
            except Exception as e:
                print('Error: Failed on matching {0} to {1}; error was {2}'.format(i, j, e))
                ERROR_LOG.write(
                    'ERROR: InChI={inchi}, molids={molids}, msg="{msg}"\n'.format(
                        inchi=mol1.inchi,
                        msg=e,
                        molids=[mol1.molid, mol2.molid],
                    ),
                )

                if debug:
                    # This will throw errors outside of the try block in debug mode
                    raise
                else:
                    continue

            matrix[i, j] = alignment_score
            if alignment_score <= DELETION_THRESHOLD:
                if not mol1 in to_delete_molecules:
                    to_delete_molecules.append(mol1)
                    print('Will delete {0} (canonical_molid: {1})'.format(mol1.molid, mol2.molid))

                    if alignment_score <= TINY_RMSD_SHOULD_DELETE and mol1 not in to_delete_NOW_molecules:
                        to_delete_NOW_molecules.append(mol1)

                    pymol_files.append(FILE_TEMPLATE.format(molecule_name=molecule_name, version='{0}_aligned_on_{1}'.format(i, j), extension='pdb'))
                    pymol_files.append(FILE_TEMPLATE.format(molecule_name=molecule_name, version='{0}'.format(j), extension='pdb'))

            with open(aligned_pdb_file, 'w') as fh:
                fh.write(aligned_pdb_str)
        if ONLY_DO_ONE_ROW:
            break

    if not exists(matrix_log_file):
        numpy.savetxt(matrix_log_file, matrix, fmt='%4.3f')

    print('Debug these results by running: "pymol {0} {1}"'.format(
        FILE_TEMPLATE.format(molecule_name=molecule_name, version='0', extension='pdb'),
        FILE_TEMPLATE.format(molecule_name=molecule_name, version='*_aligned_on_0', extension='pdb'),
    ))
    print(matrix)
    print('View these results online at URL: {0}'.format(
        join(
            ATB_HOST,
            'index.py?tab=blind_rmsd_fit&molids={0}'.format(str([m.molid for m in molecules]).replace(' ', '')),
        ),
    ))

    # Write the list of molids to delete in a file
    deletion_file = SCHEDULED_FOR_DELETION_MOLECULES_FILE.format(molecule_name=molecule_name)
    to_delete_molids = [m.molid for m in to_delete_molecules]
    with open(deletion_file, 'w') as fj:
        fj.write('''
        pymol -M {pymol_files}
        echo "Do you want to continue?(yes/no)"
        read continue
        if [ "$continue" != "yes" ]; then
            exit
        fi
        '''.format(pymol_files=' '.join(pymol_files)))
        for molid in to_delete_molids:
            fj.write(' wget "{HOST}/api/current/molecules/delete_duplicate.py?molid={molid}&confirm=true"\n'.format(HOST=api.host, molid=molid))
    #print "Could delete following molids: {0} (indexes: {1})".format(to_delete_molids, [i for i, mol in enumerate(molecules) if mol.molid in to_delete_molids])
    #print 'To do so, run: "chmod +x {deletion_file} && ./{deletion_file}"'.format(deletion_file=deletion_file)
    if not (no_delete or debug): 
        try:
            print("Deleting NOW the following molids: {0}".format([mol.delete_duplicate() for mol in to_delete_NOW_molecules]))
        except urllib.error.HTTPError as e:
            print('Something went wrong while trying to delte duplicate. Error was: ')
            print(e)
    else:
        print('Running in no_delete / debug mode. Otherwise, would have deleted molids: {0}'.format([mol.molid for mol in to_delete_NOW_molecules]))

    print(NEXT_TEST_STR)

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
    parser.add_argument('--verbosity', help="Verbosity level. Zero is silent.", type=int)
    parser.add_argument('--auto', help="Get the inchis from the API", action='store_true')
    parser.add_argument('--nodelete', help="Do not delete molecules", action='store_true')
    parser.add_argument('--max-matrix-size', help="Maximum size of the distance matrix.", dest='max_matrix_size', default=None, type=int)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    from re import sub
    args = parse_command_line()

    if args.auto:
        test_molecules = api.Molecules.duplicated_inchis(offset=0, limit=1000, min_n_atoms=0, max_n_atoms=1000)
        for i, mol in enumerate(test_molecules):
            if not mol['molecule_name'] or mol['molecule_name'] == '':
                mol['molecule_name'] = 'unknown_mol_{n}'.format(n=i)
            mol['molecule_name'] = sub('[\[\]\(\)~`$]', '', mol['molecule_name']).replace(' ', '_').replace('/', '_')
            if len(mol['molecule_name']) >= 75:
                mol['molecule_name'] = mol['molecule_name'][0:75] + '_' + str(i)
    else:
        with open(TEST_DATA_FILE) as fh: test_molecules = load(fh.read())

    print("Test data is:\n{0}\n".format(dump(test_molecules)))
    for test_datum in test_molecules:
        if args.only and test_datum['molecule_name'] not in args.only: continue
        if 'id1' in test_datum and 'id2' in test_datum:
            test = molecule_test_alignment_generator(test_datum, verbosity=0)
            setattr(Test_RMSD, "test_" + test_datum['molecule_name'], test)
        get_distance_matrix(
            test_datum,
            verbosity=args.verbosity,
            debug=args.debug,
            no_delete=args.nodelete,
            max_matrix_size=args.max_matrix_size,
        )

    print('Faulty inchis')
    print(set(faulty_inchis))

    suite = unittest.TestLoader().loadTestsFromTestCase(Test_RMSD)
    unittest.TextTestRunner(verbosity=4).run(suite)
