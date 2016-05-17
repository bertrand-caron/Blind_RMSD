from Blind_RMSD.pdb import pdb_data_for, align_pdb_on_pdb
from os.path import join

DATA_DIR = 'data'

def pdb_data_file(file_name):
    return join(DATA_DIR, str(file_name) + '.pdb')

def pdb_str(test_ID):
    with open(pdb_data_file(test_ID)) as fh:
        return fh.read()

TEST_PAIRS = (
    (1, 2),
    (3, 4),
    (5, 6),
    (7, 8),
    (9, 10),
    (11, 12),
#    (15, 16),
    (17, 17),
    (18, 18),
)

if __name__ == '__main__':
    for test_pair in TEST_PAIRS:
        print 'Running test: {0}'.format(test_pair)

        pdb_data = [pdb_str(test_ID) for test_ID in test_pair]
        aligned_pdb, alignment_score, alignment_results = align_pdb_on_pdb(
            reference_pdb_str=pdb_data[0],
            other_pdb_str=pdb_data[1],
            silent=False,
            soft_fail=False,
            assert_is_isometry=True,
        )

        fit_file_name = '{other}_on_{reference}'.format(reference=test_pair[0], other=test_pair[1])
        with open(pdb_data_file(fit_file_name), 'w') as fh:
            fh.write(aligned_pdb)

        print alignment_score
        print 'Debug this alignment by running "pymol {0}"'.format(
            ' '.join(
                map(
                    pdb_data_file,
                    test_pair + (fit_file_name,),
                )
            ),
        )
        print
