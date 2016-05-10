from Blind_RMSD.pdb import pdb_data_for, align_pdb_on_pdb

def pdb_str(test_ID):
    with open('data/{0}.pdb'.format(test_ID)) as fh:
        return fh.read()

TEST_PAIRS = (
    (1, 2),
    (3, 4),
    (5, 6),
    (7, 8),
)

if __name__ == '__main__':
    for test_pair in TEST_PAIRS:
        pdb_data = [pdb_str(test_ID) for test_ID in test_pair]
        aligned_pdb, alignment_score = align_pdb_on_pdb(
            reference_pdb_str=pdb_data[0],
            other_pdb_str=pdb_data[1],
            silent=False,
            soft_fail=False,
        )

        with open('data/2_on_1.pdb', 'w') as fh:
            fh.write(aligned_pdb)

        print alignment_score
