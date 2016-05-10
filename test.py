from Blind_RMSD.pdb import pdb_data_for, align_pdb_on_pdb

def pdb_str(test_ID):
    with open('data/1.pdb') as fh:
        return fh.read()

if __name__ == '__main__':
    pdb_data = [pdb_str(test_ID) for test_ID in (1, 2)]
    aligned_pdb, alignment_score = align_pdb_on_pdb(
        reference_pdb_str=pdb_data[0],
        other_pdb_str=pdb_data[1],
        silent=False,
        soft_fail=True,
    )

    with open('data/2_on_1.pdb', 'w') as fh:
        fh.write(aligned_pdb)
