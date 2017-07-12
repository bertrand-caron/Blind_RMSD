from argparse import ArgumentParser, Namespace

from Blind_RMSD.pdb import pdb_data_for, align_pdb_on_pdb
from Blind_RMSD.align import DEFAULT_FLAVOURED_KABSCH_MIN_N_UNIQUE_POINTS

def parse_args() -> Namespace:
    parser = ArgumentParser()
    parser.add_argument('--reference')
    parser.add_argument('--other')
    parser.add_argument('--united-atom-fit', action='store_true')
    parser.add_argument('--N', type=int, default=DEFAULT_FLAVOURED_KABSCH_MIN_N_UNIQUE_POINTS)
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    other_pdb_data, reference_pdb_data = map(
        lambda pdb_filepath: pdb_data_for(
            open(pdb_filepath).read(),
            exception_searching_keywords=[],
            united_atom_fit=args.united_atom_fit,
        ),
        (args.other, args.reference)
    )

    print(
        align_pdb_on_pdb(
            other_pdb_data=other_pdb_data,
            reference_pdb_data=reference_pdb_data,
            soft_fail=False,
            debug=True,
            verbosity=100,
            flavoured_kabsch_min_n_unique_points=args.N,
        )
    )
