from traceback import format_exc
from collections import namedtuple
from os.path import abspath, join, dirname, exists
from os import mkdir
from scipy.spatial.distance import squareform
from functools import reduce
from typing import NamedTuple, Any, List, Optional, Tuple, Dict

from Blind_RMSD.helpers.log import log
from Blind_RMSD.helpers.moldata import flavour_list, point_list, aligned_pdb_str, united_hydrogens_point_list
from Blind_RMSD.align import pointsOnPoints, FAILED_ALIGNMENT, NULL_PDB_WRITING_FCT
from Blind_RMSD.helpers.exceptions import Topology_Error, Permutation_Not_Found_Error

from chemical_equivalence.calcChemEquivalency import partial_mol_data_for_pdbstr, ALL_EXCEPTION_SEARCHING_KEYWORDS, MolDataFailure

UNITED_RMSD_FIT = True

PDB, RMSD = str, float

PDB_Data = NamedTuple(
    'PDB_Data',
    [
        ('data', Any),
        ('point_lists', Any),
        ('flavour_lists', Any),
        ('extra_points_lists', Any),
        ('pdb_str', PDB),
        ('united_atom_fit', bool),
    ],
)

FAILED_ALIGNMENT_PDB_STR = ''

Alignment_Results = NamedTuple(
    'Alignment_Results',
    [
        ('success', bool),
        ('error_exc', Any),
        ('permutation', Dict[int, int]),
        ('atom_names_permutation', Dict[int, int]),
        ('total_fit_points', int),
        ('kabsch_fit_point', int),
    ],
)

DEBUG_DIR = join(dirname(abspath(__file__)), 'debug')

def pdb_data_for(
    pdb_str: str,
    exception_searching_keywords: List[str] = ALL_EXCEPTION_SEARCHING_KEYWORDS,
    united_atom_fit: bool = UNITED_RMSD_FIT,
    enforce_single_molecule: bool = True,
) -> PDB_Data:
    data = partial_mol_data_for_pdbstr(
        pdb_str,
        exception_searching_keywords=exception_searching_keywords,
        enforce_single_molecule=enforce_single_molecule,
    ).__dict__

    return PDB_Data(
        data=data,
        point_lists=point_list(data, united_atom_fit),
        flavour_lists=flavour_list(data, united_atom_fit),
        extra_points_lists=united_hydrogens_point_list(data, united_atom_fit),
        pdb_str=pdb_str,
        united_atom_fit=united_atom_fit,
    )

def align_pdb_on_pdb(
    reference_pdb_str: Optional[str] = None,
    other_pdb_str: Optional[str] = None,
    reference_pdb_data: Optional[PDB_Data] = None,
    other_pdb_data: Optional[PDB_Data] = None,
    io: Any = None,
    soft_fail: bool = True,
    assert_is_isometry: bool = False,
    verbosity: int = 0,
    debug: bool = False,
    test_id: str = '',
    united_atom_fit: bool = UNITED_RMSD_FIT,
    exception_searching_keywords: List[str] = ALL_EXCEPTION_SEARCHING_KEYWORDS,
    **kwargs: Dict[str, Any]
) -> Tuple[PDB, RMSD, Alignment_Results]:
    assert reference_pdb_str is not None or reference_pdb_data is not None
    if reference_pdb_data is None:
        reference_pdb_data = pdb_data_for(reference_pdb_str, united_atom_fit=united_atom_fit, exception_searching_keywords=exception_searching_keywords)

    assert other_pdb_str is not None or other_pdb_data is not None
    if other_pdb_data is None:
        other_pdb_data = pdb_data_for(other_pdb_str, united_atom_fit=united_atom_fit, exception_searching_keywords=exception_searching_keywords)

    assert len(set([pdb_data.united_atom_fit for pdb_data in (reference_pdb_data, other_pdb_data)])) == 1, [pdb_data for pdb_data in (reference_pdb_data, other_pdb_data)]

    if debug:
        def pdb_writing_fct(alignment, file_name):
            pdb_path = join(DEBUG_DIR, test_id, file_name)

            if not exists(dirname(pdb_path)):
                mkdir(dirname(pdb_path))

            log.debug('PDB Writing Function: Dumping alignment to {0} (score={1})'.format(file_name, alignment.score))
            with open(pdb_path, 'w') as fh:
                fh.write(
                    aligned_pdb_str(
                        other_pdb_data.data,
                        alignment,
                        reference_pdb_data.united_atom_fit,
                    )
                )
    else:
        pdb_writing_fct = NULL_PDB_WRITING_FCT

    try:
        alignment = pointsOnPoints(
            [other_pdb_data.point_lists, reference_pdb_data.point_lists],
            flavour_lists=[other_pdb_data.flavour_lists, reference_pdb_data.flavour_lists],
            extra_points=other_pdb_data.extra_points_lists,
            verbosity=verbosity,
            soft_fail=soft_fail,
            assert_is_isometry=assert_is_isometry,
            pdb_writing_fct=pdb_writing_fct,
            **kwargs,
        )
    except (Topology_Error, AssertionError) as e:
        raise
#    except Exception as e:
#        alignment = FAILED_ALIGNMENT
#        if io:
#            print('<pre>{0}</pre>'.format(e), file=io)
#        if not soft_fail:
#            raise

    if alignment.aligned_points is None: # pragma: no cover
        if io:
            print('<pre>{0}</pre>'.format(alignment), file=io)
        final_aligned_pdb_str = other_pdb_data.pdb_str
        success = False
    else:
        try:
            final_aligned_pdb_str = aligned_pdb_str(other_pdb_data.data, alignment, reference_pdb_data.united_atom_fit)
            success = True
        except AssertionError:
            final_aligned_pdb_str = other_pdb_data.pdb_str
            success = False

    return (
        final_aligned_pdb_str,
        alignment.score,
        Alignment_Results(
            success=success,
            error_exc=format_exc(),
            permutation=alignment.final_permutation,
            atom_names_permutation=(
                {
                    other_pdb_data.data['atoms'][other_atom_id + 1]['symbol']: reference_pdb_data.data['atoms'][reference_atom_id + 1]['symbol']
                    for (other_atom_id, reference_atom_id) in alignment.final_permutation
                }
                if alignment.final_permutation is not None
                else None
            ),
            total_fit_points=None,
            kabsch_fit_point=None,
        ),
    )

def rmsd_matrix_for(list_of_pdb_str: List[str]) -> Any:
    list_of_pdb_data = list(map(
        pdb_data_for,
        list_of_pdb_str,
    ))

    def get_alignment_score(reference_pdb_data, other_pdb_data):
        try:
            alignment = align_pdb_on_pdb(
                reference_pdb_data=reference_pdb_data,
                other_pdb_data=other_pdb_data,
                soft_fail=True,
            )
            return alignment[1]
        except Topology_Error:
            return float('inf')

    rmsds = [
        [
            get_alignment_score(
                pdb_data_1,
                pdb_data_2,
            )
            for pdb_data_2 in list_of_pdb_data[i + 1:]
        ]
        for (i, pdb_data_1) in enumerate(list_of_pdb_data)
    ]

    return squareform(reduce(lambda acc, e: acc + e, rmsds, []))
