from traceback import format_exc
from collections import namedtuple
from os.path import abspath, join, dirname, exists
from os import mkdir
from scipy.spatial.distance import squareform

from chemEquivalency.calcChemEquivalency import partial_mol_data_for_pdbstr
from Blind_RMSD.helpers.moldata import flavour_list, point_list, aligned_pdb_str, united_hydrogens_point_list
from Blind_RMSD.align import pointsOnPoints, FAILED_ALIGNMENT, NULL_PDB_WRITING_FCT
from Blind_RMSD.helpers.exceptions import Topology_Error

UNITED_RMSD_FIT = True

PDB_Data = namedtuple('PDB_Data', 'data, point_lists, flavour_lists, extra_points_lists, pdb_str')

FAILED_ALIGNMENT_PDB_STR = ''

Alignment_Results = namedtuple('Alignment_Results', 'success, error_exc')

DEBUG_DIR = join(dirname(abspath(__file__)), 'debug')

def pdb_data_for(pdb_str):
    data = partial_mol_data_for_pdbstr(pdb_str).__dict__

    return PDB_Data(
        data=data,
        point_lists = point_list(data, UNITED_RMSD_FIT),
        flavour_lists = flavour_list(data, UNITED_RMSD_FIT),
        extra_points_lists = united_hydrogens_point_list(data, UNITED_RMSD_FIT),
        pdb_str = pdb_str,
    )

def align_pdb_on_pdb(reference_pdb_str=None, other_pdb_str=None, reference_pdb_data=None, other_pdb_data=None, io=None, soft_fail=True, assert_is_isometry=False, verbosity=0, debug=False, test_id=''):
    assert reference_pdb_str is not None or reference_pdb_data is not None
    if reference_pdb_data is None:
        reference_pdb_data = pdb_data_for(reference_pdb_str)

    assert other_pdb_str is not None or other_pdb_data is not None
    if other_pdb_data is None:
        other_pdb_data = pdb_data_for(other_pdb_str)

    if debug:
        def pdb_writing_fct(alignment, file_name):
            pdb_path = join(DEBUG_DIR, test_id, file_name)

            if not exists(dirname(pdb_path)):
                mkdir(dirname(pdb_path))

            print '    PDB Writing Function: Dumping alignment to {0} (score={1})'.format(file_name, alignment.score)
            with open(pdb_path, 'w') as fh:
                fh.write(
                    aligned_pdb_str(
                        other_pdb_data.data,
                        alignment,
                        UNITED_RMSD_FIT,
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
        )
    except (Topology_Error, AssertionError) as e:
        raise
#    except Exception as e:
#        alignment = FAILED_ALIGNMENT
#        if io:
#            print >> io, '<pre>{0}</pre>'.format(e)
#        if not soft_fail:
#            raise

    if alignment.aligned_points is None: # pragma: no cover
        if io:
            print >> io, '<pre>{0}</pre>'.format(alignment)
        final_aligned_pdb_str = other_pdb_data.pdb_str
        success = False
    else:
        try:
            final_aligned_pdb_str = aligned_pdb_str(other_pdb_data.data, alignment, UNITED_RMSD_FIT)
            success = True
        except AssertionError:
            final_aligned_pdb_str = other_pdb_data.pdb_str
            success = False

    return (
        final_aligned_pdb_str,
        alignment.score,
        Alignment_Results(
            success,
            format_exc(),
        ),
    )

def rmsd_matrix_for(list_of_pdb_str):
    list_of_pdb_data = map(
        pdb_data_for,
        list_of_pdb_str,
    )

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
