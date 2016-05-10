from collections import namedtuple

from chemEquivalency.calcChemEquivalency import partial_mol_data_for_pdbstr
from Blind_RMSD.helpers.moldata import flavour_list, element_list, point_list, aligned_pdb_str, united_hydrogens_point_list
from Blind_RMSD.align import pointsOnPoints, FAILED_ALIGNMENT

UNITED_RMSD_FIT = True

PDB_Data = namedtuple('PDB_Data', 'data, point_lists, element_lists, flavour_lists, extra_points_lists, pdb_str')

FAILED_ALIGNMENT_PDB_STR = ''

def pdb_data_for(pdb_str):
    data = partial_mol_data_for_pdbstr(pdb_str).__dict__

    return PDB_Data(
        data=data,
        point_lists = point_list(data, UNITED_RMSD_FIT),
        element_lists = element_list(data, UNITED_RMSD_FIT),
        flavour_lists = flavour_list(data, UNITED_RMSD_FIT),
        extra_points_lists = united_hydrogens_point_list(data, UNITED_RMSD_FIT),
        pdb_str = pdb_str,
    )

def align_pdb_on_pdb(reference_pdb_str=None, other_pdb_str=None, reference_pdb_data=None, other_pdb_data=None, io=None, silent=True, soft_fail=True):
    assert reference_pdb_str is not None or reference_pdb_data is not None
    if reference_pdb_data is None:
        reference_pdb_data = pdb_data_for(reference_pdb_str)

    assert other_pdb_str is not None or other_pdb_data is not None
    if other_pdb_data is None:
        other_pdb_data = pdb_data_for(other_pdb_str)

    try:
        alignment = pointsOnPoints(
            [other_pdb_data.point_lists, reference_pdb_data.point_lists],
            element_lists=[other_pdb_data.element_lists, reference_pdb_data.element_lists],
            flavour_lists=[other_pdb_data.flavour_lists, reference_pdb_data.flavour_lists],
            extra_points=other_pdb_data.extra_points_lists,
            silent=silent,
            soft_fail=soft_fail,
        )
    except Exception, e:
        alignment = FAILED_ALIGNMENT
        if io:
            print >> io, '<pre>{0}</pre>'.format(e)
        if not soft_fail:
            raise

    if alignment.aligned_points is None: # pragma: no cover
        if io:
            #print >> io, '<p>Aligment Failed for PDB {0}.</p>'.format(i+1)
            print >> io, '<pre>{0}</pre>'.format(alignment)
        final_aligned_pdb_str = other_pdb_data.pdb_str
    else:
        final_aligned_pdb_str = aligned_pdb_str(other_pdb_data.data, alignment, UNITED_RMSD_FIT)

    return (
        final_aligned_pdb_str,
        alignment.score,
    )
