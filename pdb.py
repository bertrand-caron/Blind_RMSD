from collections import namedtuple

from chemEquivalency.calcChemEquivalency import partial_mol_data_for_pdbstr
from Blind_RMSD.helpers.moldata import flavour_list, element_list, point_list, aligned_pdb_str, united_hydrogens_point_list
from Blind_RMSD.align import pointsOnPoints, FAILED_ALIGNMENT

UNITED_RMSD_FIT = True

PDB_Data = namedtuple('PDB_Data', 'data, point_lists, element_lists, flavour_lists, extra_points_lists')

def pdb_data_for(pdb_str):
    data = partial_mol_data_for_pdbstr(pdb_str).__dict__

    return PDB_Data(
        data=data,
        point_lists = point_list(data, UNITED_RMSD_FIT),
        element_lists = element_list(data, UNITED_RMSD_FIT),
        flavour_lists = flavour_list(data, UNITED_RMSD_FIT),
        extra_points_lists = united_hydrogens_point_list(data, UNITED_RMSD_FIT),
    )

def align_pdb_on_pdb(reference_pdb_str=None, other_pdb_str=None, reference_pdb_data=None, other_pdb_data=None, io=None):
    assert reference_pdb_str is not None or reference_pdb_data is not None
    if reference_pdb_data is None:
        reference_pdb_data = pdb_data_for(reference_pdb_str)

    assert other_pdb_str is not None or other_pdb_data is not None
    if other_pdb_data is None:
        other_pdb_data = pdb_data_for(other_pdb_str)

    try:
        alignment = pointsOnPoints(
            (other_pdb_data.point_lists, reference_pdb_data.point_lists),
            element_lists=(other_pdb_data.element_lists, reference_pdb_data.element_lists),
            flavour_lists=(other_pdb_data.flavour_lists, reference_pdb_data.flavour_lists),
            soft_fail=True,
            extra_points=other_pdb_data.extra_points_lists,
        )
    except Exception, e:
        alignment = FAILED_ALIGNMENT
        if io:
            print >> io, '<pre>{0}</pre>'.format(e)
    if alignment.aligned_points is None: # pragma: no cover
        if io:
            #print >> io, '<p>Aligment Failed for PDB {0}.</p>'.format(i+1)
            print >> io, '<pre>{0}</pre>'.format(alignment)
        raise Exception('Alignment failed.')
    else:
        return (
            aligned_pdb_str(other_pdb_data.data, alignment, UNITED_RMSD_FIT),
            alignment.score,
        )
