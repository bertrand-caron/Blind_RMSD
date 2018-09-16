[![DOI](https://zenodo.org/badge/148723006.svg)](https://zenodo.org/badge/latestdoi/148723006)

# Dependencies

* `bertrand-caron/biopython`: Fork of Biopython (only used for its Vector class).
* `ATB-UQ/chemical_equivalence`: Predicting chemical equivalence of atoms in molecules using graph automorphisms (nauty) and chemical symmetry constraints.
* `bertrand-caron/atb_outputs`: Automated Topology Builder (ATB) molecule data object and output modules (PDB, YML, PICKLE, LGF, GRAPH)
* `bertrand-caron/chemistry_helpers`: Various chemistry related helpers.
* `GNU Make`: Included with most Linux distributions.

# Installation

* Run `make install` to automatically download all dependencies.
Note that `chemical_equivalence`, `atb_outputs` and `chemistry_helpers` will be cloned in the parent directory of `Blind_RMSD`.
Please add the parent directory to your shell's `PYTHONPATH` variable.
If you are using `bash` (Bourne Again Shell), run the following line:
`export PYTHONPATH=$(dirname $PWD):$PYTHONPATH`

# Usage

The `test.py` script contains code to run the Blind RMSD algorithm on a series of test cases found in the `data/` directory.
You can run it by running: `make test` (once you have successfully set up your `PYTHONPATH` variable).

# Troubleshooting

* If you get a `ModuleNotFoundError`, this most likely mean there is a problem with your PYTHONPATH.
Consider what happens in the following example when setting the `PYTHONPATH` variable to an empty string:

```
bertrand-caron@iMac Blind_RMSD  [11:13:47] $ PYTHONPATH="" make test
echo 'Installed'
Installed
PYTHONPATH="" python3 test.py
Traceback (most recent call last):
  File "test.py", line 3, in <module>
    from Blind_RMSD.pdb import pdb_data_for, align_pdb_on_pdb, ALL_EXCEPTION_SEARCHING_KEYWORDS
ModuleNotFoundError: No module named 'Blind_RMSD'
make: *** [test] Error 1
```

# Citation / Attribution

To cite this work, please use the following [Zenodo DOI](https://zenodo.org/badge/latestdoi/148723006).
