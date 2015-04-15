test: Vector.py charnley_rmsd/kabsch.py pmx
	python test_rmsd.py

test-pymol:
	pymol -M testing/methanol*

Vector.py: biopython
	ln -s biopython/Bio/PDB/$@ $@
	touch $@

biopython:
	git clone https://github.com/bertrand-caron/biopython.git

charnley_rmsd/kabsch.py: charnley_rmsd
	if [[ -e $@ ]]; then unlink $@; fi
	ln -s calculate_rmsd $@

charnley_rmsd:
	git clone https://github.com/charnley/rmsd.git $@
	touch $@/__init__.py

testing:
	mkdir $@

pmx:
		echo -e "Please install pmx\n\n"
.PHONY: pmx
