install: Vector.py charnley_rmsd/kabsch.py pmx testing
	echo 'Installed'

test: Vector.py charnley_rmsd/kabsch.py pmx testing
	make clean
	python test_rmsd.py --auto | tee log.out

test-pymol: test
	pymol -M testing/ethanol/*.pdb

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

clean:
	rm -rf testing/*
.PHONY: clean
