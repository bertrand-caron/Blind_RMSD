test: Vector.py charnley_rmsd/kabsch.py pmx testing ../atb_py
	make clean
	python test_rmsd.py --auto

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

../atb_py:
	cd .. && git clone ssh://git@scmb-gitlab.biosci.uq.edu.au:2023/ATB/atb_py.git && cd atb_py && python setup.py install --user
