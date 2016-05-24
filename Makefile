PYTHONPATH = PYTHONPATH="/home/$$USER/ATB:/home/$$USER/ATB-Dependencies"

PYTHON_EXEC = $(PYTHONPATH) python

install: helpers/Vector.py lib/charnley_rmsd/kabsch.py pmx testing
	echo 'Installed'

clean_atb_duplicates: helpers/Vector.py lib/charnley_rmsd/kabsch.py pmx testing
	make clean
	$(PYTHON_EXEC) tasks/$@.py --auto --verbosity 1 | tee log_2.out
.PHONY: clean_atb_duplicates

test: helpers/Vector.py lib/charnley_rmsd/kabsch.py pmx testing
	$(PYTHON_EXEC) test.py
.PHONY: test

test-pymol: test
	pymol -M testing/ethanol/*.pdb

helpers/Vector.py: lib/biopython
	ln -s ../$</Bio/PDB/Vector.py $@
	touch $@

lib/biopython:
	cd lib && git clone https://github.com/bertrand-caron/biopython.git

lib/charnley_rmsd/kabsch.py: lib/charnley_rmsd
	if [[ -e $@ ]]; then unlink $@; fi
	ln -s calculate_rmsd $@

lib/charnley_rmsd:
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

errors:
	$(PYTHONPATH) pylint -E *.py helpers/*.py tasks/*.py
.PHONY: errors
