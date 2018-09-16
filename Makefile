PYTHON_EXEC = PYTHONPATH=$(PYTHONPATH) python3

DEPENDENCIES = ../chemistry_helpers ../chemical_equivalence ../atb_outputs helpers/Vector.py testing

../chemistry_helpers:
	cd .. && git clone https://github.com/bertrand-caron/chemistry_helpers.git

../chemical_equivalence:
	cd .. && git clone https://github.com/ATB-UQ/chemical_equivalence.git

../atb_outputs:
	cd .. && git clone https://github.com/bertrand-caron/atb_outputs

install: $(DEPENDENCIES)
	echo 'Installed'

clean_atb_duplicates: install
	make clean
	$(PYTHON_EXEC) tasks/$@.py --auto --verbosity 1 | tee log_2.out
.PHONY: clean_atb_duplicates

test: install
	$(PYTHON_EXEC) test.py
.PHONY: test

debug: install
	$(PYTHON_EXEC) test_interactive.py --reference data/a.pdb --other data/c.pdb --N 10
.PHONY: debug

test-pymol: test
	pymol -M testing/ethanol/*.pdb

helpers/Vector.py:
	make lib/biopython
	ln -s ../$</Bio/PDB/Vector.py $@
	touch $@

lib/biopython:
	cd lib && git clone https://github.com/bertrand-caron/biopython.git

testing:
	mkdir $@

clean:
	rm -rf testing/*
.PHONY: clean

errors:
	$(PYTHONPATH) pylint -E *.py helpers/*.py tasks/*.py
.PHONY: errors
