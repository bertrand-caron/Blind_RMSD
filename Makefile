test: Vector.py charnley_rmsd/kabsch.py testing/benzene1.pdb testing/benzene2.pdb testing/ethanol1.pdb testing/ethanol2.pdb testing/methanol1.pdb testing/methanol2.pdb pmx testing/sulfuric_acid1.pdb testing/sulfuric_acid2.pdb
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

testing/benzene1.pdb: testing
	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=21364&outputType=top&dbfile=pdb_fromuser" -O $@

# This molecule is actually restricted access only ...
#testing/benzene2.pdb: testing
#	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=8479&outputType=top&dbfile=pdb_fromuser" -O $@

testing/ethanol1.pdb: testing
	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=22354&outputType=top&dbfile=pdb_fromuser" -O $@

testing/ethanol2.pdb: testing
	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=23009&outputType=top&dbfile=pdb_fromuser" -O $@

testing/sulfuric_acid1.pdb: testing
	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=1745&outputType=top&dbfile=pdb_fromuser" -O $@

testing/sulfuric_acid2.pdb: testing
	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=20681&outputType=top&dbfile=pdb_fromuser" -O $@

#testing/methanol1.pdb: testing
#	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=15607&outputType=top&dbfile=pdb_fromuser" -O $@

#testing/methanol2.pdb: testing
#	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=19901&outputType=top&dbfile=pdb_fromuser" -O $@

pmx:
		echo -e "Please install pmx\n\n"
.PHONY: pmx
