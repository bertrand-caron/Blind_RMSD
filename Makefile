test: Vector.py testing/benzene1.pdb testing/benzene2.pdb testing/ethanol1.pdb testing/ethanol2.pdb testing/methanol1.pdb testing/methanol2.pdb pmx
	python test_rmsd.py
	#pymol -M testing/methanol*

Vector.py: biopython
	ln -s biopython/Bio/PDB/$@ $@
	touch $@

biopython:
	git clone https://github.com/bertrand-caron/biopython.git

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

testing/methanol1.pdb: testing
	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=15607&outputType=top&dbfile=pdb_fromuser" -O $@

testing/methanol2.pdb: testing
	wget "http://compbio.biosci.uq.edu.au/atb/download.py?molid=19901&outputType=top&dbfile=pdb_fromuser" -O $@

pmx:
		echo -e "Please install pmx\n\n"
.PHONY: pmx
