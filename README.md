# Evolution of human-virus PPI interfaces

Building a human-virus structural interaction network (SIN), resolving protein-protein interaction (PPI) interfaces, mapping and categorizing interfaces, and computing the evolutionary rates of different PPI interfaces.

## Getting started
1. Clone this repository.  
`$ git clone 'https://github.com/wanchunsu/HV_PPI_interface_evolution.git'`

2. Navigate to the downloaded folder.  
`$ cd HV_PPI_interface_evolution`

3. Install requirements listed in [`dependencies.txt`](https://github.com/wanchunsu/HV_PPI_interface_evolution/blob/main/dependencies.txt). 

## Usage
1. Navigate to the `scripts` directory.  
`$ cd scripts`

2. Follow the instructions in the [`instructions.bash`](https://github.com/wanchunsu/HV_PPI_interface_evolution/blob/main/scripts/instructions.bash) script. Below are the main steps in this pipeline:
	1. Build a human-virus SIN and resolve interaction interfaces.
	2. Map interfacial residues back onto human target proteins and [`categorize them`](https://github.com/wanchunsu/HV_PPI_interface_evolution/blob/main/scripts/categorize_interfacial_residues.py).
	3. Determine the evolutionary rates of different PPI interfaces.

**Note**: The working directory is `scripts` and all data and results are stored in the `data` directory. 

## Software
* **Python 3.9.12**
* **BLAST 2.13.0** (see https://www.ncbi.nlm.nih.gov/books/NBK52640/ for installation guide)
* **DSSP (mkdssp)** (installed in the [`install_dssp.bash`](https://github.com/wanchunsu/HV_PPI_interface_evolution/blob/main/scripts/install_dssp.bash) script)
* **MAFFT v7.490** (installed in the [`construct_msas.bash`](https://github.com/wanchunsu/HV_PPI_interface_evolution/blob/main/scripts/construct_msas.bash) script)
* **PAML** (installed in the [`get_binned_dnds.bash`](https://github.com/wanchunsu/HV_PPI_interface_evolution/blob/main/scripts/get_binned_dnds.bash) script)
* **PAL2NAL v14** (installed in the [`calc_site_specific_dnds.bash`](https://github.com/wanchunsu/HV_PPI_interface_evolution/blob/main/scripts/calc_site_specific_dnds.bash) script)
* **HyPhy 2.5.40** (installed in the [`calc_site_specific_dnds.bash`](https://github.com/wanchunsu/HV_PPI_interface_evolution/blob/main/scripts/calc_site_specific_dnds.bash) script)

## Databases
* **IntAct** (https://www.ebi.ac.uk/intact/home)
* **PDB** (https://www.rcsb.org/)
* **UniProt** (https://www.uniprot.org/)
* **Ensembl** (https://useast.ensembl.org/index.html)
* **RefSeq** (https://www.ncbi.nlm.nih.gov/refseq/)
* **TimeTree** (http://timetree.org/)

## Contact
Wan-Chun Su wan.su@mail.mcgill.ca









