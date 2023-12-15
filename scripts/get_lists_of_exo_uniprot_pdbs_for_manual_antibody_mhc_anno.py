'''

Obtain list of human target proteins with exogenous PPI data 
and their associated PDB structures 
============================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp

def get_lists_of_pdb_and_uniprot_for_manual_annotation(exogenous_ppis_file, file_to_annotate_pdb, file_to_annotate_uniprot_ids):
	list_of_pdbs_to_anno = []
	list_of_uniprot_to_anno = []
	with open(exogenous_ppis_file) as exo:
		next(exo)
		for line in exo:
			pdb = line.split("\t")[2]
			uniprot = line.split("\t")[0]
			if pdb not in list_of_pdbs_to_anno:
				list_of_pdbs_to_anno.append(pdb)
			if uniprot not in list_of_uniprot_to_anno:
				list_of_uniprot_to_anno.append(uniprot)

	with open(file_to_annotate_pdb, 'w') as o_pdb:
		for p in list_of_pdbs_to_anno:
			o_pdb.write(p +"\t\n")

	with open(file_to_annotate_uniprot_ids, 'w') as o_uniprot:
		for u in list_of_uniprot_to_anno:
			o_uniprot.write(u +"\t\n")


def main():
	script_dir = osp.dirname(__file__)
	antibody_and_mhc_complexes_dir = osp.join(script_dir, '..', 'data', 'antibody_mhc_complex_annotations')
	if not osp.exists(antibody_and_mhc_complexes_dir):
		os.mkdir(antibody_and_mhc_complexes_dir)
	exogenous_ppis_file = osp.join(script_dir, '..', 'data', 'combined_homology_and_pdb_anno_human_virus_interface_residues.tsv')
	file_to_annotate_pdb = osp.join(antibody_and_mhc_complexes_dir, 'exogenous_ppi_pdbs_with_antibody_mhc_annotations.tsv')
	file_to_annotate_uniprot_ids = osp.join(antibody_and_mhc_complexes_dir, 'exogenous_human_targs_with_antibody_mhc_annotations.tsv')
	get_lists_of_pdb_and_uniprot_for_manual_annotation(exogenous_ppis_file, file_to_annotate_pdb, file_to_annotate_uniprot_ids)

if __name__ == '__main__':
	main()