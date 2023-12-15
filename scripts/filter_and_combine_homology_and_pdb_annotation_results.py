'''

Filter off ppis-pdb-pdbchain combinations that are also seen in the homology modeling results
then combine results from homology modeling and PDB annotations.
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp
from get_interfacial_residues_using_sasa import get_stats_for_interface_residues_exogenous, get_stats_for_interface_residues_endogenous
import argparse

def filter_off_overlap_between_homology_and_pdb_annotations(homology_interfacial_residues_fi, pdb_annotation_residues_fi, remaining_fi, combined_fi):
	""" Filters off the pdb annotations that already exist in the homology modeling results (left with pdb interaction annotations that were not found from homology modeling)
	
	Args: 
		homology_interfacial_residues_fi: interfacial residue info from homology modeling results
		pdb_annotation_residues_fi: interfacial residue info from pdb annotation results
		remaining_fi: name of output file (to store the interfacial residues for pdb annotations after filtering off those that already exist in the homology results)
		combined_fi: name of output file with unique interfacial residues from homology modeling and pdb annotations combined
	"""
	list_of_homology_results = [] #store the homology modeling results
	l_check = []
	with open(combined_fi, 'wt') as combined:
		#target_prot	partner_prot	PDBid	target_chainid	partner_chainid	target_interfacial_residues
		combined.write("target_prot" + "\t" + "partner_prot" + "\t" + "PDBid" + "\t" + "target_chainid" + "\t" + "partner_chainid" + "\t" + "target_interfacial_residues" + "\n")
		with open(homology_interfacial_residues_fi, 'r') as hom_int_res:
			next(hom_int_res)
			for line in hom_int_res:
				info = line.split("\t")
				to_store = info[0:5]#store the target_prot	partner_prot	PDBid	target_chainid	partner_chainid
				to_store_first_three = info[0:3]
				l_check.append(to_store_first_three)
				list_of_homology_results.append(to_store)
				combined.write(line)

		with open(pdb_annotation_residues_fi, 'r') as pdb_anno_int_res:
			next(pdb_anno_int_res)
			with open(remaining_fi, 'w') as out_fi:
				for line in pdb_anno_int_res:
					line = line.strip("\n")
					info = line.split("\t")
					to_check = info[0:5]
					if to_check not in list_of_homology_results:
						#if info[0:3] in l_check:
							#print(line)
						out_fi.write(line +"\n")	
						combined.write(line + "\n")	




'''
Examples on how to run this script
### To run for the exogenous pdb anno interfacial residues: ###

python3 filter_and_combine_homology_and_pdb_annotation_results.py \
-m ../data/exo_structural_annotation/interface_residues/exo_hv_interface_residues.tsv \
-p ../data/pdb_annotations/interface_residues/pdb_hv_annotations_interface_residues.tsv \
-f ../data/pdb_annotations/interface_residues/filtered_pdb_hv_annotations_interface_residues.tsv \
-c ../data/combined_homology_and_pdb_anno_human_virus_interface_residues.tsv \
-t exo

homology_interfacial_residues_fi = osp.join(script_dir, '..', 'data', 'exo_structural_annotation', 'interface_residues', 'exo_hv_interface_residues.tsv')
pdb_interfacial_residues_fi = osp.join(path_to_pdb_annotations, 'interface_residues', 'pdb_hv_annotations_interface_residues.tsv')
pdb_interfacial_residues_filtered_fi = osp.join(path_to_pdb_annotations, 'interface_residues', 'filtered_pdb_hv_annotations_interface_residues.tsv')
combined_interfacial_residues_fi = osp.join(script_dir, '..', 'data', 'combined_homology_and_pdb_anno_human_virus_interface_residues.tsv')


### To run for the endogenous pdb anno interfacial residues: ###

python3 filter_and_combine_homology_and_pdb_annotation_results.py \
-m ../data/endo_structural_annotation/interface_residues/endo_human_interface_residues.tsv \
-p ../data/pdb_annotations/interface_residues/pdb_human_human_annotations_interface_residues.tsv \
-f ../data/pdb_annotations/interface_residues/filtered_pdb_human_human_annotations_interface_residues.tsv \
-c ../data/combined_homology_and_pdb_anno_human_human_interface_residues.tsv \
-t endo

homology_interfacial_residues_fi = osp.join(script_dir, '..', 'data', 'endo_structural_annotation', 'interface_residues', 'endo_human_interface_residues.tsv')
pdb_interfacial_residues_fi = osp.join(path_to_pdb_annotations, 'interface_residues', 'pdb_human_human_annotations_interface_residues.tsv')
pdb_interfacial_residues_filtered_fi = osp.join(path_to_pdb_annotations, 'interface_residues', 'filtered_pdb_human_human_annotations_interface_residues.tsv')
combined_interfacial_residues_fi = osp.join(script_dir, '..', 'data', 'combined_homology_and_pdb_anno_human_human_interface_residues.tsv')
'''
def main():
	script_dir = osp.dirname(__file__)
	
	parser = argparse.ArgumentParser()

	#read in pdb_chain_pairs file and output file as arguments! 
	#So that we can use this script for other input files too!
	parser.add_argument('-m', '--homology_modeling_int_res_fi')
	parser.add_argument('-p', '--pdb_int_res_fi')
	parser.add_argument('-f', '--filtered_pdb_int_res_fi')
	parser.add_argument('-c', '--combined_hom_and_pdb_int_res_fi')
	parser.add_argument('-t', '--endo_or_exo')

	args  = parser.parse_args()
	
	homology_interfacial_residues_fi = args.homology_modeling_int_res_fi
	pdb_interfacial_residues_fi = args.pdb_int_res_fi
	pdb_interfacial_residues_filtered_fi = args.filtered_pdb_int_res_fi
	combined_interfacial_residues_fi = args.combined_hom_and_pdb_int_res_fi
	endo_or_exo = args.endo_or_exo	
	
	#filter off overlap
	filter_off_overlap_between_homology_and_pdb_annotations(homology_interfacial_residues_fi , pdb_interfacial_residues_fi, pdb_interfacial_residues_filtered_fi, combined_interfacial_residues_fi)
	print("##### Stats for pdb annotations (after filtering off homology modeling overlap) #####")
	
	if endo_or_exo == 'exo':
		get_stats_for_interface_residues_exogenous(pdb_interfacial_residues_filtered_fi) #these are the stats AFTER filtering off the overlap btwn homology modeling and pdb annotations
	elif endo_or_exo == 'endo':
		get_stats_for_interface_residues_endogenous(pdb_interfacial_residues_filtered_fi) #these are the stats AFTER filtering off the overlap btwn homology modeling and pdb annotations
	
	#After combining homology modeling results and pdb annotations (unique only)
	print("##### Stats for filtered pdb annotations + homology modeling (combined) #####")
	
	if endo_or_exo == 'exo':
		get_stats_for_interface_residues_exogenous(combined_interfacial_residues_fi) #these are the stats for the combined filtered pdb annotation and homology modeling results
	elif endo_or_exo == 'endo':
		get_stats_for_interface_residues_endogenous(combined_interfacial_residues_fi) #these are the stats for the combined filtered pdb annotation and homology modeling results

	


if __name__ == '__main__':
	main()