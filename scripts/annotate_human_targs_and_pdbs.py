'''

Automatically annotate PDB structures where all mapped 
human targets are non-antibody/mhc
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp

def make_dict_of_pdbs_with_corresponding_human_targs(all_exo_ppis):
	dict_of_pdbs_with_corr_targs = {}
	with open(all_exo_ppis) as exo:
		next(exo)
		for line in exo:
			human_targ = line.split("\t")[0]
			pdb = line.split("\t")[2]
			if pdb not in dict_of_pdbs_with_corr_targs:
				dict_of_pdbs_with_corr_targs[pdb] = [human_targ]
			else:
				if human_targ not in dict_of_pdbs_with_corr_targs[pdb]:
					dict_of_pdbs_with_corr_targs[pdb].append(human_targ)

	return dict_of_pdbs_with_corr_targs


def make_dict_of_human_targ_antibody_mhc_status(annotated_human_targs_file):
	dict_of_human_targ_anno = {}
	with open(annotated_human_targs_file) as anno:
		for line in anno:
			line = line.strip("\n")
			info = line.split("\t")

			uniprotid = info[0]
			annot = info[1]

			# verify that annotation was done correctly
			if annot == '': #check for missing annotation
				print(line)
			if annot != 'T' and annot != 'F': #check for wrong annotation characters
				print(line)

			dict_of_human_targ_anno[uniprotid] = annot
	return dict_of_human_targ_anno


def annotate_pdbs_for_non_anitbody_mhc_targ_prots(dict_of_pdbs_with_corr_targs, dict_of_human_targ_anno,  output_anno):
	
	dict_of_anno_pdb = {}

	for pdb in dict_of_pdbs_with_corr_targs:
		
		list_of_human_targ_annos_for_this_pdb = [dict_of_human_targ_anno[t] for t in dict_of_pdbs_with_corr_targs[pdb]]
		if all(x == list_of_human_targ_annos_for_this_pdb[0] for x in list_of_human_targ_annos_for_this_pdb):
			if list_of_human_targ_annos_for_this_pdb[0] == 'F': #if all human targets in this pdb are non-antibody/mhc, then this pdb must also be non antibody/mhc
				dict_of_anno_pdb[pdb] = 'F'
			else: #all human targets in this pdb are antibody/mhc proteins
				dict_of_anno_pdb[pdb] = ''
				continue # leave this open for manual annotation
				"""
				could be 1 of 2 cases
				1) all human targets are targetting virus protein (annotate with 'T')
				2) virus protein is targetting human antibody/mhc proteins (annotate with 'F')
				"""
		else: # human target proteins mapping to this pdb have different annotations (some are antibody/mhc, some are non-antibody/mhc)
			dict_of_anno_pdb[pdb] = ''
			print(pdb)
			continue #leave this open for manual annotation (annotate as 'D' or 'F')
			# 'D': Found that for some of these cases, it is the virus protein targetting a host receptor and host proteins binding to the virus protein at the same time
			# OR 
			# 'F': The virus protein is targeting both an antibody/MHC (e.g. for degradation) and a non-antibody/non-MHC protein

	with open(output_anno, 'w') as o:
		for p in dict_of_anno_pdb:
			o.write(p + "\t" + dict_of_anno_pdb[p] + "\n")

def main():
	script_dir = osp.dirname(__file__)
	all_exo_ppis = osp.join(script_dir, '..', 'data', 'combined_homology_and_pdb_anno_human_virus_interface_residues.tsv')
	annotated_human_targs_file = osp.join(script_dir, '..', 'data', 'antibody_mhc_complex_annotations',  'exogenous_human_targs_with_antibody_mhc_annotations.tsv')
	output_anno = osp.join(script_dir, '..', 'data', 'antibody_mhc_complex_annotations',  'exogenous_ppi_pdbs_with_antibody_mhc_annotations.tsv')

	# make dictionary to store pdbs (key) and their corresponding human targets in a list (value)
	dict_of_pdbs_with_corr_targs = make_dict_of_pdbs_with_corresponding_human_targs(all_exo_ppis)

	#make dictionary to store the antibody/mhc annotation (value) of each human target protein (key)
	dict_of_human_targ_anno = make_dict_of_human_targ_antibody_mhc_status(annotated_human_targs_file)
	
	#output annotations for the false anitbody/mhc pdb complexes
	annotate_pdbs_for_non_anitbody_mhc_targ_prots(dict_of_pdbs_with_corr_targs, dict_of_human_targ_anno, output_anno)
		
	# Next step is to manually annotate the rest

if __name__ == '__main__':
	main()