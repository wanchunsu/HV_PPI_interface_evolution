'''

Get % mismatch between human target proteins and their orthologs on a per-protein basis
=======================================================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp
import json
import argparse

def get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, residue_results):
	""" Get the fraction divergence between the target protein and it ortholog (store in dict)
		
	Args:
		alignment_matches_btwn_human_and_close_species_orthologs: dictionary of alignments between target protein and its ortholog 
		residue_results: file of residues (either exo-specific, endo-specific, mimicry, surface, or buried)
	Returns:
		divergence_dict: fraction divergence per target protein

	"""
	divergence_dict = {}
	with open(residue_results) as int_res:
		next(int_res)
		for line in int_res:
			line = line.strip("\n")
			targ_prot = line.split("\t")[0]
			if line.split("\t")[1] == '': #this could happen for the surface residues (some proteins don't have surface residues, so will be an empty string '')
				continue
			int_res = line.split("\t")[1].split(",")
			#int_res = [int(r) for r in int_res]
			total_num_int_res = len(int_res)
			total_mismatches = 0
			total_gaps = 0
			total_substitutions = 0
			for i in int_res:

				if i in alignment_matches_btwn_human_and_close_species_orthologs[targ_prot]:
					if alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == 0: #mismatch
						total_mismatches += 1
						total_substitutions += 1
					elif alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == -1: #maps to a gap
						total_mismatches += 1
						total_gaps += 1
				else: # targ prot res not in alignment (i.e. doesn't map to anything on the ortholog -- treat as a gap too)
					total_mismatches += 1
					total_gaps += 1

			fraction_divergence = total_mismatches / total_num_int_res 
			fraction_gaps = total_gaps / total_num_int_res 
			fraction_substitutions = total_substitutions / total_num_int_res 

			divergence_dict[targ_prot] = {"divergence": round(fraction_divergence*100, 3), "substitutions": round(fraction_substitutions*100, 3), "gaps": round(fraction_gaps*100, 3), "total_num_int_res": total_num_int_res}

	return divergence_dict


def get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, residue_results):
	""" Get list of 0 and 1s representing substitutions for residues in a residue category, or representing gaps for residues in a residue category
		Used for boostrap method for calculating standard error for average divergence values for each residue type
	"""
	dict_of_binary_vals = {} 
	with open(residue_results) as int_res:
		next(int_res)
		for line in int_res:
			line = line.strip("\n")

			targ_prot = line.split("\t")[0]
			if line.split("\t")[1] == '': #this could happen for the surface residues (some proteins don't have surface residues, so will be an empty string '')
				continue
			int_res = line.split("\t")[1].split(",")
			#int_res = [int(r) for r in int_res]

			list_of_gaps_and_non_gaps = [] #store 1 if gap else store 0
			list_of_subs_and_non_subs = [] #store 1 if substituion else store 0
			list_of_div_and_non_divs = [] #store 1 if is a difference (either gap/substitution) else store 0

			for i in int_res:

				if i in alignment_matches_btwn_human_and_close_species_orthologs[targ_prot]:
					if alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == 0: #maps to a mismatch
						
						list_of_subs_and_non_subs.append(1) #substitution so add a 1 to substitution list
						list_of_gaps_and_non_gaps.append(0) # not a gap so add a 0 to gap list
						list_of_div_and_non_divs.append(1) #is a difference so add 1 to divergence list

					elif alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == -1: #maps to a gap
						list_of_subs_and_non_subs.append(0) # not substitution so add a 0 to substitution list
						list_of_gaps_and_non_gaps.append(1) # gap so add a 1 to gap list
						list_of_div_and_non_divs.append(1) #is a difference so add 1 to divergence list
					
					else: #match, neither gap nor substitution, so add 0 to both lists
						list_of_subs_and_non_subs.append(0) # not a substitution so add a 0 to substitution list
						list_of_gaps_and_non_gaps.append(0) # not gap so add a 0 to gap list
						list_of_div_and_non_divs.append(0) # not a difference so add 0 to divergence list

				else: # targ prot res not in alignment (i.e. doesn't map to anaything on the ortholog -- treat as a gap too)
					list_of_subs_and_non_subs.append(0) # not substitution so add a 0 to substitution list
					list_of_gaps_and_non_gaps.append(1) # gap so add a 1 to gap list
					list_of_div_and_non_divs.append(1) #is a difference so add 1 to divergence list

			dict_of_binary_vals[targ_prot] = {"div_list": list_of_div_and_non_divs, "subs_list": list_of_subs_and_non_subs, "gaps_list": list_of_gaps_and_non_gaps}
	return  dict_of_binary_vals

def get_per_protein_avg_divergence(list_of_divergence_values, subs_or_gaps):
	dict_of_relevant_data = {}
	residue_types = ['exo', 'endo', 'mimicry', 'surface_res', 'buried_res', 'all_exo', 'all_endo']
	for r in residue_types:
		for org in list_of_divergence_values:
			for targ_prot in list_of_divergence_values[org][r]:
				if targ_prot not in dict_of_relevant_data:
					dict_of_relevant_data[targ_prot] = {}
				if r not in dict_of_relevant_data[targ_prot]:
					dict_of_relevant_data[targ_prot][r] = []

				dict_of_relevant_data[targ_prot][r].append(list_of_divergence_values[org][r][targ_prot][subs_or_gaps])

	dict_of_averages = {}
	for targ_prot in dict_of_relevant_data:
		dict_of_averages[targ_prot] = {}
		for r in dict_of_relevant_data[targ_prot]:
			avg_div = sum(dict_of_relevant_data[targ_prot][r])/len(dict_of_relevant_data[targ_prot][r])
			dict_of_averages[targ_prot][r] = avg_div

	return dict_of_averages


def main():

	script_dir = osp.dirname(__file__)
	
	# Get command line arguments for specific folder and file
	
	parser = argparse.ArgumentParser()

	parser.add_argument('-f', '--folder')
	
	args = parser.parse_args()
	folder_name = args.folder #'v_target_h' 
	
	all_exo = osp.join(script_dir, '..', 'data', 'categorized_interfacial_residues', folder_name, 'exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv')
	
	# Files with exo-specific, endo-specific, and mimicry interfacial residues as well as surface and buried residues
	exo_specific_file =  osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_specific_interfacial_residues.tsv")
	endo_specific_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_specific_interfacial_residues.tsv")
	mimicry_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "mimicked_interfacial_residues.tsv")
	surface_res_file = osp.join(script_dir, '..', 'data', "surface_residues", folder_name, "valid_surface_residues_on_uniprot_seq.tsv")
	buried_res_file = osp.join(script_dir, '..', 'data', "buried_residues", folder_name, "valid_buried_residues_on_uniprot_seq.tsv")

	# Files for all exo and all endo
	all_exo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")
	all_endo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")


	output_fraction_divergence_json = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', folder_name, 'per_protein_divergence_with_close_species.json')
	
	if not osp.exists(osp.dirname(output_fraction_divergence_json)):
		os.mkdir(osp.dirname(output_fraction_divergence_json))

	list_of_divergence_values = {}
	
	output_average_divergence_json = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', folder_name, 'per_protein_avg_divergence.json')

	print(f"###   Divergence values for {folder_name}: ###  \n")
	for blast_res in os.listdir(osp.join(script_dir, '..', 'data', 'blast', 'results', 'reciprocal_best_hit')):
		
		if blast_res.endswith("against_human.blast.out"): 
			
			species_name = blast_res.split("_")[0]
			print(f"### Calculating divergence of interfacial and non-interfacial residues ({species_name} comparison) ###")

			
			orthologs_folder = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', folder_name, 'rbh', species_name)

			alignment_matches_out_fi = osp.join(orthologs_folder, "rbh_ortholog_alignment_matches.json")
			
			#Get alignment between the human and close_species RBH orthologs
			alignment_matches_btwn_human_and_close_species_orthologs = {}
			with open(alignment_matches_out_fi) as a:
				alignment_matches_btwn_human_and_close_species_orthologs = json.load(a)

			list_of_divergence_values[species_name] = {}

			############## Calculating divergence for each protein individually (each protein will have its own divergence value) ##############
			list_of_divergence_values[species_name]['exo'] = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, exo_specific_file)
			list_of_divergence_values[species_name]['endo'] = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, endo_specific_file)
			list_of_divergence_values[species_name]['mimicry'] = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, mimicry_file)
			list_of_divergence_values[species_name]['surface_res'] = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, surface_res_file)
			list_of_divergence_values[species_name]['buried_res'] = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, buried_res_file)

			list_of_divergence_values[species_name]['all_exo'] = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, all_exo_file)
			list_of_divergence_values[species_name]['all_endo'] = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, all_endo_file)

			################## Get list of 0/1s to store whether a residues in a category are non-subs/subs OR non-gap/gap #############
			# This will be for calculating standard errors using bootstrap
			# Get the list of binary 0/1s representing non-substitution/substitution or non-gap/gap at each residue in a residue category
			exo_div_lists = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, exo_specific_file)
			endo_div_lists = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, endo_specific_file)
			mimicry_div_lists = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, mimicry_file)
			surface_res_div_lists = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, surface_res_file)
			buried_res_div_lists = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, buried_res_file)

			all_exo_div_lists = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, all_exo_file)
			all_endo_div_lists = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, all_endo_file)

			list_of_divergence_values[species_name]['exo_lists'] = exo_div_lists
			list_of_divergence_values[species_name]['endo_lists'] = endo_div_lists
			list_of_divergence_values[species_name]['mimicry_lists'] = mimicry_div_lists
			list_of_divergence_values[species_name]['surface_res_lists'] = surface_res_div_lists
			list_of_divergence_values[species_name]['buried_res_lists'] = buried_res_div_lists

			list_of_divergence_values[species_name]['all_exo_lists'] = all_exo_div_lists
			list_of_divergence_values[species_name]['all_endo_lists'] = all_endo_div_lists


	#print(list_of_divergence_values)
	with open(output_fraction_divergence_json, "w") as outfile:
		json.dump(list_of_divergence_values, outfile)	
	
	# Also output the avg divergence for each protein 
	average_divergence_dict = {}
	subs_or_gaps_list = ['substitutions', 'gaps', 'divergence']
	for subs_or_gaps in subs_or_gaps_list:
		average_divergence_dict[subs_or_gaps] = get_per_protein_avg_divergence(list_of_divergence_values, subs_or_gaps)

	with open(output_average_divergence_json, 'w') as out_avg:
		json.dump(average_divergence_dict, out_avg)

if __name__ == '__main__':
	main()

