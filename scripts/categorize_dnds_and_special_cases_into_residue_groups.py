'''

Categorize dN/dS values into residue categories and output
overall scores for each category.
==========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import json
import os
import os.path as osp
import argparse
import numpy as np

def parse_all_dnds_on_targ_prots_into_a_dict(dnds_on_targ_prots_folder):
	""" Group all dnds values on targ prots into a single nested dict (with targ prot ids as outer key)

	"""
	all_dnds_on_targ_prots_dict = {}

	for dnds_targ_prot in os.listdir(dnds_on_targ_prots_folder):
		targ_prot_id = dnds_targ_prot.split('.')[0]
		with open(osp.join(dnds_on_targ_prots_folder, dnds_targ_prot)) as dnds_tp:
			targ_prot_dnds_dict = json.load(dnds_tp)
			all_dnds_on_targ_prots_dict[targ_prot_id] = targ_prot_dnds_dict

	return all_dnds_on_targ_prots_dict

def get_scores_for_category(residues_file, all_dnds_on_targ_prots_dict):
	# Store all dn/ds values for a certain residue type (exo, endo, mimicry, surface, or buried) from all targ proteins into a list
	list_of_dnds_vals_across_proteins = []
	special_cases_across_proteins = 0
	res_without_dnds_val = 0
	total_num_res = 0

	# Also store dict with the per protein residues
	dict_of_dnds_vals_on_protein_residues = {}
	with open(residues_file) as res:
		next(res)
		for line in res:
			line = line.strip("\n")
			targ_prot = line.split("\t")[0]

			if targ_prot not in all_dnds_on_targ_prots_dict: #case where protein is an immunoglobulin (so doesn't have cds sequence -- thus ignored for dnds analysis)
				continue
			if line.split("\t")[1] == '':
				continue
			residues = line.split("\t")[1].split(",") #split up residues separated by "," into a list of residues
			if residues == []:
				continue
			total_num_res += len(residues)

			well_defined_dnds_vals_on_prot_residues = [all_dnds_on_targ_prots_dict[targ_prot][pos][0] for pos in all_dnds_on_targ_prots_dict[targ_prot] if pos in residues and all_dnds_on_targ_prots_dict[targ_prot][pos][1]=='normal']

			special_cases = [pos for pos in all_dnds_on_targ_prots_dict[targ_prot] if pos in residues and all_dnds_on_targ_prots_dict[targ_prot][pos][1]!='normal']#dS==0 , either no synoynmous substitutions or all but one residue at that position is a gap

			no_dnds_val = [pos for pos in residues if pos not in all_dnds_on_targ_prots_dict[targ_prot]] #residues in this category that don't have dN/dS values (b/c didn't map to a cds region, so couldn't get its dN/dS)
			#There are 0 residues w/out dnds values, so no_dnds_val should be empty

			list_of_dnds_vals_across_proteins.extend(well_defined_dnds_vals_on_prot_residues) #list of well-defined dnds values across all proteins
			special_cases_across_proteins += len(special_cases) #list of residues that are special cases across all proteins
			res_without_dnds_val += len(no_dnds_val)  #list of residues that don't have dnds values across all proteins(so far none at all -- all protein residues have a dnds value)
			#print(len(dnds_vals_on_prot_residues))

			# Store per protein site-specific dnds (don't store if there are no residues on this protein that belong in the given category)
			if len(well_defined_dnds_vals_on_prot_residues) != 0:
				avg_dnds_on_prot_for_res_category = np.mean(well_defined_dnds_vals_on_prot_residues)
				if len(well_defined_dnds_vals_on_prot_residues) == 1:
					std_on_prot_for_res_category = 0
					stderr_on_prot_for_res_category = 0
				else:
					std_on_prot_for_res_category = np.std(well_defined_dnds_vals_on_prot_residues, ddof=1)
					stderr_on_prot_for_res_category = std_on_prot_for_res_category/np.sqrt(len(well_defined_dnds_vals_on_prot_residues))

				dict_of_dnds_vals_on_protein_residues[targ_prot] = {"dnds": well_defined_dnds_vals_on_prot_residues, "mean": avg_dnds_on_prot_for_res_category, "std": std_on_prot_for_res_category, "stderr": stderr_on_prot_for_res_category, "special": len(special_cases), "no_dnds": len(no_dnds_val), "total_res": len(residues)}
			#print(len(dict_of_dnds_vals_on_protein_residues[targ_prot]))
	# print(total_num_res)
	# print(f'{len(list_of_dnds_vals_across_proteins)}, {special_cases_across_proteins}, {res_without_dnds_val}')
	info_across_proteins = {"dnds": list_of_dnds_vals_across_proteins, "special": special_cases_across_proteins, "no_dnds": res_without_dnds_val, "total_res": total_num_res}

	return info_across_proteins, dict_of_dnds_vals_on_protein_residues
	
def main():
	script_dir = osp.dirname(__file__)
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--hyphy_folder')

	args = parser.parse_args()

	hyphy_folder = args.hyphy_folder

	folder_name = 'v_target_h' #we're only working with virus_target_human residues #args.folder #'all', 'h_target_v', 'v_target_h'


	if 'srv' in hyphy_folder:
		dnds_on_targ_prots_folder = osp.join(script_dir, '..', 'data', 'dnds', 'dnds_on_targ_prots_srv')
		output_folder = osp.join(script_dir, '..', 'data', 'dnds', 'dnds_for_residue_categories_srv')
		if not osp.exists(output_folder):
			os.mkdir(output_folder)
	elif 'species_tree' in hyphy_folder:
		dnds_on_targ_prots_folder = osp.join(script_dir, '..', 'data', 'dnds', 'dnds_on_targ_prots_species_tree')
		output_folder = osp.join(script_dir, '..', 'data', 'dnds', 'dnds_for_residue_categories_species_tree')
		if not osp.exists(output_folder):
			os.mkdir(output_folder)
	else:
		dnds_on_targ_prots_folder = osp.join(script_dir, '..', 'data', 'dnds', 'dnds_on_targ_prots')
		output_folder = osp.join(script_dir, '..', 'data', 'dnds', 'dnds_for_residue_categories')
		if not osp.exists(output_folder):
			os.mkdir(output_folder)

	outfile_residue_types = osp.join(output_folder, folder_name + '.dnds_for_residues.json')
	outfile_split_by_prots = osp.join(output_folder, folder_name + '.dnds_for_residues_split_by_targ_prot.json')

	# Files with exo-specific, endo-specific, and mimicry interfacial residues as well as surface and buried residues
	exo_specific_file =  osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_specific_interfacial_residues.tsv")
	endo_specific_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_specific_interfacial_residues.tsv")
	mimicry_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "mimicked_interfacial_residues.tsv")
	surface_res_file = osp.join(script_dir, '..', 'data', "surface_residues", folder_name, "valid_surface_residues_on_uniprot_seq.tsv")
	buried_res_file = osp.join(script_dir, '..', 'data', "buried_residues", folder_name, "valid_buried_residues_on_uniprot_seq.tsv")

	all_exo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")
	all_endo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")

	all_dnds_on_targ_prots_dict = parse_all_dnds_on_targ_prots_into_a_dict(dnds_on_targ_prots_folder)

	
	dnds_to_output_res_types = {} #dictionary to store dn/ds values split by residue type (i.e. {'exo':{"dnds": well_defined_dnds_vals_on_prot_residues, "special": num_special_cases, "no_dnds": num_no_dnds_val, "total_res": num_residues}, 'endo': {..}, ...})
	dnds_split_by_prots = {} #dictionary to store dn/ds values split by protein then by residue type (i.e. {'exo': {'targprot1': {"dnds": well_defined_dnds_vals_on_prot_residues, "special": num_special_cases, "no_dnds": num_no_dnds_val, "total_res": num_residues}, 'targprot2': {}, ...}, 'endo': {...},...})

	dnds_to_output_res_types['exo'], dnds_split_by_prots['exo'] = get_scores_for_category(exo_specific_file, all_dnds_on_targ_prots_dict)
	dnds_to_output_res_types['endo'], dnds_split_by_prots['endo'] = get_scores_for_category(endo_specific_file, all_dnds_on_targ_prots_dict)
	dnds_to_output_res_types['mimicry'], dnds_split_by_prots['mimicry'] = get_scores_for_category(mimicry_file, all_dnds_on_targ_prots_dict)
	dnds_to_output_res_types['surface'], dnds_split_by_prots['surface'] = get_scores_for_category(surface_res_file, all_dnds_on_targ_prots_dict)
	dnds_to_output_res_types['buried'], dnds_split_by_prots['buried'] = get_scores_for_category(buried_res_file, all_dnds_on_targ_prots_dict)

	dnds_to_output_res_types['all_exo'], dnds_split_by_prots['all_exo'] = get_scores_for_category(all_exo_file, all_dnds_on_targ_prots_dict)
	dnds_to_output_res_types['all_endo'], dnds_split_by_prots['all_endo'] = get_scores_for_category(all_endo_file, all_dnds_on_targ_prots_dict)

	with open(outfile_residue_types, 'w') as out_res_types, open(outfile_split_by_prots, 'w') as out_split_by_prots:
		json.dump(dnds_to_output_res_types, out_res_types)
		json.dump(dnds_split_by_prots, out_split_by_prots)



if __name__ == '__main__':
	main()