'''

Normalize raw Rate4Site scores and compute overall evolutionary scores 
======================================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp
import argparse
from Bio import SeqIO
import numpy as np
import json

'''
Two different score normalization techniques:
1) Z-score normalization (mean = 1, std = 0) Default rate4site output. Will be using parse_rate4site_scores_into_dict with file_ext = '.scores'. This will parse z-score-normalized scores into dict.
2) [THIS SCRIPT] Get raw scores from rate4site and normalize by mean on that protein. Will be using compute_mean_normalized_rate4site_scores with file_ext = '.rawscores'. This will normalize by mean and parse normalized scores into a dict
'''

def compute_mean_normalized_rate4site_scores(unnormalized_rate4site_results_folder, file_ext):
	""" Normalize raw rate4site scores by mean i.e. calculate normalized scores by dividing each site's unnormalized score by the mean score on the whole protein sequence 
	Parse conservation scores into dictionary with targ_prot id as outer key, residue number as inner key, and rate4site score for that residue as value
	
	Args:
		unnormalized_rate4site_results_folder: folder storing unnormalized rate4site scores (same folder as rate4site pre-normalized)
		file_ext: file extension specifiying the type of scores we want to use (.rawscores)
	Return:
		dict_of_rate4site_scores: dictionary storing scores normalized by mean 

	"""
	
	dict_of_rate4site_scores ={}
	for targ_prot_rate4site_results in os.listdir(unnormalized_rate4site_results_folder):
		if targ_prot_rate4site_results.endswith(file_ext):
			targ_prot = targ_prot_rate4site_results.split('.')[0]
			dict_of_rate4site_scores[targ_prot] = {}
			mean_score_for_prot = 0
			with open(osp.join(unnormalized_rate4site_results_folder, targ_prot_rate4site_results)) as res:
				#print(targ_prot_rate4site_results)
				for line in res:
					if line.startswith("#") and 'Average' in line:
						line = line.strip("\n")
						info = line.split()
						mean_score_for_prot = float(info[2])
			with open(osp.join(unnormalized_rate4site_results_folder, targ_prot_rate4site_results)) as res:
				unnormalized_scores = []
				for line in res:
					#print(line)
					if not line.startswith("#") and not line == "\n" : #these are the headers and empty lines respectively
						line = line.strip("\n")
						info = line.split()
						#print(info)
						res_pos = int(info[0])
						res_score = float(info[2])
						mean_normalized_res_score = res_score/mean_score_for_prot
						# print(unnormalized_rate4site_results_folder, targ_prot, res_pos, mean_normalized_res_score) #for verfication purposes -- verified
						dict_of_rate4site_scores[targ_prot][res_pos] = mean_normalized_res_score
						

	#print(len(dict_of_rate4site_scores))
	return dict_of_rate4site_scores

def calculate_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, residue_results):
	scores_for_residue_category = []

	with open(residue_results) as res:
		next(res)
		for line in res:
			line = line.strip("\n")
			targ_prot = line.split("\t")[0]
			if targ_prot not in dict_of_rate4site_scores:
				continue
			if line.split("\t")[1] == '':
				continue
			residues = line.split("\t")[1].split(",") #split up residues separated by "," into a list of residues
			residues = [int(r) for r in residues] #convert to int
			score_for_this_prot = [dict_of_rate4site_scores[targ_prot][r] for r in residues]
			 
			
			scores_for_residue_category.extend(score_for_this_prot)
			
	total_num_residues_in_category = len(scores_for_residue_category)
	mean_score_for_category = np.mean(scores_for_residue_category)
	std_for_category = np.std(scores_for_residue_category, ddof=1)
	stderr_for_category = std_for_category/np.sqrt(total_num_residues_in_category)

	#print(scores_for_residue_category)
	return mean_score_for_category, std_for_category, stderr_for_category, scores_for_residue_category


def calculate_per_protein_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, residue_results):
	scores_for_residue_category = {}

	with open(residue_results) as res:
		next(res)
		for line in res:
			line = line.strip("\n")
			targ_prot = line.split("\t")[0]
			if targ_prot not in dict_of_rate4site_scores:
				continue
			if line.split("\t")[1] == '':
				continue
			residues = line.split("\t")[1].split(",") #split up residues separated by "," into a list of residues
			residues = [int(r) for r in residues] #convert to int
			score_for_this_prot = [dict_of_rate4site_scores[targ_prot][r] for r in residues]
			
			total_num_residues_in_category = len(score_for_this_prot)
			mean_score_for_category = np.mean(score_for_this_prot)
			if len(score_for_this_prot) == 1: #only on value so std and stderr are both 0
				std_for_category = 0
				stderr_for_category = 0
			else:
				std_for_category = np.std(score_for_this_prot, ddof=1)
				stderr_for_category = std_for_category/np.sqrt(total_num_residues_in_category)


			scores_for_residue_category[targ_prot] = {"mean": mean_score_for_category, "std": std_for_category, "stderr": stderr_for_category, "scores": score_for_this_prot}
			
	

	#print(scores_for_residue_category)
	return scores_for_residue_category


def main():
	script_dir = osp.dirname(__file__)
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--folder')
	parser.add_argument('-r', '--rate4site_result_type')
	parser.add_argument('-o', '--overall_scores_folder')
	parser.add_argument('-a', '--all_residues_scores_folder')
	parser.add_argument('-p', '--per_protein_scores_folder')
	parser.add_argument('-l', '--ortholog_folder')
	args = parser.parse_args()

	folder_name = args.folder #'all', 'h_target_v', 'v_target_h'
	rate4site_results = args.rate4site_result_type # 'rate4site_autogenerate_tree' folder or 'rate4site_species_tree' folder	
	ortholog_folder_label = args.ortholog_folder
	rate4site_results_folder = osp.join('..', 'data', 'homologs_and_conservation', ortholog_folder_label, folder_name, rate4site_results)

	output_overall_scores_folder_name = args.overall_scores_folder #'rate4site_overall_scores' (for all orthologs), 'rate4site_overall_scores_human_mouse', rate4site_overall_scores_chimpanzee_cow_dog_gorilla_macaque_mouse_orangutan_pig_rabbit_rat (w/out chicken)
	output_all_residue_scores_folder_name = args.all_residues_scores_folder #'rate4site_all_residue_scores' (for all orthologs), 'rate4site_all_residue_scores_human_mouse', rate4site_all_residue_scores_chimpanzee_cow_dog_gorilla_macaque_mouse_orangutan_pig_rabbit_rat (w/out chicken)
	output_per_protein_scores_folder_name = args.per_protein_scores_folder #'rate4site_per_protein_scores'

	# Files with exo-specific, endo-specific, and mimicry interfacial residues as well as surface and buried residues
	exo_specific_file =  osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_specific_interfacial_residues.tsv")
	endo_specific_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_specific_interfacial_residues.tsv")
	mimicry_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "mimicked_interfacial_residues.tsv")
	surface_res_file = osp.join(script_dir, '..', 'data', "surface_residues", folder_name, "valid_surface_residues_on_uniprot_seq.tsv")
	buried_res_file = osp.join(script_dir, '..', 'data', "buried_residues", folder_name, "valid_buried_residues_on_uniprot_seq.tsv")

	all_exo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")
	all_endo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")

	#output folder for the overall ragte4site scores across the residue categories
	output_folder = osp.join('..', 'data', 'homologs_and_conservation', ortholog_folder_label, output_overall_scores_folder_name)
	if not osp.exists(output_folder):
		os.mkdir(output_folder)

	rate4site_overall_scores_for_residue_categories_mean_normalized = osp.join(output_folder, folder_name + "_" + rate4site_results + ".mean_normalized.json")
	
	#output folder to store all the scores for all the residues in a given category
	rate4site_all_residue_category_scores_folder = osp.join('..', 'data', 'homologs_and_conservation', ortholog_folder_label, output_all_residue_scores_folder_name)
	if not osp.exists(rate4site_all_residue_category_scores_folder):
		os.mkdir(rate4site_all_residue_category_scores_folder)
	
	rate4site_all_scores_for_residue_categories_mean_normalized = osp.join(rate4site_all_residue_category_scores_folder, folder_name + "_" + rate4site_results + "_all_residue_scores.mean_normalized.json")


	#output folder to store all the scores on a per protein basis in each residue category
	rate4site_per_protein_folder = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', ortholog_folder_label, output_per_protein_scores_folder_name)
	if not osp.exists(rate4site_per_protein_folder):
		os.mkdir(rate4site_per_protein_folder)

	rate4site_per_protein_scores_mean_normalized = osp.join(rate4site_per_protein_folder, folder_name + '_' + rate4site_results + '_per_protein_scores.mean_normalized.json')


	############ Using mean-normalized scores --we use raw rate4site scores and normalize each raw score by the mean score on that protein ############
	print('############### Computing and using mean-normalized rate4site scores ###############') 
	print(f"Getting scores for {folder_name}, {rate4site_results} . . . \n")

	#normalize and format scores for all residues on the target proteins into a dictionary
	dict_of_rate4site_scores = compute_mean_normalized_rate4site_scores(rate4site_results_folder, '.rawscores') #we want to use the raw (unormalized) rate4site scores. This function will normalize the scores by the mean score on that protein and then calculate the overall score
	
	scores = {'exo':{}, 'endo':{}, 'mimicry':{}, 'surface_res':{}, 'buried_res':{}, 'all_exo': {}, 'all_endo': {}} #dict to store the mean, std and standard error for the residues in the specified category
	all_residue_scores = {'exo':[], 'endo':[], 'mimicry':[], 'surface_res':[], 'buried_res':[], 'all_exo': [], 'all_endo': []}
	per_protein_scores = {'exo':{}, 'endo':{}, 'mimicry':{}, 'surface_res':{}, 'buried_res':{}, 'all_exo': {}, 'all_endo': {}}

	# Get the scores for the residues in a given category and calculate the mean, std, and standard error
	exo_mean, exo_std, exo_stderr, exo_scores = calculate_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, exo_specific_file)
	endo_mean, endo_std, endo_stderr, endo_scores = calculate_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, endo_specific_file)
	mimicry_mean, mimicry_std, mimicry_stderr, mimicry_scores = calculate_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, mimicry_file)
	surface_res_mean, surface_res_std, surface_res_stderr, surface_res_scores= calculate_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, surface_res_file)
	buried_res_mean, buried_res_std, buried_res_stderr, buried_res_scores = calculate_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, buried_res_file)

	all_exo_mean, all_exo_std, all_exo_stderr, all_exo_scores = calculate_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, all_exo_file)
	all_endo_mean, all_endo_std, all_endo_stderr, all_endo_scores = calculate_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, all_endo_file)

	#Store mean, std and standard error under corresponding key
	scores['exo']['mean'], scores['exo']['std'], scores['exo']['stderr'] = exo_mean, exo_std, exo_stderr
	scores['endo']['mean'], scores['endo']['std'],scores['endo']['stderr'] = endo_mean, endo_std, endo_stderr
	scores['mimicry']['mean'], scores['mimicry']['std'], scores['mimicry']['stderr'] = mimicry_mean, mimicry_std, mimicry_stderr
	scores['surface_res']['mean'], scores['surface_res']['std'],scores['surface_res']['stderr'] = surface_res_mean, surface_res_std, surface_res_stderr
	scores['buried_res']['mean'], scores['buried_res']['std'],scores['buried_res']['stderr'] = buried_res_mean, buried_res_std, buried_res_stderr

	scores['all_exo']['mean'], scores['all_exo']['std'], scores['all_exo']['stderr'] = all_exo_mean, all_exo_std, all_exo_stderr
	scores['all_endo']['mean'], scores['all_endo']['std'],scores['all_endo']['stderr'] = all_endo_mean, all_endo_std, all_endo_stderr

	#Store all the residue scores for each category
	all_residue_scores['exo'] = exo_scores
	all_residue_scores['endo'] = endo_scores
	all_residue_scores['mimicry'] = mimicry_scores
	all_residue_scores['surface_res'] = surface_res_scores
	all_residue_scores['buried_res'] = buried_res_scores

	all_residue_scores['all_exo'] = all_exo_scores
	all_residue_scores['all_endo'] = all_endo_scores


	# Get and store the per protein scores in each category
	per_protein_scores['exo'] = calculate_per_protein_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, exo_specific_file)
	per_protein_scores['endo'] = calculate_per_protein_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, endo_specific_file)
	per_protein_scores['mimicry'] = calculate_per_protein_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, mimicry_file)
	per_protein_scores['surface_res'] = calculate_per_protein_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, surface_res_file)
	per_protein_scores['buried_res'] = calculate_per_protein_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, buried_res_file)

	per_protein_scores['all_exo'] = calculate_per_protein_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, all_exo_file)
	per_protein_scores['all_endo'] = calculate_per_protein_overall_score_for_the_given_residue_category(dict_of_rate4site_scores, all_endo_file)

	#output the scores dictionary
	with open(rate4site_overall_scores_for_residue_categories_mean_normalized, 'w') as output_fi:
		json.dump(scores, output_fi)

	#output all the residue scores dict
	with open(rate4site_all_scores_for_residue_categories_mean_normalized, 'w') as all_residue_scores_output_fi:
		json.dump(all_residue_scores, all_residue_scores_output_fi)

	#output per protein scores for each residue category
	with open(rate4site_per_protein_scores_mean_normalized, 'w') as per_protein_scores_output_fi:
		json.dump(per_protein_scores, per_protein_scores_output_fi)

if __name__ == '__main__':
	main()