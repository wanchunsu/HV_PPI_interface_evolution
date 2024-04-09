'''

Bootstrapping to get standard error of sequence mismatch avgs
=============================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
import os
import json
import numpy as np
import argparse
""" 
For specific residue type r with length n do, for each of the 19 organisms, 
	1. Get list of 0/1s representing no gap/gap or no substitution/substitution at each target protein residue when comapred to organism ortholog
	2. Sample with replacement to get a list of length n 0/1s
	3. Compute % divergence (i.e. [sum(0/1s)/n] * 100%)
Take the avg of the 19 % divergences (1 per organism) => this will be one divergence value

Repeat this 1000 times to get a distribution of divergence values.
Get the standard deviation of this distribution. This will be the standard error for the average divergence value.

"""

##TODO:
'''
Go through list_of_residue_categories, for each res type r
	list_to_store_repeats = []
	repeat 1000 times:
		list_to_store_divergence_vals
		for each org in dict_of_divergence_values:
			input dict_of_divergence_values[org]['subs_list'][r] into fn sample_with_replacement(subs_or_gaps_list)
			store in list_to_store_divergence_vals

		get avg of list_to_store_divergence_vals  and store in list_to_store_repeats

	Get the standard deviation of this distribution. This will be the standard error for the average divergence value 
	that we computed by taking the avg of the original divergence values (stored in dict_of_divergence_values[org]['all_real_mismatches'] and dict_of_divergence_values[org]['all_gaps'])




'''

def sample_with_replacement(subs_or_gaps_list):
	''' Function that samples with replacement len(subs_or_gaps_list) # of times from the subs_or_gaps_list  
	Outputs the sum of the list over the len of the simulated list. (This is essentially a fraction; the # of 1s in the list over the length of the list)
	
	Return:
		simulated_list: simulated list
		simualted_div: sum of the simulated list elements (i.e. # of 1s) over the len of the simulated list
		simualted_percent_div: simualted_div * 100 (percentage format)
	'''
	len_to_sample = len(subs_or_gaps_list)
	simulated_list = np.random.choice(subs_or_gaps_list, size=len_to_sample, replace=True)

	simulated_div = sum(simulated_list)/len(simulated_list)
	simulated_percent_div = simulated_div * 100 #get this number in percentage format

	return simulated_list, simulated_div, simulated_percent_div


def get_bootstrap_standard_err(dict_of_divergence_values, residue_type, gaps_or_subs, repeats = 1000):
	''' For thie given residue type (residue_type) and divergence type (gaps_or_subs), this function
	performs boostrap with `repeats` repeats, each time taking the avg of the simulated divergenes across the 19 organisms
	
	Then outputs the standard deviation of the 1000 repeat results (this will be the standard error for the original observed avg divergence, for plotting purposes)
	
	Return: 
		standard_err_avg_percent_div: standard deviation of the 1000 repeat results in percentage format to accomodate our current divergences in percent format

	'''
	list_to_store_repeats_div = [] #store avg div values for the 1000 repeats
	list_to_store_repeats_percent_div = [] #store avg percent div values for the 1000 repeats

	if gaps_or_subs == 'gaps':
		list_type = 'gaps_list'
	elif gaps_or_subs == 'subs':
		list_type = 'subs_list'
	elif gaps_or_subs == 'all_div':
		list_type = 'div_list'

	for i in range(repeats):
		list_to_store_curr_repeat_div_values = [] #store div values for the 19 organisms in this current repeat
		list_to_store_curr_repeat_percent_div_values = [] #store percent div values for the 19 organisms in this current repeat

		for org in dict_of_divergence_values:

			binary_list = dict_of_divergence_values[org][list_type][residue_type]
			

			simulated_list, simulated_div, simulated_percent_div = sample_with_replacement(binary_list)
			

			list_to_store_curr_repeat_div_values.append(simulated_div)
			list_to_store_curr_repeat_percent_div_values.append(simulated_percent_div)

		# if len(list_to_store_curr_repeat_div_values)!= 19 or len(list_to_store_curr_repeat_percent_div_values)!= 19:
		# 	print("Incorrect number of organisms!")
		
		curr_repeat_avg_div = sum(list_to_store_curr_repeat_div_values) / len(list_to_store_curr_repeat_div_values)
		curr_repeat_avg_percent_div = sum(list_to_store_curr_repeat_percent_div_values)/ len(list_to_store_curr_repeat_percent_div_values)

		list_to_store_repeats_div.append(curr_repeat_avg_div)
		list_to_store_repeats_percent_div.append(curr_repeat_avg_percent_div)

	# print(len(list_to_store_repeats_div), len(list_to_store_repeats_percent_div))  # should be 1000 for both since we have 1000 repeats

	standard_err_avg_div =  np.std(list_to_store_repeats_div, ddof=1) #standard deviation of these 1000 avg divergence values will be the standard error for the original observed avg div value
	standard_err_avg_percent_div = np.std(list_to_store_repeats_percent_div, ddof=1)
	#print(standard_err_avg_div, standard_err_avg_percent_div)

	#we'll just return the percentage-wise standard error since our plots are in percentages.
	return standard_err_avg_percent_div




def main():
	script_dir = osp.dirname(__file__)
	# read in fraction_divergence_with_close_species.json for all, v_target_h, h_target_v
	parser = argparse.ArgumentParser()

	parser.add_argument('-f', '--folder')
	parser.add_argument('-o', '--ortholog_folder')

	args = parser.parse_args()
	folder_name = args.folder #'all', 'h_target_v', 'v_target_h'
	ortholog_folder_label = args.ortholog_folder

	print(f"Getting standard error for avg divergence values in {folder_name}")
	divergence_values = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', ortholog_folder_label, folder_name, 'fraction_divergence_with_close_species.json')
	
	dict_of_divergence_values = {}
	with open(divergence_values) as div_fi:
		dict_of_divergence_values = json.load(div_fi)

	#output file to store standard errors for each residue types gaps and subs
	output_standard_errs_json = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', ortholog_folder_label, folder_name, 'standard_errors.json')
	

	list_of_residue_categories = ['exo', 'endo', 'mimicry', 'surface_res', 'buried_res', 'all_exo', 'all_endo']

	dict_of_bootstrap_standard_errs = {'subs': {}, 'gaps':{}, 'all_div':{}}

	for residue_type in list_of_residue_categories:
		res_type_subs_standard_err = get_bootstrap_standard_err(dict_of_divergence_values, residue_type, 'subs', repeats = 1000)
		res_type_gaps_standard_err =  get_bootstrap_standard_err(dict_of_divergence_values, residue_type, 'gaps', repeats = 1000)
		res_type_div_standard_err =  get_bootstrap_standard_err(dict_of_divergence_values, residue_type, 'all_div', repeats = 1000)

		dict_of_bootstrap_standard_errs['subs'][residue_type] = res_type_subs_standard_err
		dict_of_bootstrap_standard_errs['gaps'][residue_type] = res_type_gaps_standard_err
		dict_of_bootstrap_standard_errs['all_div'][residue_type] = res_type_div_standard_err

	#print(dict_of_bootstrap_standard_errs)

	with open(output_standard_errs_json, 'w') as o:
		json.dump(dict_of_bootstrap_standard_errs, o)

	

if __name__ == '__main__':
	main()
