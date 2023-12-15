'''

Categorize interfacial residues into exogenous-specific, 
endogenous-specific, and mimicked.
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
import os
import argparse

def convert_ppis_to_dict(ppi_file_with_uniprot_int_res):
	dict_of_int_res_on_uniprot = {}

	with open(ppi_file_with_uniprot_int_res) as ppi_file:
		header = next(ppi_file)
		
		for line in ppi_file:
			line = line.strip("\n")
			info = line.split("\t")
			uniprot_id = info[0]
			int_res_on_uniprot = info[5].split(",")
			int_res_on_uniprot = [int(r) for r in int_res_on_uniprot] #for sorting purposes (easier to look at sorted residues)

			if uniprot_id not in dict_of_int_res_on_uniprot:
				dict_of_int_res_on_uniprot[uniprot_id] = []

			dict_of_int_res_on_uniprot[uniprot_id].extend(int_res_on_uniprot)
			sorted_and_unique = sorted(list(set(dict_of_int_res_on_uniprot[uniprot_id])))
			dict_of_int_res_on_uniprot[uniprot_id] = sorted_and_unique 

	return dict_of_int_res_on_uniprot

def group_int_res_by_uniprot_id(dict_of_int_res_on_uniprot, uniprot_int_res_grouped_by_uniprot_id):
	""" Make file with interfacial residues (based on uniprot seq residues) grouped by uniprot id
	Args:
		dict_of_int_res_on_uniprot: dict of uniprot ids and their interfacial residues (based on uniprot seq)
		uniprot_int_res_grouped_by_uniprot_id: output file with interfacial residues grouped by uniprot id

	"""
	
	with open( uniprot_int_res_grouped_by_uniprot_id, 'w') as o:
		if "exogenous" in uniprot_int_res_grouped_by_uniprot_id:
			o.write("target_prot" +"\t" + "all_exo_interfacial_residues_on_uniprot_seq" + "\n")
		else:
			o.write("target_prot" +"\t" + "all_endo_interfacial_residues_on_uniprot_seq" + "\n")

		for uniprotid in dict_of_int_res_on_uniprot:
			int_res_converted_to_str = [str(i) for i in dict_of_int_res_on_uniprot[uniprotid]] #convert to string to print

			o.write(uniprotid +"\t" + ",".join(int_res_converted_to_str) + "\n")



def categorize_exo_and_endo_only_and_mimicry(endo_dict_of_int_res_on_uniprot, exo_dict_of_int_res_on_uniprot, exo_specific_file, endo_specific_file, mimicry_file):
	with open(exo_specific_file, 'w') as exo_only_out, open(endo_specific_file, 'w') as endo_only_out, open(mimicry_file, 'w') as mimicry_out:
		#write headers
		exo_only_out.write("target_prot" +"\t" + "exogenous_only_interfacial_residues"+"\n")
		endo_only_out.write("target_prot" +"\t" + "endogenous_only_interfacial_residues"+"\n")
		mimicry_out.write("target_prot" +"\t" + "mimicked_interfacial_residues"+"\n")

		for uniprot_id in exo_dict_of_int_res_on_uniprot: 

			exo_int_res = exo_dict_of_int_res_on_uniprot[uniprot_id]
			exo_int_res = [str(i) for i in exo_int_res] #convert to string for printing

			if uniprot_id not in endo_dict_of_int_res_on_uniprot: #this uniprot seq only has exogenous interfacial residues (no endogenous)
				exo_only_out.write(uniprot_id + "\t" + ",".join(exo_int_res) +"\n")
				continue

			endo_int_res = endo_dict_of_int_res_on_uniprot[uniprot_id]
			endo_int_res = [str(i) for i in endo_int_res] #convert to string for printing

			exo_specific = [r for r in exo_int_res if r not in endo_int_res] #save only exogenous-specific int res
			endo_specific = [r for r in endo_int_res if r not in exo_int_res] #save only endogenous-specific int res
			mimicry = [r for r in endo_int_res if r in exo_int_res] #save interfacial residues that are both exogenous and endogenous (i.e. mimicked residues)

			
			if exo_specific != []:
				exo_only_out.write(uniprot_id + "\t" +  ",".join(exo_specific) +"\n")

			if endo_specific != []:
				endo_only_out.write(uniprot_id + "\t" + ",".join(endo_specific) +"\n")

			if mimicry != []:
				mimicry_out.write(uniprot_id + "\t" + ",".join(mimicry) + "\n")
			



#### Printing stats ####


def count_residues(int_res_grouped_by_uniprot_seq):
	total_num_int_res = 0
	with open(int_res_grouped_by_uniprot_seq) as int_res:
		next(int_res) 
		for line in int_res:
			line = line.strip("\n")
			info = line.split("\t")
			if info[1] == '': # this targ prot has not residues in this category
				continue
			int_res = info[1].split(",")
			num_int_res = len(int_res)
			total_num_int_res +=num_int_res
	return total_num_int_res

def output_list_of_all_mapped_human_targ_prots(exo_uniprot_int_res_grouped_by_uniprot_id, list_of_hum_targs_fi):
	list_of_human_targs = []
	with open(exo_uniprot_int_res_grouped_by_uniprot_id) as exo_fi:
		next(exo_fi)
		for line in exo_fi:
			list_of_human_targs.append(line.split("\t")[0])

	with open(list_of_hum_targs_fi, 'w') as o:
		for human_targ in list_of_human_targs:
			o.write(human_targ + "\n")


def main():

	script_dir = osp.dirname(__file__)
	categorized_exo_endo_folder = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues")
	

	# Get command line arguments for specific folder and exogenous file
	# We have three types of exogenous files
	#1. all exogenous ppis (this is essentially 2. + 3.)
	#2. exogenous ppis where the virus protein targets the human protein
	#3. exogenous ppis where the human protein targets the virus protein 
	parser = argparse.ArgumentParser()

	parser.add_argument('-x', '--exo_file_uniprot_int_res')
	parser.add_argument('-n', '--endo_file_uniprot_int_res')
	parser.add_argument('-f', '--folder')
	
	
	args = parser.parse_args()
	exo_results = args.exo_file_uniprot_int_res #'exogenous_ppis_with_int_res_on_uniprot_seq.tsv', 'h_targeting_v_exogenous_ppis_with_int_res_on_uniprot_seq.tsv', or 'v_targeting_h_exogenous_ppis_with_int_res_on_uniprot_seq.tsv'
	endo_results = args.endo_file_uniprot_int_res #'endogenous_ppis_with_int_res_on_uniprot_seq.tsv', 'h_targeting_v_endogenous_ppis_with_int_res_on_uniprot_seq.tsv', or 'v_targeting_h_endogenous_ppis_with_int_res_on_uniprot_seq.tsv'
	folder_name = args.folder #'all', 'human_target_virus', 'virus_target_human'
	
	exo_ppi_file_with_uniprot_int_res = osp.join(script_dir, '..', 'data', 'exo_and_endo', exo_results)
	endo_ppi_file_with_uniprot_int_res = osp.join(script_dir, '..', 'data', 'exo_and_endo', endo_results)
	folder_path = osp.join(categorized_exo_endo_folder, folder_name)

	endo_uniprot_int_res_grouped_by_uniprot_id = osp.join(folder_path, "endogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")
	exo_uniprot_int_res_grouped_by_uniprot_id = osp.join(folder_path, "exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")

	exo_specific_file =  osp.join(folder_path, "exogenous_specific_interfacial_residues.tsv")
	endo_specific_file = osp.join(folder_path, "endogenous_specific_interfacial_residues.tsv")
	mimicry_file = osp.join(folder_path, "mimicked_interfacial_residues.tsv")

	if folder_name == 'all': #save list of all human target proteins that map to uniprot seqs 
		list_of_hum_targs_fi = osp.join(categorized_exo_endo_folder, 'list_of_human_targs.txt')
		output_list_of_all_mapped_human_targ_prots(exo_uniprot_int_res_grouped_by_uniprot_id, list_of_hum_targs_fi)
		
	if not osp.exists(folder_path):
		os.makedirs(folder_path)

		
	endo_dict_of_int_res_on_uniprot = convert_ppis_to_dict(endo_ppi_file_with_uniprot_int_res)
	exo_dict_of_int_res_on_uniprot = convert_ppis_to_dict(exo_ppi_file_with_uniprot_int_res)

	#Make file with interfacial residues (based on uniprot seq residues) grouped by uniprot id 
	group_int_res_by_uniprot_id(endo_dict_of_int_res_on_uniprot, endo_uniprot_int_res_grouped_by_uniprot_id)
	group_int_res_by_uniprot_id(exo_dict_of_int_res_on_uniprot, exo_uniprot_int_res_grouped_by_uniprot_id)

	#Categorize endo-specific, exo-specific, and mimicry interfacial residues:
	categorize_exo_and_endo_only_and_mimicry(endo_dict_of_int_res_on_uniprot, exo_dict_of_int_res_on_uniprot, exo_specific_file, endo_specific_file, mimicry_file)


	#Printing stats:
	print(f'\n\n### Printing stats for {folder_name}: ###')
	
	num_endo_specific_int_res = count_residues(endo_specific_file)
	num_exo_specific_int_res = count_residues(exo_specific_file)
	num_mimicked_int_res = count_residues(mimicry_file)
	num_all_endo_int_res = num_endo_specific_int_res + num_mimicked_int_res #essentially same as count_residues(endo_uniprot_int_res_grouped_by_uniprot_id)
	num_all_exo_int_res = num_exo_specific_int_res + num_mimicked_int_res #essentially same as count_residues(exo_uniprot_int_res_grouped_by_uniprot_id)

	print(f'Total number of endo interfacial residues on human target proteins (uniprot seq): {num_all_endo_int_res}')
	print(f'Total number of exo interfacial residues on human target proteins (uniprot seq): {num_all_exo_int_res}')

	print("___________________________________________________________________________________\n")


	print(f'Total number of endogenous-specific interfacial residues on human target proteins (uniprot seq): {num_endo_specific_int_res}')
	print(f'Total number of exogenous-specific interfacial residues on human target proteins (uniprot seq): {num_exo_specific_int_res}')
	print(f'Total number of mimicked interfacial residues on human target proteins (uniprot seq): {num_mimicked_int_res}')

	print("___________________________________________________________________________________\n")
	num_all_int_res_no_repeats = num_endo_specific_int_res + num_exo_specific_int_res + num_mimicked_int_res
	print(f'Total number of non-repeating interfacial residues: {num_all_int_res_no_repeats}')
	

	
if __name__ == '__main__':
	main()

