'''

Output lists of potential pathogenic and immune exogenous 
PPIs for verification
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
import os

def get_list_of_potential_ppis_with_selected_annos(potential_ppis_fi):
	list_of_potential_ppis_patho = []
	list_of_potential_ppis_immune = []
	with open(potential_ppis_fi) as p:

		for line in p:
			line = line.strip("\n")
			info = line.split("\t")
			anno = info[-1]
			ppi_to_store = info[0] + "," + info[1]
			if anno == 'F' and ppi_to_store not in list_of_potential_ppis_patho:
				list_of_potential_ppis_patho.append(ppi_to_store)
			if anno == 'T' and ppi_to_store not in list_of_potential_ppis_immune:
				list_of_potential_ppis_immune.append(ppi_to_store)

	return list_of_potential_ppis_patho, list_of_potential_ppis_immune

def get_list_of_known_ppis(known_ppis_fi):
	list_of_known_ppis = []
	with open(known_ppis_fi) as k:
		next(k) #skip header
		for line in k:
			line = line.strip("\n")
			info = line.split("\t")
			ppi_to_store = info[0] + "," + info[1]
			if ppi_to_store not in list_of_known_ppis:
				list_of_known_ppis.append(ppi_to_store)
	return list_of_known_ppis


def output_potential_exo_ppis_that_are_not_found_in_known(potential_ppis_fi, known_ppis_fi, to_check_patho_outfi, to_check_immune_outfi):

	list_of_potential_ppis_patho, list_of_potential_ppis_immune = get_list_of_potential_ppis_with_selected_annos(potential_ppis_fi)
	list_of_known_ppis = get_list_of_known_ppis(known_ppis_fi)
	to_check_patho = []
	to_check_immune = []

	for ppi_pair in list_of_potential_ppis_patho:
		if ppi_pair not in list_of_known_ppis:
			if ppi_pair not in to_check_patho:
				to_check_patho.append(ppi_pair)
	with open(to_check_patho_outfi, 'w') as to_check_fi_patho:
		for pp in to_check_patho:
			to_check_fi_patho.write(pp.split(",")[0] + "\t" + pp.split(",")[1] + "\t\n")
	print(f"Need to verify {len(to_check_patho)} pathogenic exogenous PPIs")	

	# We only want to look at pathogenic PPIs (so no need to do anything with immune -- can just drop them for further analyses)
	# for ppi_pair in list_of_potential_ppis_immune:
	# 	if ppi_pair not in list_of_known_ppis:
	# 		if ppi_pair not in to_check_immune:
	# 			to_check_immune.append(ppi_pair)
	# with open(to_check_immune_outfi, 'w') as to_check_fi_immune:
	# 	for pp in to_check_immune:
	# 		to_check_fi_immune.write(pp.split(",")[0] + "\t" + pp.split(",")[1] + "\t\n")
	# print(f"Need to verify {len(to_check_immune)} immune exogenous PPIs")	


def main():
	script_dir = osp.dirname(__file__)
	exogenous_potential_ppis_fi = osp.join('..', 'data', 'annotated_combined_homology_and_pdb_anno_human_virus_interface_residues.tsv')
	exogenous_known_ppis_fi = osp.join('..', 'data', 'IntAct_human_virus', 'intact_exogenous_human_ppis_with_valid_ids.tsv')
	
	#Files to store list of ppis to verify manually
	exogenous_to_check_patho_outfi = osp.join('..', 'data', 'pdb_annotations', 'verify_ppis', 'pathogenic_exogenous_ppis_to_verify.tsv')
	exogenous_to_check_immune_outfi = osp.join('..', 'data', 'pdb_annotations', 'verify_ppis', 'immune_exogenous_ppis_to_verify.tsv')

	if not osp.exists(osp.dirname(exogenous_to_check_patho_outfi)):
		os.mkdir(osp.dirname(exogenous_to_check_patho_outfi))


	output_potential_exo_ppis_that_are_not_found_in_known(exogenous_potential_ppis_fi, exogenous_known_ppis_fi, exogenous_to_check_patho_outfi, exogenous_to_check_immune_outfi)
	
	# Update: No need to do this for endogenous b/c we don't use pdb annotations for endo PPIs (already have enough from known PPIs)
	#output_potential_ppis_that_are_not_found_in_known(endogenous_potential_ppis_fi, endogenous_known_ppis_fi, endogenous_to_check_outfi, 'endo')

if __name__ == '__main__':
	main()