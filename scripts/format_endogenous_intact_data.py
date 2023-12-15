'''

Process within-human (endogenous) PPI data from IntAct 
======================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
from uniprot_tools import get_uniprot_sequences, write_uniprot_seqs_to_fasta, obsolete_or_alt_ids, get_unreviewed
import os


def format_intact(input_fi, output_fi):
	""" Format and process intact data so that we only have human-human ppis that have physical and direct associations (keep the ones with physical association, association, and direct interaction)
	Look into 'association'

	Args:
		input_fi: original downloaded intact data
		output_fi; processed and formatted intact data
	"""
	with open(input_fi) as infi:
		next(infi) #skip header line
		with open(output_fi, 'w') as outfi:
			outfi.write("interactorA\tinteractorB\ttaxidA\ttaxidB\tdetection_method\tinteraction_type\tpublication_id\n")
			for line in infi:
				#info = line.split("\t")

				
				line = line.strip("\n")
				info = line.split("\t")

				interactorA = info[0].replace("uniprotkb:", "")
				interactorB = info[1].replace("uniprotkb:", "")

				detection_method = info[6]
				publication_id  = info[8]
				taxidA = info[9]
				taxidB = info[10]
				interaction_type = info[11].split("(")[1].replace(")", "")

				interactor_typeA = info[20]
				interactor_typeB = info[21]
				if interactor_typeA == '-' or interactor_typeB == '-':
					continue
				interactor_typeA = info[20].split("(")[1].replace(")", "")
				interactor_typeB = info[21].split("(")[1].replace(")", "")

				#check that both binding partners are human proteins and that they are physically interacting (we will take 'direct interaction', 'physical association', and 'association')
				#still need to check if 'association' is a valid interaction type to be considered
				if 'taxid:9606' in taxidA and 'taxid:9606' in taxidB: 
					
					if interaction_type == "direct interaction" or interaction_type == "physical association" or interaction_type == 'association': #double check 'assocaition' by looking at publications
						#get just the taxid
						taxidA = taxidA.split("|")[0].split("(")[0]
						taxidB = taxidB.split("|")[0].split("(")[0]
						
						if interactor_typeA == 'protein' and interactor_typeB == 'protein': #must be protein-protein interaction
							# print(interactorA + '\t' + interactorB + '\t' + taxidA + '\t' + taxidB + '\t' + detection_method + '\t' + interaction_type +  '\t' + interactor_typeA + '\t' + interactor_typeB +'\t' + publication_id +"\n")
							
							outfi.write(interactorA + '\t' + interactorB + '\t' + taxidA + '\t' + taxidB + '\t' + detection_method + '\t' + interaction_type + '\t'+  publication_id +"\n")

def get_exogenous_target_proteins(exogenous_ppis, list_of_human_target_proteins):
	""" Get list of unique human proteins that are involved in exogenous human-virus ppis
	Args:
		exogenous_ppis: file with exogenous ppis (i.e. human-virus protein pairs with interfacial residues)
		list_of_human_proteins: output file with the list of human target proteins that are involved in exogenous human-virus ppis (this list will be used for getting endogenous ppis)
	Return:
		list of human target proteins
	"""
	list_of_targ_prots = []
	with open(exogenous_ppis) as infi:
		next(infi)
		for line in infi:
			human_prot = line.split("\t")[0]
			if human_prot not in list_of_targ_prots:
				list_of_targ_prots.append(human_prot)
	with open(list_of_human_target_proteins, 'w') as out_list_fi:				
		out_list_fi.write('\n'.join(list_of_targ_prots))
	return list_of_targ_prots

# def get_only_human_ppis_with_at_least_2_pubs(all_human_ppis_fi):
# 	list_of_2_hit_ppis = []
# 	dict_of_ppi_hits = {}
# 	with open(all_human_ppis_fi, 'r') as all_human_ppis:
# 		for line in all_human_ppis:
# 			line = line.strip("\n")
# 			info = line.split("\t")
# 			interactorA = info[0].split("-")[0] #take canonical protein rather than isoform or segment (if applicable) (i.e. delete everything including and after the "-")
# 			interactorB = info[1].split("-")[0] #take canonical protein rather than isoform or segment (if applicable) (i.e. delete everything including and after the "-")
# 			pubmedid = info[6]
# 			if (interactorA, interactorB) not in dict_of_ppi_hits and (interactorB, interactorA) not in dict_of_ppi_hits:
# 				dict_of_ppi_hits[(interactorA, interactorB)] = [pubmedid]
# 			elif (interactorA, interactorB) in dict_of_ppi_hits:
# 				dict_of_ppi_hits[(interactorA, interactorB)].append(pubmedid)
# 			elif (interactorB, interactorA) in dict_of_ppi_hits:
# 				dict_of_ppi_hits[(interactorB, interactorA)].append(pubmedid)
# 	for ppi in dict_of_ppi_hits:
# 		pubmedids_unique = list(set(dict_of_ppi_hits[ppi]))
# 		if len(pubmedids_unique) > 1:
# 			list_of_2_hit_ppis.append(ppi)
# 	return list_of_2_hit_ppis


def get_endogenous_ppis(exogenous_ppis, list_of_human_target_proteins, all_human_ppis_fi, endogenous_ppis_with_target_proteins): #(exogenous_ppis, list_of_human_target_proteins, list_of_2_hit_ppis, all_human_ppis_fi, endogenous_ppis_with_target_proteins):
	list_of_human_partner_pairs = [] #order matters(we want to store the target protein as the first protein-- interactorA), 
	#E.g. if protA and protB are both target proteins and interact, we will store protA-protB (in terms of protA) and protB-protA (in terms of protB) as different entries

	list_of_targ_prots = get_exogenous_target_proteins(exogenous_ppis, list_of_human_target_proteins)

	list_of_targ_prots_with_intact_endo = [] #keep track of number of target proteins that have endogenous interactions
	with open(all_human_ppis_fi, 'r') as all_human_ppis:
		with open(endogenous_ppis_with_target_proteins, 'w') as outfi:
			for line in all_human_ppis:
				line = line.strip("\n")
				info = line.split("\t")
				interactorA = info[0].split("-")[0] #take canonical protein rather than isoform or segment (if applicable) (i.e. delete everything including and after the "-")
				interactorB = info[1].split("-")[0] #take canonical protein rather than isoform or segment (if applicable) (i.e. delete everything including and after the "-")
				
				if interactorA == interactorB: #not including self-interactions
					continue

				# #Only take PPIs that have at least 2 publication hits
				# if (interactorA, interactorB) not in list_of_2_hit_ppis and (interactorB, interactorA) not in list_of_2_hit_ppis:
				# 	continue

				if len(interactorA) != 6 or len(interactorB) != 6: #not a uniprot id
					continue
				if interactorA in list_of_targ_prots:
					if (interactorA, interactorB) not in list_of_human_partner_pairs:
						list_of_human_partner_pairs.append((interactorA, interactorB))
						outfi.write(interactorA + "\t" + interactorB + "\n")

						#keep track of number of target proteins that have endogenous interactions
						# if interactorA not in list_of_targ_prots_with_intact_endo:
						# 	list_of_targ_prots_with_intact_endo.append(interactorA)
						

				if interactorB in list_of_targ_prots:
					if (interactorB, interactorA) not in list_of_human_partner_pairs:
						list_of_human_partner_pairs.append((interactorB, interactorA))
						outfi.write(interactorB + "\t" + interactorA + "\n")

						#keep track of number of target proteins that have endogenous interactions
						# if interactorB not in list_of_targ_prots_with_intact_endo:
						# 	list_of_targ_prots_with_intact_endo.append(interactorB)

	#print("\n".join(list_of_targ_prots_with_intact_endo))

	
	


def get_list_of_human_proteins_from_ppis(endogenous_ppis_with_target_proteins, unique_human_proteins_fi):
	""" Get list of the human proteins from ppis file (these are the proteins that we want to blast against pdb)

	Args: 
		endogenous_ppis_with_target_proteins (tsv file): tsv file with ppi pairs 

	Return: 
		list of all the unique proteins in the given ppi file (regardless of whether they appear in column1/column2)
	"""
	list_of_unique_proteins = []
	with open(endogenous_ppis_with_target_proteins) as ppis:
		for line in ppis:
			line = line.strip("\n")
			info = line.split("\t")

			protA = info[0]
			protB = info[1]

			if protA not in list_of_unique_proteins:
				list_of_unique_proteins.append(protA)
			if protB not in list_of_unique_proteins:
				list_of_unique_proteins.append(protB)
	print(f'Number of proteins to map to uniprot: {len(list_of_unique_proteins)}')

	#2023--OUTPUT the list!
	#Then feed into uniprot to get the sequence and reviewed status
	with open(unique_human_proteins_fi, 'w') as human_fi:
		human_fi.write("\n".join(list_of_unique_proteins))
	






def main():
	script_dir = osp.dirname(__file__)
	intact_fi = osp.join(script_dir, '..', 'data', 'IntAct_human_human', 'intact.txt')
	all_human_ppis_fi = osp.join(script_dir, '..', 'data', 'IntAct_human_human', 'intact_all_physical_human_ppis.tsv')

	exogenous_ppis = osp.join(script_dir, '..', 'data', 'combined_homology_and_pdb_anno_human_virus_interface_residues.tsv') #using all target proteins regardless of whether it is part of an antibody/mhc complex
	
	list_of_human_target_proteins = osp.join(script_dir, '..', 'data', 'IntAct_human_human', 'list_of_target_human_proteins.txt')

	endogenous_ppis_with_target_proteins = osp.join(script_dir, '..', 'data', 'IntAct_human_human', 'intact_endogenous_human_ppis.tsv')
	

	unique_human_proteins_fi = osp.join(script_dir, '..', 'data', 'IntAct_human_human', "unique_human_list.txt")


	print("##########    Processing IntAct human ppi data     ##########")
	if not osp.exists(all_human_ppis_fi):
		format_intact(intact_fi, all_human_ppis_fi)

	#print("##########    Make list of PPIs with at least 2 publication hits     ##########")
	#list_of_2_hit_ppis = get_only_human_ppis_with_at_least_2_pubs(all_human_ppis_fi)

	print("##########    Extracting endogenous ppi data for target proteins     ##########")
	#get_endogenous_ppis(exogenous_ppis, list_of_human_target_proteins, list_of_2_hit_ppis, all_human_ppis_fi, endogenous_ppis_with_target_proteins)
	get_endogenous_ppis(exogenous_ppis, list_of_human_target_proteins, all_human_ppis_fi, endogenous_ppis_with_target_proteins)
	
	print("##########    Getting list of human proteins to retrieve Uniprot info for     ##########")
	get_list_of_human_proteins_from_ppis(endogenous_ppis_with_target_proteins, unique_human_proteins_fi)
	




if __name__ == '__main__':
	main()


