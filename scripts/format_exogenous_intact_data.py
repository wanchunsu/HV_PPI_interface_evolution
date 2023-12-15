'''

This script processes human-virus PPI data from IntAct.
=======================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp

def format_intact(intact_fi, output_intact_formatted):
	""" Format intact file and isolate only columns of interest
	"""
	with open(intact_fi, 'r') as intact:
		next(intact)

		with open(output_intact_formatted, 'wt') as output_fi:
			output_fi.write("interactorA\tinteractorB\ttaxidA\ttaxidB\tdetection_method\tinteraction_type\n")
			for line in intact:
			
				line = line.strip("\n")
				info = line.split("\t")

				interactorA = info[0].replace("uniprotkb:", "")
				interactorB = info[1].replace("uniprotkb:", "")

				detection_method = info[6]
				taxidA = info[9]
				taxidB = info[10]
				
				interaction_type = info[11].split("(")[1].replace(")", "")


				output_fi.write(interactorA + '\t' + interactorB + '\t' + taxidA + '\t' + taxidB + '\t' + detection_method + '\t' + interaction_type + "\n")

def get_human_physical_interactions(formatted_intact, fully_filtered_intact):
	'''
	Get only interactions where one partner is human and is a "direct interaction", "physical association" or "association"
	'''
	with open(formatted_intact, 'r') as intact_formatted:
		next(intact_formatted)
		with open(fully_filtered_intact, 'wt') as output_fully_filtered:
			output_fully_filtered.write("interactorA\tinteractorB\ttaxidA\ttaxidB\tdetection_method\tinteraction_type\n")
			for line in intact_formatted:
				line = line.strip("\n")
				info = line.split("\t")
				interactorA = info[0]
				interactorB = info[1]
				taxidA = info[2]
				taxidB = info[3]
				detection_method = info[4]
				interaction_type = info[5]
				A_or_B_is_human = "invalid"
				if interaction_type == "direct interaction" or interaction_type == "physical association" or interaction_type == 'association': #also added in "association"
					if 'taxid:9606' in taxidA and  'taxid:9606' not in taxidB: #not both human (not even sure why this 'virus' dataset has human-human interactions)
						A_or_B_is_human = 'A'
					elif 'taxid:9606' in taxidB and 'taxid:9606' not in taxidA:
						A_or_B_is_human = 'B'

					if 'taxid:-1' in taxidA or 'taxid:-2'in taxidA  or 'taxid:-1'in taxidB  or 'taxid:-2'in taxidB : #these are synthesized proteins, not the real thing, so filter off
						continue
					if 'intact' in interactorA or 'intact' in interactorB: #filter off ones that don't have uniprot ID (i.e. they only have intact id and are not proteins)
						continue

				if A_or_B_is_human == 'A': #human comes first, so output as is
					output_fully_filtered.write(line + "\n")

				elif A_or_B_is_human == 'B': #virus comes first, so put human first
					output_fully_filtered.write(interactorB + "\t" + interactorA + "\t" + taxidB + "\t" + taxidA + "\t" +  detection_method + '\t' + interaction_type + "\n")
				
				else: #not physical or either both or none of the partners are human, so filter off
					continue

def get_list_of_non_virus_binding_partners(intact_formatted_filtered, list_of_all_virus_taxids):
	""" There are some false positive virus data (e.g. human-mouse ppi annotated as a virus ppi entry), so compile a list of the human binding partners that aren't viruses
	Args: 
		intact_formatted_filtered: processed intact data with possible false positive virus data
		list_of_all_virus_taxids: list of all taxonomy ids of all viruses
	Return: 
		list of taxonomy info that don't correspond to viruses
	"""
	tax_v = open(list_of_all_virus_taxids).readlines()
	tax_v = [t.strip("\n") for t in tax_v]
	
	bp_not_virus = []
	with open(intact_formatted_filtered, 'r') as input_fi:
		next(input_fi) #skip header
	
		for line in input_fi:
			taxidB = line.split("\t")[3]
			taxidB_short = taxidB.split("|")[0].split("(")[0].replace("taxid:", "") #get just the number

			if taxidB_short not in tax_v:#this taxid does not correspond to a virus
				if taxidB not in bp_not_virus:  #append if not already in list
					bp_not_virus.append(taxidB)
	return bp_not_virus
						


def filter_off_non_virus_binding_partners(intact_formatted_filtered, list_of_non_virus, final_filtered_data):
	""" Make new intact data with ppis involving human and non-virus partners filtered off

	Args:
		intact_formatted_filtered: processed intact data with possible false positive virus data
		list_of_non_virus: list of taxonomy info that don't correspond to viruses
		final_filtered_data: output data file
	"""
	pairs_seen = set() # holds lines already seen (need to filter off duplicate lines)
	with open(intact_formatted_filtered, 'r') as input_fi:
		next(input_fi) #skip header
		
		with open(final_filtered_data, 'wt') as final_filtered_o:

			final_filtered_o.write("interactorA\tinteractorB\ttaxidA\ttaxidB\n")
			for line in input_fi:
				pair_of_interactors = line.split("\t")[0] + "|" + line.split("\t")[1]
				if pair_of_interactors not in pairs_seen: # not a duplicate

					pairs_seen.add(pair_of_interactors)
					taxidB = line.split("\t")[3]
					#print(taxidB)
					if taxidB in list_of_non_virus: #this means the interacting partner is not a virus (filter off)
						# print(taxidB)
						continue
					else:
						interactorA = line.split("\t")[0]
						interactorB = line.split("\t")[1]
						taxidA = line.split("\t")[2]
						taxidB = line.split("\t")[3]

						final_filtered_o.write(interactorA + "\t" + interactorB + "\t" + taxidA + "\t" + taxidB +"\n")


def get_exogenous_ppis(final_filtered_data, exogenous_ppis_out_fi):
	""" Make file with unique human-virus ppis where isoforms and proteins segments are replaced by the canonical protein and entries with non-uniprot ids are deleted
	Args:
		final_filtered_data: processed data file
		exogenous_ppis_out_fi: output file to store the exogenous ppis
	"""
	list_of_human_virus_pairs = [] #to store unique ppis, so that we don't have duplicates
	with open(final_filtered_data, 'r') as final_filtered_d:
		next(final_filtered_d)
		with open(exogenous_ppis_out_fi, 'w') as outfi:
			for line in final_filtered_d:
				human_prot = line.split("\t")[0]
				virus_prot = line.split("\t")[1]
				
				#replace isoform ids with canonical id, and for protein segements (ending with -PRO...), just take the whole protein id
				human_prot = human_prot.split("-")[0]
				virus_prot = virus_prot.split("-")[0]
				
				if len(human_prot) != 6 or len(virus_prot) != 6: #not a uniprot id
					continue
				if (human_prot, virus_prot) not in list_of_human_virus_pairs:
					list_of_human_virus_pairs.append((human_prot, virus_prot))
					outfi.write(human_prot + "\t" + virus_prot + "\n")

def get_list_of_human_and_virus_proteins_from_ppis(exogenous_ppis_out_fi, unique_human_proteins_fi, unique_virus_proteins_fi):
	""" Get lists of the human and virus proteins from ppis file (these are the proteins that we want to blast against pdb)

	Args: 
		exogenous_ppis_out_fi (tsv file): tsv file with ppi pairs 

	Return: 
		list of all the unique human proteins in the given ppi file 
		list of all the unique virus proteins in the given ppi file
	"""
	list_of_unique_human_proteins = []
	list_of_unique_virus_proteins = []
	with open(exogenous_ppis_out_fi) as ppis:
		for line in ppis:
			line = line.strip("\n")
			info = line.split("\t")

			human_p = info[0]
			virus_p = info[1]

			if human_p not in list_of_unique_human_proteins:
				list_of_unique_human_proteins.append(human_p)
			if virus_p not in list_of_unique_virus_proteins:
				list_of_unique_virus_proteins.append(virus_p)
	print(f'Number of human proteins to map to uniprot: {len(list_of_unique_human_proteins)}')
	print(f'Number of virus proteins to map to uniprot: {len(list_of_unique_virus_proteins)}')
	
	# Output the two lists
	# Will then feed into uniprot to get the sequence and reviewed status
	with open(unique_human_proteins_fi, 'w') as human_fi, open(unique_virus_proteins_fi, 'w') as virus_fi:
		human_fi.write("\n".join(list_of_unique_human_proteins))
		virus_fi.write("\n".join(list_of_unique_virus_proteins))
	

def main():
	script_dir = os.path.dirname(__file__)
	intact_data_dir = os.path.join(script_dir, '..', 'data', 'IntAct_human_virus')

	intact_data = osp.join(intact_data_dir, "IntAct_virus_ppis.txt")
	intact_formatted_fi = osp.join(intact_data_dir, "IntAct_virus_ppis.formatted.txt")
	intact_formatted_filtered_fi = osp.join(intact_data_dir, "IntAct_virus_ppis.formatted.filtered.txt")

	list_of_all_virus_taxids = osp.join(intact_data_dir, "all_virus_taxonomy.txt")
	intact_final_filtered_data = osp.join(intact_data_dir, "IntAct_virus_ppis.formatted.filtered_final.txt")

	exogenous_ppis_out_fi = osp.join(intact_data_dir, "intact_exogenous_human_ppis.tsv")

	
	unique_human_proteins_fi = osp.join(intact_data_dir, "unique_human_list.txt")
	unique_virus_proteins_fi = osp.join(intact_data_dir, "unique_virus_list.txt")


	format_intact(intact_data, intact_formatted_fi) #first format intact data 
	get_human_physical_interactions(intact_formatted_fi, intact_formatted_filtered_fi) #get only physical interactions
	list_of_non_virus = get_list_of_non_virus_binding_partners(intact_formatted_filtered_fi, list_of_all_virus_taxids) #get list of binding partners with human proteins (for some reason some of these are not viruses even though we downloaded these interactions form the virus database)
	filter_off_non_virus_binding_partners(intact_formatted_filtered_fi, list_of_non_virus, intact_final_filtered_data)
	get_exogenous_ppis(intact_final_filtered_data, exogenous_ppis_out_fi)
	
	get_list_of_human_and_virus_proteins_from_ppis(exogenous_ppis_out_fi, unique_human_proteins_fi, unique_virus_proteins_fi)
	

	

if __name__ == "__main__":
	main()


