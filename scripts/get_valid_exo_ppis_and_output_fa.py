'''

Get valid exogenous PPIs (both interacting partners need to be reviewed and valid)
and output their fasta sequences.
==================================================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''
import os.path as osp
import os
from uniprot_tools import load_uniprot_info, obsolete_or_alt_ids2, get_unreviewed2, write_uniprot_seqs_to_fasta

def load_proteins_list(protein_list_fi):
	protein_ids = []
	with open(protein_list_fi) as p:
		for line in p:
			line = line.strip("\n")
			protein_ids.append(line)
	return protein_ids

def retrieve_uniprot_seqs(human_protein_fi, virus_protein_fi, human_uniprot_tsv, virus_uniprot_tsv, exogenous_ppis_out_fi, valid_ids_exogenous_ppis_out_fi, out_fa_human, out_fa_virus):
	""" Retrieve the uniprot sequences for all the human proteins involved in the human ppis of interest

	Args:
		exogenous_ppis_out_fi (tsv file): tsv file with ppi pairs
		valid_ids_exogenous_ppis_out_fi: new tsv file without any obsolete or alternative ids
		out_fa_human (file): human proteins fasta file to output
		out_fa_virus (file): virus proteins fasta file to output

	"""
	#get list of uniprot ids
	print("##########    Getting list of uniprot ids     ##########")
	human_uniprot_ids = load_proteins_list(human_protein_fi)
	virus_uniprot_ids = load_proteins_list(virus_protein_fi)
	
	print("##########    Retrieving info for uniprot ids     ##########")
	# loading mapping df from uniprot info files
	mapping_df_human = load_uniprot_info(human_uniprot_tsv) 
	mapping_df_virus = load_uniprot_info(virus_uniprot_tsv)
	
	# Get obsolete or alternative ids (will delete interactions containing them)
	obs_or_alt_ids_human = obsolete_or_alt_ids2(human_uniprot_ids, mapping_df_human)
	obs_or_alt_ids_virus = obsolete_or_alt_ids2(virus_uniprot_ids, mapping_df_virus)
	print(f'obs_or_alt_humans: {obs_or_alt_ids_human}')
	print(f'obs_or_alt_virus: {obs_or_alt_ids_virus}')

	#get unreviewed uniprot ids (we don't want unreviewed sequences, so will remove)
	unreviewed_human = get_unreviewed2(mapping_df_human)
	unreviewed_virus = get_unreviewed2(mapping_df_virus)
	print(f'unreviewed_humans: {unreviewed_human}')
	print(f'unreviewed_virus: {unreviewed_virus}')

	list_of_original_exogenous = [] #make list to store the exogenous ppis first
	with open(exogenous_ppis_out_fi) as exo_ppis:
		for line in exo_ppis:
			line = line.strip("\n")
			info = line.split("\t")
			list_of_original_exogenous.append((info[0], info[1])) 

	#new list to remove or append to b/c we'll be traversing the list and can't modify it while traversing
	list_of_original_exogenous_new = list_of_original_exogenous.copy() 
	
	for pair in list_of_original_exogenous:
		human_prot = pair[0]
		virus_prot = pair[1]

		#2023 update, instead of mapping alternate/obsolete ids to new ones, we'll just remove entries with obsolete/alternate ids --sequence & status will be denoted NaN
		if human_prot in obs_or_alt_ids_human or virus_prot in obs_or_alt_ids_virus:
			list_of_original_exogenous_new.remove(pair)
			continue

		#delete any interactions with unreviewed partners
		if human_prot in unreviewed_human or virus_prot in unreviewed_virus: #if either protein is unreviewed, we remove this pair
			list_of_original_exogenous_new.remove(pair)
			continue

	unique_human = []
	unique_virus = []
	with open(valid_ids_exogenous_ppis_out_fi, 'w') as new_ppis_fi:
		new_ppis_fi.write("human_prot" + "\t" + "virus_prot" + "\n")
		for p in list_of_original_exogenous_new:
			human_prot = p[0]
			virus_prot = p[1]
			if human_prot not in unique_human:
				unique_human.append(human_prot)
			if virus_prot not in unique_virus:
				unique_virus.append(virus_prot)

			new_ppis_fi.write(human_prot + "\t" + virus_prot + "\n")

	# make list of unique human prot
	# write the contents of the pandas dataframes into fasta files (making sure to delete any entries with unreviewed/obsolete/alternate proteins)
	mapping_df_relevant_human = mapping_df_human[mapping_df_human.From.isin(unique_human)] #drop invalid (unreviewed) entries, as well as those whose partners are all invalid
	mapping_df_relevant_virus = mapping_df_virus[mapping_df_virus.From.isin(unique_virus)] #drop invalid (unreviewed) entries, as well as those whose partners are all invalid

	mapping_df_relevant_human = mapping_df_relevant_human[mapping_df_relevant_human.From == mapping_df_relevant_human.Entry] #drop obsolete/alternate entries
	mapping_df_relevant_virus = mapping_df_relevant_virus[mapping_df_relevant_virus.From == mapping_df_relevant_virus.Entry] #drop obsolete/alternate entries

	mapping_df_relevant_human = mapping_df_relevant_human[["From", "Sequence"]] #don't need entry & reviewed column anymore because we've already taken into account any alternative or obsolete ids and unreviewed
	mapping_df_relevant_virus = mapping_df_relevant_virus[["From", "Sequence"]] #don't need entry & reviewed column anymore because we've already taken into account any alternative or obsolete ids and unreviewed

	mapping_df_relevant_human = mapping_df_relevant_human.drop_duplicates()
	mapping_df_relevant_virus = mapping_df_relevant_virus.drop_duplicates()

	print(f'Number of human proteins: {len(unique_human)}')
	print(f'Number of virus proteins: {len(unique_virus)}')

	print("##########    Writing human sequences to fasta file     ##########")
	write_uniprot_seqs_to_fasta(mapping_df_relevant_human, out_fa_human)
	print("##########    Writing virus sequences to fasta file     ##########")
	write_uniprot_seqs_to_fasta(mapping_df_relevant_virus, out_fa_virus)
	print("Done!")
	
def main():
	script_dir = os.path.dirname(__file__)
	intact_data_dir = os.path.join(script_dir, '..', 'data', 'IntAct_human_virus')

	human_protein_fi = osp.join(intact_data_dir, "unique_human_list.txt") #list of unique human proteins to map
	virus_protein_fi = osp.join(intact_data_dir, "unique_virus_list.txt") #list of unique virus proteins to map

	human_uniprot_tsv = osp.join(intact_data_dir, "human_uniprot.tsv")# uniprot info downloaded for human_protein_fi
	virus_uniprot_tsv = osp.join(intact_data_dir, "virus_uniprot.tsv")# uniprot info downloaded for virus_protein_fi
	
	exogenous_ppis_out_fi = osp.join(intact_data_dir, "intact_exogenous_human_ppis.tsv")

	valid_ids_exogenous_ppis_out_fi = osp.join(intact_data_dir, "intact_exogenous_human_ppis_with_valid_ids.tsv")
	out_fa_human =  osp.join(intact_data_dir, 'fasta_files', 'human_uniprot_sequences.fasta')
	out_fa_virus = osp.join(intact_data_dir, 'fasta_files', 'virus_uniprot_sequences.fasta')
	if not osp.exists(osp.dirname(out_fa_human)):
		os.mkdir(osp.dirname(out_fa_human))

	retrieve_uniprot_seqs(human_protein_fi, virus_protein_fi, human_uniprot_tsv, virus_uniprot_tsv, exogenous_ppis_out_fi, valid_ids_exogenous_ppis_out_fi, out_fa_human, out_fa_virus)
	
if __name__ == '__main__':
	main()
