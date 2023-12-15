'''

Get valid endgenous ppis and their corresponding sequences
=========================================================
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


def retrieve_uniprot_seqs(human_protein_fi, human_uniprot_tsv, endogenous_ppis_with_target_proteins, valid_ids_endogenous_ppis_with_target_proteins, out_fa):
	""" Retrieve the uniprot sequences for all the human proteins involved in the human ppis of interest

	Args:
		endogenous_ppis_with_target_proteins (tsv file): tsv file with ppi pairs
		valid_ids_endogenous_ppis_with_target_proteins: new tsv file without any obsolete or alternative ids
		out_fa (file): fasta file to output

	"""
	#get list of uniprot ids
	print("##########    Getting list of uniprot ids     ##########")
	uniprot_ids = load_proteins_list(human_protein_fi)
	#print("\n".join(uniprot_ids))
	
	
	#import function from uniprot_tools to get the sequences using requests
	print("##########    Retrieving sequences from uniprot ids     ##########")
	mapping_df = load_uniprot_info(human_uniprot_tsv)
	#print(mapping_df)
	

	# Get obsolete or alternative ids (will delete interactions containing them)
	obs_or_alt_ids = obsolete_or_alt_ids2(uniprot_ids, mapping_df) #get obsolete or alternative id mappings
	print(f'obs_or_alt: {obs_or_alt_ids}')

	#get unreviewed uniprot ids (we don't want unreviewed sequences, so will remove)
	unreviewed = get_unreviewed2(mapping_df)
	print(f'unreviewed: {unreviewed}')

	list_of_original_endogenous =[]

	with open(endogenous_ppis_with_target_proteins, 'r') as endo_ppis:
		for line in endo_ppis:
			line = line.strip("\n")
			info = line.split("\t")

			list_of_original_endogenous.append((info[0], info[1]))

	#new list to remove or append to b/c we'll be traversing the list and can't modify it while traversing
	list_of_original_endogenous_new = list_of_original_endogenous.copy() #new list to remove or append to
	
	for pair in list_of_original_endogenous:
		prot1 = pair[0]
		prot2 = pair[1]



		
		#2023 update, instead of mapping alternate/obsolete ids to new ones, we'll just remove entries with obsolete/alternate ids --sequence & status will be denoted NaN
		if prot1 in obs_or_alt_ids or prot2 in obs_or_alt_ids:
			list_of_original_endogenous_new.remove(pair)
			continue
		
		
		#delete any interactions with unreviewed partners
		if prot1 in unreviewed or prot2 in unreviewed: #if either protein is unreviewed, we remove this pair
			list_of_original_endogenous_new.remove(pair)
			continue
	
	count_unique_human_targs = []
	count_unique_human_partners = []

	with open(valid_ids_endogenous_ppis_with_target_proteins, 'w') as new_ppis_fi:
		new_ppis_fi.write("target_prot" + "\t" + "partner_prot" + "\n")
		for p in list_of_original_endogenous_new:
			prot1 = p[0]
			prot2 = p[1]

			new_ppis_fi.write(prot1 + "\t" + prot2 + "\n")

			#count unique number of human targets and human partners
			if prot1 not in count_unique_human_targs:
				count_unique_human_targs.append(prot1)
			if prot2 not in count_unique_human_partners:
				count_unique_human_partners.append(prot2)
	print(f'Unique number of human target proteins: {len(count_unique_human_targs)}')
	print(f'Unique number of human partner proteins: {len(count_unique_human_partners)}')
	
	# write the contents of the pandas dataframe into a fasta file (making sure to delete any entries with unreviewed/obsolete/alternate proteins)
	mapping_df_relevant = mapping_df[mapping_df.Reviewed == 'reviewed'] #drop unreviewed entries

	mapping_df_relevant = mapping_df_relevant[mapping_df_relevant.From == mapping_df_relevant.Entry] #drop obsolete/alternate entries

	mapping_df_relevant = mapping_df_relevant[["From", "Sequence"]] #don't need entry & reviewed column anymore because we've already taken into account any alternative or obsolete ids and unreviewed

	mapping_df_relevant = mapping_df_relevant.drop_duplicates()

	print(f'Total number of human proteins: {len(mapping_df_relevant.index)}')
	
	print("##########    Writing sequences to fasta file     ##########")
	write_uniprot_seqs_to_fasta(mapping_df_relevant, out_fa)

	print("Done!")

def main():
	script_dir = osp.dirname(__file__)
	intact_data_dir = osp.join(script_dir, '..', 'data', 'IntAct_human_human')
	
	human_protein_fi = osp.join(intact_data_dir, "unique_human_list.txt") #list of unique human proteins to map
	human_uniprot_tsv = osp.join(intact_data_dir, "human_uniprot.tsv") # uniprot info downloaded for human_protein_fi

	endogenous_ppis_with_target_proteins = osp.join(intact_data_dir, 'intact_endogenous_human_ppis.tsv')

	valid_ids_endogenous_ppis_with_target_proteins = osp.join(intact_data_dir, 'intact_endogenous_human_ppis_with_valid_ids.tsv')
	fasta_outfile = osp.join(intact_data_dir, 'fasta_files', 'human_uniprot_sequences.fasta')
	if not osp.exists(osp.dirname(fasta_outfile)):
		os.mkdir(osp.dirname(fasta_outfile))

	retrieve_uniprot_seqs(human_protein_fi, human_uniprot_tsv, endogenous_ppis_with_target_proteins, valid_ids_endogenous_ppis_with_target_proteins, fasta_outfile)


if __name__ == '__main__':
	main()