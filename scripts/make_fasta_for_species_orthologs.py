'''

Make individual fasta files for orthologs
===========================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp
import json 
from fasta_tools import parse_fasta_file_into_dict

def get_orthologs(list_of_targ_prots_of_interest, rbh_folder): 
	# make dict of corresponding orthologs for each target protein (target protein:{horse: horse_ortholog_id, cow: cow_ortholog_id...})
	organism_list = os.listdir(rbh_folder)
	target_prots_and_their_orthologs = {}
	for org in organism_list:
		json_file_with_orthologs_and_aligned_pos_seq = osp.join(rbh_folder, org, 'rbh_orthologs_with_aligned_pos_and_seq.json')
		json_fi_with_orthologs = {}

		with open(json_file_with_orthologs_and_aligned_pos_seq) as ij:
			json_fi_with_orthologs = json.load(ij)

		for targ_prot in json_fi_with_orthologs:
			if targ_prot not in list_of_targ_prots_of_interest: #we only want to get the orthologs for the targ_prots of interest (i.e. targ_prots in v_target_h that aren't immunoglobulins)
				continue
			ortholog = json_fi_with_orthologs[targ_prot][0]

			if targ_prot not in target_prots_and_their_orthologs:
				target_prots_and_their_orthologs[targ_prot] = {}
			target_prots_and_their_orthologs[targ_prot][org] = ortholog

	return target_prots_and_their_orthologs

def retrieve_orthologous_species_proteomes_dict(proteomes_folder, organism_name):
	""" Format specified organism proteome fasta file into a dicitonary (keeping only uniprot id as header -- dict key and sequence as value)
	Args:
		proteomes_folder = folder containing proteomes 
		organism_name: name of organism who's proteome we want to retrieve

	"""
	organism_proteome_fi = osp.join(proteomes_folder, organism_name + '_proteome_all.fasta')
	proteome_fasta_dict = parse_fasta_file_into_dict(organism_proteome_fi, change_header=True, split_header_by='|', ind_to_keep = 1)
	return proteome_fasta_dict


def make_fasta_files(proteomes_folder, target_prots_and_their_orthologs, organism_name, species_orthologs_folder):
	""" Output fasta seqs for the orthologous proteins in our specified species
	proteomes_folder: dir to retrieve proteomes from
	target_prots_and_their_orthologs: dict of target proteins and their corresponding orthologs in different species
	organism_name: common name of species that we're interested in
	species_orthologs_folder: output folder for our orthologous protein seqs (fasta file) from our species of interest


	"""
	species_proteomes = retrieve_orthologous_species_proteomes_dict(proteomes_folder, organism_name)

	with open(osp.join(species_orthologs_folder, organism_name + '_orthologs.fasta'), 'w') as outfi:
		#>targ_prot|species_ortholog (note we want to keep the target protein id there for later analyses)
		for targ_prot in target_prots_and_their_orthologs:
			ortholog = target_prots_and_their_orthologs[targ_prot][organism_name]
			ortholog_seq = species_proteomes[ortholog]
			outfi.write('>'+ targ_prot + '|' + ortholog + "\n") #header is the human target protein and the orthologous protein id (sep = '|')
			outfi.write(ortholog_seq + "\n")

def main():
	script_dir = osp.dirname(__file__)
	folder_with_targ_prots_of_interest = osp.join(script_dir, '..', 'data', 'dnds', 'human_targ_prot_fasta')
	rbh_folder = osp.join('..', 'data', 'homologs_and_conservation', 'orthologs', 'v_target_h', 'rbh')
	proteomes_folder = osp.join('..', 'data', 'homologs_and_conservation', 'orthologs', 'proteomes')
	species_orthologs_folder = osp.join('..', 'data', 'homologs_and_conservation', 'orthologs', 'species_ortholog_prots_fasta')
	
	if not osp.exists(species_orthologs_folder):
		os.mkdir(species_orthologs_folder)


	targ_prots_of_interest_fis = os.listdir(folder_with_targ_prots_of_interest)##list of targ prots in v_target_h that aren't immunoglobulins
	list_of_targ_prots_of_interest = [targ_prots_fi.split('.')[0] for targ_prots_fi in targ_prots_of_interest_fis]
	
	target_prots_and_their_orthologs = get_orthologs(list_of_targ_prots_of_interest, rbh_folder)

	for org in os.listdir(rbh_folder):
		
		#output orthologous protein seqs into a fasta file (one file per species)
		print(f'Making orthologous proteins fasta file for {org} . . .')
		make_fasta_files(proteomes_folder, target_prots_and_their_orthologs, org, species_orthologs_folder) 

if __name__ == '__main__':
	main()