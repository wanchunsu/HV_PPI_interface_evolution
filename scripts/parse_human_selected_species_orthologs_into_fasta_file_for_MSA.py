'''

Parse human target proteins and orthologs into fasta files for creating MSAs
============================================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''


import os
import os.path as osp
import argparse
import json 
from Bio import SeqIO
from fasta_tools import parse_fasta_file_into_dict


def parse_ortholog_aligned_seqs(dir_for_organisms):
	""" Parse the specified orthologs for each protein and store them into a nested dictionary with targ_prot as the outer key 
	and the ortholog_species_name as the value in a list (one element per species ortholog)

	Args:
		dir_for_organisms: directory of species data (where orthologs are stored)
		
	Return:
		target_prots_and_their_orthologs: dictionary storing the above mentioned information
		


	"""
	

	organism_list = os.listdir(dir_for_organisms)
	
	target_prots_and_their_orthologs = {} # in the following format ({targ_prot: {mouse_ortholog_aligned_pos: aligned_seq_mouse, chimp_orthologs, ...}, targ_prot2: [mouse_ortholog, chimp_ortholog...], ...})
	for org in organism_list:
		json_file_with_orthologs_and_aligned_pos_seq = osp.join(dir_for_organisms, org, 'rbh_orthologs_with_aligned_pos_and_seq.json')
		json_fi_with_orthologs = {}

		with open(json_file_with_orthologs_and_aligned_pos_seq) as ij:
			json_fi_with_orthologs = json.load(ij)

		for targ_prot in json_fi_with_orthologs:
			ortholog = json_fi_with_orthologs[targ_prot][0]
			position = str(json_fi_with_orthologs[targ_prot][1][0]) + '_' + str(json_fi_with_orthologs[targ_prot][1][1])
			fasta_header = ortholog + '_' + org + '_' + position
			fasta_aligned_seq = json_fi_with_orthologs[targ_prot][2].replace('-', '') #remove any gaps in the aligned seqs 
			
			if targ_prot not in target_prots_and_their_orthologs:
				target_prots_and_their_orthologs[targ_prot] = {}
			target_prots_and_their_orthologs[targ_prot][fasta_header] = fasta_aligned_seq
			

	#print(target_prots_and_their_orthologs)
	#print(dict_of_ortholog_proteome)
	return target_prots_and_their_orthologs


def output_human_orthologs_fasta_fi_to_align(targ_prot_fasta_folder, target_prots_and_their_orthologs, output_dir):
	""" Output the target prot and ortholog (region aligned to targ prot) to perform MSA on

	Args:
		targ_prot_fasta_folder: folder containing target protein fasta seqs
		target_prots_and_their_orthologs: dictionary storing target proteins and their orthologs (as produced from parse_ortholog_aligned_seqs function above)
		output_dir: directory for outputting the sequences to align via MSA (we'll have one fasta file per target protein)

	"""
	for targ_prot in target_prots_and_their_orthologs:
		fasta = SeqIO.parse(open(osp.join(targ_prot_fasta_folder, targ_prot + '.fasta')), 'fasta')
		targ_prot_sequence = ''
		for fa in fasta:
		  	name, targ_prot_sequence = fa.id, str(fa.seq)


		fasta_file_name = osp.join(output_dir, targ_prot + '_to_align.fasta')
		with open(fasta_file_name, 'w') as fa_out:
			fa_out.write(">" + targ_prot + "\n" + targ_prot_sequence + "\n") #write the targ prot as first fasta sequence (this will be the reference for when we compute the rate4site)
			for ortholog_species_pos in target_prots_and_their_orthologs[targ_prot]:
				ortholog_aligned_seq = target_prots_and_their_orthologs[targ_prot][ortholog_species_pos]
				fa_out.write(">" + ortholog_species_pos + "\n" + ortholog_aligned_seq + "\n")

		print(f'Finished making fasta files for MSA for target prot: {targ_prot}')


def main():
	script_dir = osp.dirname(__file__)
	parser = argparse.ArgumentParser()

	parser.add_argument('-f', '--folder')
	# parser.add_argument('-o', '--organisms')
	parser.add_argument('-i', '--output_id')
	parser.add_argument('-o', '--ortholog_folder')
	args = parser.parse_args()

	folder_name = args.folder #'all', 'h_target_v', 'v_target_h'

	output_folder_id = args.output_id
	ortholog_folder_label = args.ortholog_folder
	# organism_names = args.organisms # 'org1_org2_....' #list of organisms whose orthologs we want
	
	print(f'\n #####   Running for {folder_name} ##### ')
	dir_for_organisms = osp.join('..', 'data', 'homologs_and_conservation', ortholog_folder_label, folder_name, 'rbh')

	
	targ_prot_fasta_folder = osp.join('..', 'data', 'exo_and_endo', 'human_target_fasta_folder')

	
	output_dir = osp.join('..', 'data', 'homologs_and_conservation', ortholog_folder_label, folder_name, 'fasta_fis_for_MSA_' + output_folder_id)

	if not osp.exists(output_dir):
		os.mkdir(output_dir)

	target_prots_and_their_orthologs = parse_ortholog_aligned_seqs(dir_for_organisms)	

	output_human_orthologs_fasta_fi_to_align(targ_prot_fasta_folder, target_prots_and_their_orthologs, output_dir)

	
if __name__ == '__main__':
	main()