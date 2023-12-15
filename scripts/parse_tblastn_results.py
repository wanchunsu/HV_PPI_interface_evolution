'''

Parse tBLASTn results to get the CDS sequences for human 
target proteins and their orthologs.
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
import os
from fasta_tools import parse_fasta_file_into_dict, retrieve_first_fasta_seq_from_fi


def parse_tblastn_results(blast_result_fi):
	""" Parse blast results file into a dictionary storing the best cds hit
	blast_result_fi: file containing blast results for protein against cds tblastn

	Return:
		dict_of_results: dctionary storing best hits for each protein

	"""
	dict_of_results = {}
	with open(blast_result_fi) as blast_res:
		for line in blast_res:
			line = line.strip("\n")
			info = line.split("\t")
			qseqid = info[0]
			start_pos_cds = int(info[8])
			end_pos_cds = int(info[9])
			cds_id = info[-1]

			start_pos_prot = int(info[6])
			end_pos_prot = int(info[7])
			qseq_aligned = info[12]
			sseq_translated_aligned = info[13]
			pident = float(info[2])

			if qseqid not in dict_of_results: #store only the first one
				dict_of_results[qseqid] = [cds_id, start_pos_prot-1, end_pos_prot, start_pos_cds-1, end_pos_cds, qseq_aligned, sseq_translated_aligned, pident] #start_pos_cds-1,end_pos_cds are the indices that will give us the substring at start_pos_cds (inclusive-->end_pos_cds substring (inclusive )
	return dict_of_results	
	

def parse_and_check_human_tblastn_results(human_blast_results_folder, blast_results_folder, human_target_proteins):
	""" Parse human against refseq cds seq tblastn results. 
	If identity is <100%, then we check if human targ against ensembl cds seqs give better results.
	If identity is == 100%, we also check if alignment length is same as the protein itself, if not, we check if ensembl cds has longer alignment (given that it also have 100% identity)
	Args:
		blast_results_folder: folder containing all tblastn results
		human_target_proteins_fasta: fasta file containing sequences of our human target proteins
		prot_fasta_folder: folder containing 

	"""
	human_against_ensembl_cds_blast_results = parse_tblastn_results(osp.join(blast_results_folder, 'human.prot_against_cds.blast.out'))

	dict_of_human_results = {}
	for res in os.listdir(human_blast_results_folder):
		with open(osp.join(human_blast_results_folder, res)) as blast_res:
			for line in blast_res:
				line = line.strip("\n")
				info = line.split("\t")
				qseqid = info[0]
				start_pos_cds = int(info[8])
				end_pos_cds = int(info[9])
				cds_id = "refseq" #this will be to distinguish between ensembl and refseq cds

				start_pos_prot = int(info[6])
				end_pos_prot = int(info[7])
				qseq_aligned = info[12]
				sseq_translated_aligned = info[13]
				pident = float(info[2])

				
				query_aligned_len = end_pos_prot - start_pos_prot + 1
				
				##For cases where the refseq cds doesn't match the protein seq entirely, check if tblastn against ensembl cds gives better coverage(pident) 
				if pident != 100.0:
					#print(f'not identical: {qseqid}')
					#check if ensembl cds identity is 100%

					if human_against_ensembl_cds_blast_results[qseqid][7] > pident:  #if yes, then choose ensembl cds
						#print(f'\tensembl gives better coverage {pident}, {human_against_ensembl_cds_blast_results[qseqid][7]}')
						dict_of_human_results[qseqid] = human_against_ensembl_cds_blast_results[qseqid]
					else: #if not better, then stick with refseq cds
						dict_of_human_results[qseqid] = [cds_id, start_pos_prot-1, end_pos_prot, start_pos_cds-1, end_pos_cds, qseq_aligned, sseq_translated_aligned, pident] #start_pos_cds-1,end_pos_cds are the indices that will give us the substring at start_pos_cds (inclusive-->end_pos_cds substring (inclusive )
				
				else: #100% identity btwn protein seq and refseq cds
					#check if the refseq cds alignment length is the same as the actual protein
					if query_aligned_len != len(human_target_proteins[qseqid]):
						#print(f'{qseqid} not same length')
						ensembl_alignment_len = human_against_ensembl_cds_blast_results[qseqid][2] - human_against_ensembl_cds_blast_results[qseqid][1]
						ensembl_alignment_pident = human_against_ensembl_cds_blast_results[qseqid][7]
						
						if ensembl_alignment_len > query_aligned_len and ensembl_alignment_pident == 100.0: #check that ensembl alignment gives same 100% identity and better coverage
							#print(f'\tensembl gives longer len: {len(human_target_proteins[qseqid])}, {ensembl_alignment_len}, {query_aligned_len}')
							dict_of_human_results[qseqid] = human_against_ensembl_cds_blast_results[qseqid]
						
						else:
							dict_of_human_results[qseqid] = [cds_id, start_pos_prot-1, end_pos_prot, start_pos_cds-1, end_pos_cds, qseq_aligned, sseq_translated_aligned, pident] #start_pos_cds-1,end_pos_cds are the indices that will give us the substring at start_pos_cds (inclusive-->end_pos_cds substring (inclusive )	
					
					else:
						dict_of_human_results[qseqid] = [cds_id, start_pos_prot-1, end_pos_prot, start_pos_cds-1, end_pos_cds, qseq_aligned, sseq_translated_aligned, pident] #start_pos_cds-1,end_pos_cds are the indices that will give us the substring at start_pos_cds (inclusive-->end_pos_cds substring (inclusive )
	return dict_of_human_results	 


def store_human_prots_and_cds_into_dict(human_blast_results_folder, blast_results_folder, human_target_proteins_fasta, refseq_cds_folder, species_cds_dir):
	""" Make dictionaries to store proteins and their cds' for human (need to do this separately from other organisms because we also have refseq cds')
	Args:
		human_blast_results_folder: folder containing all tblastn results for human proteins against refseq cds
		blast_results_folder: folder containing all tblastn results for protein against ensembl cds
		human_target_proteins_fasta: fasta file containing sequences of our human target proteins
		refseq_cds_folder: folder containing refseq cds sequences
	Return:
		human_targ_prots_cds_dict: will store the human target protein ids (keys) and lists containing the protein seq its corresponding cds seq (i.e. [prot_seq, cds_seq])
		
	"""
	print(f"Parsing tblastn results for human . . .")

	# Dictionary to store proteins and cds sequences
	human_targ_prots_cds_dict = {} # targ_protid: [targ_prot_seq, targ_prot_corresponding_cds_seq]

	#Load target proteins
	human_target_proteins = parse_fasta_file_into_dict(human_target_proteins_fasta, change_header=False, split_header_by='|', ind_to_keep = 1)


	dict_of_human_results = parse_and_check_human_tblastn_results(human_blast_results_folder, blast_results_folder, human_target_proteins)

	human_cds_fi = osp.join(species_cds_dir, 'human.cds.fa')
	human_cds_ensembl = parse_fasta_file_into_dict(human_cds_fi, change_header=False, split_header_by='|', ind_to_keep = 1)

	for targ_prot in dict_of_human_results:
			
		cds_id = dict_of_human_results[targ_prot][0]
		start_prot_ind = dict_of_human_results[targ_prot][1]
		end_prot_ind = dict_of_human_results[targ_prot][2]
		start_cds_ind = dict_of_human_results[targ_prot][3]
		end_cds_ind = dict_of_human_results[targ_prot][4]
		qseq_aligned = dict_of_human_results[targ_prot][5]
		sseq_translated_aligned = dict_of_human_results[targ_prot][6] #this is the translation of the cds (may or may not be the same as the original protein seq, will be same if percent identity in tblastn was 100%)
		
		prot_seq = sseq_translated_aligned #this will ensure that we have the exact protein seq that the best mapped nucleotide seq corresponds to once translated (this may or may not be the same as the original protein seq)
		if cds_id == 'refseq': #we're using the refseq cds sequences
			cds_wholeseq = retrieve_first_fasta_seq_from_fi(osp.join(refseq_cds_folder, targ_prot + '.cds.fasta'))
			cds_seq = cds_wholeseq[start_cds_ind:end_cds_ind]
		else:
			cds_seq= human_cds_ensembl[cds_id][start_cds_ind:end_cds_ind] #get the ensembl cds sequence corresponding to the current protein sequence
		
		cds_seq = cds_seq.replace('N', 'X')
		targ_prot_pos = targ_prot + '_' + str(start_prot_ind+1) + '_' + str(end_prot_ind)
		human_targ_prots_cds_dict[targ_prot] = [prot_seq, cds_seq]

	return human_targ_prots_cds_dict


def store_prots_and_cds_into_dicts(blast_results_folder, species_orthologs_prots_dir, species_cds_dir):
	""" Make dictionaries to store proteins and their cds' for species other than human
	Args:
		blast_results_folder: folder containing all tblastn results
		species_orthologs_prots_dir: directory containing fasta files for species orthologs
		species_cds_dir: directory containing all cds seqs for every species in our analysis
		
	Return:
		
		1) ortholog_prots_dict: will store for each target prot id (outer key), a nested dict with the organism name as key and ortholog prot seqid and prot seq as lists (i.e. targ_prot: {org_name: [ortholog_prot_id, ortholog_prot_seq]})
		2) ortholog_dna_dict: will store for each target prot id (outer key), a nested dict with the organism name as key and ortholog prot seqid and cds seq as lists (i.e. targ_prot: {org_name: [ortholog_prot_id, ortholog_dna_cds_seq]})


	"""

	# Dictionaries to store proteins and cds sequences

	ortholog_prots_dict = {} #targ_protid: {org_name: [ortholog_prot_id, ortholog_prot_seq]}
	ortholog_dna_dict = {} #targ_protid: {org_name: [ortholog_prot_id, ortholog_dna_cds_seq]}



	for res in os.listdir(blast_results_folder):
		
		org_name = res.split(".")[0]

		if org_name == 'human': #Parsing of tblastn results for human is done in a separate function (store_human_prots_and_cds_into_dict)
			continue

		print(f"Parsing tblastn results for {org_name} . . .")
		blast_result_fi = osp.join(blast_results_folder, res)
		dict_of_results = parse_tblastn_results(blast_result_fi)

		
		# if org_name != 'human': 

		#Load seqs for ortholog proteins for the current organism
		species_ortholog_prots_fi = osp.join(species_orthologs_prots_dir, org_name + '_orthologs.fasta')
		species_ortholog_prots = parse_fasta_file_into_dict(species_ortholog_prots_fi, change_header=False, split_header_by='|', ind_to_keep = 1)
	
		
		# Load cds dna seqs for the current organism
		species_cds_fi = osp.join(species_cds_dir, org_name + '.cds.fa')
		species_cds = parse_fasta_file_into_dict(species_cds_fi, change_header=False, split_header_by='|', ind_to_keep = 1)

		

		for targ_prot_ortholog in dict_of_results:
				
				
			targ_prot = targ_prot_ortholog.split("|")[0] # for human targ_prot_ortholog and targ_prot are the same thing
			
			cds_id = dict_of_results[targ_prot_ortholog][0]
			start_prot_ind = dict_of_results[targ_prot_ortholog][1]
			end_prot_ind = dict_of_results[targ_prot_ortholog][2]
			start_cds_ind = dict_of_results[targ_prot_ortholog][3]
			end_cds_ind = dict_of_results[targ_prot_ortholog][4]
			qseq_aligned = dict_of_results[targ_prot_ortholog][5]
			sseq_translated_aligned = dict_of_results[targ_prot_ortholog][6] #this is the translation of the cds (may or may not be the same as the original protein seq, will be same if percent identity in tblastn was 100%)
			
			prot_seq = sseq_translated_aligned # this will ensure that we have the exact protein seq that the best mapped nucleotide seq corresponds to once translated (this may or may not be the same as the original protein seq)
			cds_seq= species_cds[cds_id][start_cds_ind:end_cds_ind] # get the cds sequence corresponding to the current protein sequence
			cds_seq = cds_seq.replace('N', 'X')
			
			ortholog_id = targ_prot_ortholog.split("|")[1]
			
			ortholog_id = ortholog_id + '_' + str(start_prot_ind+1) + '_' + str(end_prot_ind)
			if targ_prot not in ortholog_prots_dict:
				ortholog_prots_dict[targ_prot] = {}
			if targ_prot not in ortholog_dna_dict:
				ortholog_dna_dict[targ_prot] = {}

			ortholog_prots_dict[targ_prot][org_name] = [ortholog_id, prot_seq]
			ortholog_dna_dict[targ_prot][org_name] = [ortholog_id, cds_seq]

	return ortholog_prots_dict, ortholog_dna_dict #human_targ_prots_cds_dict, ortholog_prots_dict, ortholog_dna_dict

def output_fasta(output_prots_fasta_folder, output_dna_fasta_folder, human_targ_prots_cds_dict, ortholog_prots_dict, ortholog_dna_dict):
	""" Output each target protein with its ortholog proteins in the same fasta (do the same for their cds sequences)
	Args:
		output_prots_fasta_folder: output folder to store protein fastas
		output_dna_fasta_folder: output folder to store dna fastas
		
		human_targ_prots_cds_dict: stores the human target protein ids (keys) and lists containing the protein seq its corresponding cds seq (i.e. [prot_seq, cds_seq])
		ortholog_prot_dict: stores for each target prot id (outer key), a nested dict with the organism name as key and ortholog prot seqid and prot seq as lists (i.e. targ_prot: {org_name: [ortholog_prot_id, ortholog_prot_seq]})
		ortholog_dna_dict: stores for each target prot id (outer key), a nested dict with the organism name as key and ortholog prot seqid and cds seq as lists (i.e. targ_prot: {org_name: [ortholog_prot_id, ortholog_dna_cds_seq]})



	"""
	print("### Outputting protein and cds fasta files for human targets + orthologs ###")
	for targ_prot in human_targ_prots_cds_dict: #targ_protid: [targ_prot_seq, targ_prot_corresponding_cds_seq]

		output_prot_fi = osp.join(output_prots_fasta_folder, targ_prot + '.aa.fasta')
		output_dna_fi = osp.join(output_dna_fasta_folder, targ_prot + '.dna.fasta')
		prot_seq = human_targ_prots_cds_dict[targ_prot][0]
		cds_seq = human_targ_prots_cds_dict[targ_prot][1]
		
		with open(output_prot_fi, 'w') as op, open(output_dna_fi, 'w') as od:
			op.write('>' + targ_prot + "\n")
			op.write(prot_seq + "\n")

			od.write('>' + targ_prot +"\n")
			od.write(cds_seq + "\n")

			for org in ortholog_prots_dict[targ_prot]:
				ortholog_id = ortholog_prots_dict[targ_prot][org][0]
				ortholog_prot_seq = ortholog_prots_dict[targ_prot][org][1]
				op.write('>' + ortholog_id + '_' + org  + "\n")
				op.write(ortholog_prot_seq + "\n")

			for org in ortholog_dna_dict[targ_prot]:
				ortholog_id = ortholog_dna_dict[targ_prot][org][0]
				ortholog_dna_seq = ortholog_dna_dict[targ_prot][org][1]
				od.write('>' + ortholog_id + '_' + org  + "\n")
				od.write(ortholog_dna_seq + "\n")




def main():
	script_dir = osp.dirname(__file__)
	human_blast_results_folder = osp.join(script_dir, '..', 'data', 'blast', 'results', 'human_protein_to_cds')
	refseq_cds_folder = osp.join(script_dir, '..', 'data', 'dnds', 'human_targ_prot_cds_fasta_refseq')

	blast_results_folder = osp.join(script_dir, '..', 'data', 'blast', 'results', 'protein_to_cds')
	species_orthologs_prots_dir = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'species_ortholog_prots_fasta')
	species_cds_dir = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'cds')
	human_target_proteins_fasta = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'human_target_proteins.fasta')
	
	output_prots_fasta_folder = osp.join(script_dir, '..', 'data', 'dnds', 'aa_fasta')
	output_dna_fasta_folder = osp.join(script_dir, '..', 'data', 'dnds', 'dna_fasta')

	if not osp.exists(output_prots_fasta_folder):
		os.mkdir(output_prots_fasta_folder)

	if not osp.exists(output_dna_fasta_folder):
		os.mkdir(output_dna_fasta_folder)

	human_targ_prots_cds_dict = store_human_prots_and_cds_into_dict(human_blast_results_folder, blast_results_folder, human_target_proteins_fasta, refseq_cds_folder, species_cds_dir)
	ortholog_prots_dict, ortholog_dna_dict = store_prots_and_cds_into_dicts(blast_results_folder, species_orthologs_prots_dir, species_cds_dir)
	output_fasta(output_prots_fasta_folder, output_dna_fasta_folder, human_targ_prots_cds_dict, ortholog_prots_dict, ortholog_dna_dict)
			

if __name__ == '__main__':
	main()