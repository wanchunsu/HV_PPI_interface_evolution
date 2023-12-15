'''

Process dN/dS results and keep track of special cases
where dS == 0.
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp
import json
from Bio import SeqIO
import csv
from Bio import AlignIO
from fasta_tools import parse_fasta_file_into_dict
from parse_tblastn_results import parse_tblastn_results, parse_and_check_human_tblastn_results
import argparse

'''
1. Map dN and dS values (parsed_dnds files) back onto the amino acid sequence (the translated version, not necessarily the same as the actual protein seq.)
2. Map the dN and dS positions on the amino acid seq back onto the actual protein (use alignment files -- either ensembl/refseq).
4. Calculate dN/dS -- alpha is dS, beta is dN. Take note of when we have special cases (i.e. dS==0)
	2 special cases (don't store dn/ds, just keep track of them -- plot fractions)
		1) only one non-gap residue in amino acid alignment
		2) fully conserved so dn=ds=0
3. Categorize dN/dS values into residue categories (done in next script)
'''

	

def get_translated_aa_msa_alignment_fasta(fasta_file):
	""" Get sequence of translated amino acid seq alignment (including gaps) -- this is the first seq in the msa file

	Args:
		fasta_file: fasta file with msa of translated amino acid and its orthologs
	
	Return:
		first sequence in given fasta file
	"""
	fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
	
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		break
	return sequence


def parse_dnds_positions(dnds_file, msa_aligned_aa_seq):
	""" Parse dN/dS values into a dictionary where keys are positions in the msa and values are the corresponding dN/dS values
	Args:
		dnds_file: file containing processed dN/dS output
		msa_aligned_aa_seq: sequence of translated amino acid seq alignment (including gaps) -- positions correspond exactly to the dnds positions

	Return: dicitonary storing positions on the msa and corresponding dN/dS values and residues 
	"""
	dict_of_dnds_pos = {}
	with open(dnds_file) as dnds_csv:
		csv_reader = csv.reader(dnds_csv)
		header = []
		header = next(csv_reader)
		#header: ['Site', 'alpha', 'beta', 'alpha=beta', 'LRT', 'p-value', 'Total_branch_length', 'p-asmp']
		for row, aa_pos in zip(csv_reader, msa_aligned_aa_seq): 
			site_pos = int(row[0])
			ds = float(row[1])
			dn = float(row[2])
			dict_of_dnds_pos[site_pos] = (dn, ds, aa_pos)
	
	return dict_of_dnds_pos


def get_special_cases(aa_alignment, position_to_check):
	""" Get special cases (check positions that have dS=0). They can be one of the following:
		1) No synonymous substitutions 
		2) Only one non-gapped residue at the alignment position
	Args:
		aa_alignment: parsed amino acid msa
		position_to_check: position with dS=0 to check 
	Return:
		"gap" (only one non-gapped residue at this position) or "no_substitutions" (no substitutions at the position)

	"""
	

	list_of_residues_at_pos = aa_alignment[:,position_to_check-1]

	count_non_gap = 0
	for res in list_of_residues_at_pos:
		if res != '-':
			count_non_gap +=1

	if count_non_gap == 1:
		return "gap"
	else:
		return "no_substitutions"





def map_translated_aa_residues_to_msa_aligned_residues(translated_aa_tblastn_alignment_with_orig_prot, dict_of_dnds_pos, aa_fasta_file):
	""" map dnds values to translated cds alignment positions 

	Args: 
		translated_aa_tblastn_alignment_with_orig_prot: translated cds seq aligned with orig prot -- need this to map residues and dnds vals to orig protein
		dict_of_dnds_pos: dictionary containing positions in the msa and their corresponding dN/dS values
		aa_fasta_file: file containing aa msa for the current protein (this will be for finding special cases for alignment positions -- dS==0)
	Return:
		dict_of_translated_cds_pos_and_dnds_vals: dictionary storing amino acid positions on the translated cds and their dN/dS values
		
		i.e.
		index_on_tblastn_translated_alignment (sseq in tblastn results): (dnds value, residue (amino acid/'-''), position in the msa alignment)
	"""

	#print(translated_aa_tblastn_alignment_with_orig_prot)
	dict_of_translated_cds_pos_and_dnds_vals = {}
	ind_translated_aa = 0

	#parse aa alignment file 
	aa_alignment = AlignIO.read(aa_fasta_file, "fasta")

	#Run through each residue position in the msa alignment --corresponds to positions in the dnds results
	# Match the residues to the original translated aa aligned seq from tblastn and store their indices and dnds results
	# Note there may be gaps in the msa that don't exist in the original translated aa aligned seq, so we skip those since they won't help with mapping back to the original targ prot seq
	for msa_alignment_pos in dict_of_dnds_pos:
		if ind_translated_aa < len(translated_aa_tblastn_alignment_with_orig_prot): #if we're past the last residue in the tblastn alignment, then that means everything else is a gap in the msa
			(dn, ds, msa_alignment_residue) = dict_of_dnds_pos[msa_alignment_pos]
			if translated_aa_tblastn_alignment_with_orig_prot[ind_translated_aa] == msa_alignment_residue: #residue in translated aa seq alignment corresponds to same residue in msa alignment pos
				
				if ds == 0: #special case where ds ==0 (either no syn substitutions at all, or only one non-gapped residue at position)
					case = get_special_cases(aa_alignment, msa_alignment_pos)
					if case == 'gap':
						print(f'{case}: {msa_alignment_pos}')
					dict_of_translated_cds_pos_and_dnds_vals[ind_translated_aa+1] = (0, case, msa_alignment_residue) # store (assign dn/ds=0)
				else: #dS != 0 (well-defined dN/dS)
					dNdS = dn/ds
					dict_of_translated_cds_pos_and_dnds_vals[ind_translated_aa+1] = (dNdS, "normal", msa_alignment_residue) # store 

				# if msa_alignment_residue == '-': #these are the positions that were already gaps on the tblastn alignment (need to keep these b/c correspond to residues on the alignment of the orig targ protein)
				# 	print((ind_translated_aa+ 1, dnds_val, msa_alignment_residue, msa_alignment_pos))
				ind_translated_aa += 1
				
			else: #this is a gap (doesn't correspond to the human prot)
				continue

	return dict_of_translated_cds_pos_and_dnds_vals

def parse_ensembl_tblastn_results(alignment_file):
	""" Parse blast results to store human target prot and its translated cds alignment

	Args:
		alignment_file: tblastn alignment file

	Return: 
		dict_of_alignments: dictionary storing blast results between targ prots and translated cds amino acids (storing start_pos_prot, end_pos_prot, start_pos_cds, end_pos_cds, qseq_aligned, sseq_translated_aligned)

	"""
	dict_of_alignments = {}
	with open(alignment_file) as af:
		for line in af:
			line = line.strip("\n")
			info = line.split("\t")
			qseqid = info[0]

			start_pos_prot = int(info[6])
			end_pos_prot = int(info[7])

			start_pos_cds = int(info[8])
			end_pos_cds = int(info[9])

			qseq_aligned = info[12] #original protein aligned
			sseq_translated_aligned = info[13] #translated aa aligned
			pident = float(info[2])

			if qseqid not in dict_of_alignments: #store only the first one
				# if '-' in sseq_translated_aligned:
				# 	print(line)
				dict_of_alignments[qseqid] = [start_pos_prot, end_pos_prot, start_pos_cds, end_pos_cds, qseq_aligned, sseq_translated_aligned, pident]
	
	return dict_of_alignments


def parse_human_tblastn_results(human_blast_results_folder, blast_results_folder, human_target_proteins):
	""" Parse human against refseq cds seq tblastn results. 
	If identity is <100%, then we check if human targ against ensembl cds seqs give better results.
	If identity is == 100%, we also check if alignment length is same as the protein itself, if not, we check if ensembl cds has longer alignment (given that it also have 100% identity)
	Args:
		blast_results_folder: folder containing all tblastn results
		human_target_proteins_fasta: fasta file containing sequences of our human target proteins
		prot_fasta_folder: folder containing 

	"""
	human_against_ensembl_cds_blast_results = parse_ensembl_tblastn_results(osp.join(blast_results_folder, 'human.prot_against_cds.blast.out'))

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

					if human_against_ensembl_cds_blast_results[qseqid][6] > pident:  #if yes, then choose ensembl cds
						#print(f'\tensembl gives better coverage {pident}, {human_against_ensembl_cds_blast_results[qseqid][7]}')
						dict_of_human_results[qseqid] = human_against_ensembl_cds_blast_results[qseqid]
					else: #if not better, then stick with refseq cds
						dict_of_human_results[qseqid] = [start_pos_prot, end_pos_prot, start_pos_cds, end_pos_cds, qseq_aligned, sseq_translated_aligned] #start_pos_cds-1,end_pos_cds are the indices that will give us the substring at start_pos_cds (inclusive-->end_pos_cds substring (inclusive )
				
				else: #100% identity btwn protein seq and refseq cds
					#check if the refseq cds alignment length is the same as the actual protein
					if query_aligned_len != len(human_target_proteins[qseqid]):
						#print(f'{qseqid} not same length')
						ensembl_alignment_len = human_against_ensembl_cds_blast_results[qseqid][2] - human_against_ensembl_cds_blast_results[qseqid][1]
						ensembl_alignment_pident = human_against_ensembl_cds_blast_results[qseqid][6]
						
						if ensembl_alignment_len > query_aligned_len and ensembl_alignment_pident == 100.0: #check that ensembl alignment gives same 100% identity and better coverage
							#print(f'\tensembl gives longer len: {len(human_target_proteins[qseqid])}, {ensembl_alignment_len}, {query_aligned_len}')
							dict_of_human_results[qseqid] = human_against_ensembl_cds_blast_results[qseqid]
						
						else:
							dict_of_human_results[qseqid] = [start_pos_prot, end_pos_prot, start_pos_cds, end_pos_cds, qseq_aligned, sseq_translated_aligned] #start_pos_cds-1,end_pos_cds are the indices that will give us the substring at start_pos_cds (inclusive-->end_pos_cds substring (inclusive )	
					
					else:
						dict_of_human_results[qseqid] = [start_pos_prot, end_pos_prot, start_pos_cds, end_pos_cds, qseq_aligned, sseq_translated_aligned] #start_pos_cds-1,end_pos_cds are the indices that will give us the substring at start_pos_cds (inclusive-->end_pos_cds substring (inclusive )
	return dict_of_human_results	 


def map_dnds_to_targ_prot_residue_positions(dict_of_alignments, aa_msa_fasta_folder, dnds_folder):
	''' Map dn/ds values (originally in terms of msa positions) back onto original target seq
	Will need to do mappings from msa positions to tblastn subject positions (i.e. translated cds aa seq) then to the positions on the original target protein
	For each targ protein (outer key) store the aligned position indices on the original targ prot seq (key) and dN/dS value (value)

	Args:
		dict_of_alignments: dictionary storing tblastn alignments
		aa_fasta_folder: folder storing amino acid msas
		dnds_folder: folder containing dnds results

	Return: 
		dict_of_targ_prot_dnds_values: Dictionary where each targ protein (outer key) stores the aligned position indices on the original targ prot seq (key) and dN/dS value (value)
	'''

	#for each targ protein store the aligned positions and mapped amino acid from the translated cds seq
	print('Mapping dn/ds values onto the original target protein seq . . .')
	dict_of_targ_prot_dnds_values = {}

	for targ_prot in dict_of_alignments:

		print(f'Done: {targ_prot}')
		
		dict_of_targ_prot_dnds_values[targ_prot] = {}


		#we want to keep indices in terms of targ seq
		
		pos = dict_of_alignments[targ_prot][0]

		targ_prot_aligned  = dict_of_alignments[targ_prot][4]
		translated_aa_aligned = dict_of_alignments[targ_prot][5]

		# parse and get necessary mapping dicts btwn translated aa aligned seq and msa and dn, ds positions and values
		aa_fasta_file = osp.join(aa_msa_fasta_folder, targ_prot + '_msa.fasta')
		dnds_file = osp.join(dnds_folder, targ_prot + '.parsed_dnds.csv') # we want the unprocessd dN and dS values (b/c we'll be calculating dN/dS ourselves taking into account special cases where dS==0)

		msa_aligned_aa_seq = get_translated_aa_msa_alignment_fasta(aa_fasta_file) #get msa alignment seq (this is the seq from msa alignment btwn the translated cds aa seq with its orthologs)
		dict_of_dnds_pos = parse_dnds_positions(dnds_file, msa_aligned_aa_seq) #parse dn/ds results for the current targ prot (note positions are in terms of the msa alignment and not the original translated cds aa seq)
		dict_of_translated_cds_pos_and_dnds_vals = map_translated_aa_residues_to_msa_aligned_residues(translated_aa_aligned, dict_of_dnds_pos, aa_fasta_file) # dictionary storing amino acid positions on the translated cds and their dN/dS values


		for targ_prot_res, translated_aa_res, ind_translated_aa in zip(targ_prot_aligned, translated_aa_aligned, dict_of_translated_cds_pos_and_dnds_vals):
			
			(dnds_val, case, msa_alignment_residue) = dict_of_translated_cds_pos_and_dnds_vals[ind_translated_aa]
			if targ_prot_res == '-': #gap in targ prot alignment, so skip b/c no residue to store
				continue
			else:
				#targ prot mapped to a gap in translated aa (store '-') or mapped to an amino acid (store the amino acid letter)
				if msa_alignment_residue != translated_aa_res: #just to verify that things worked
					print('incorrect')
				dict_of_targ_prot_dnds_values[targ_prot][pos] = [dnds_val, case, translated_aa_res]
				pos += 1
	return dict_of_targ_prot_dnds_values

def output_dnds_results_on_targ_prot(dict_of_targ_prot_dnds_values, output_folder):
	for targ_prot in dict_of_targ_prot_dnds_values:
		with open(osp.join(output_folder, targ_prot + '.dnds.json'), 'w') as of:
			json.dump(dict_of_targ_prot_dnds_values[targ_prot], of)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--hyphy_folder')
	parser.add_argument('-o', '--output_folder')
	args = parser.parse_args()

	hyphy_folder = args.hyphy_folder 
	output_f = args.output_folder

	script_dir = osp.dirname(__file__)
	human_blast_results_folder = osp.join(script_dir, '..', 'data', 'blast', 'results', 'human_protein_to_cds')
	blast_results_folder = osp.join(script_dir, '..', 'data', 'blast', 'results', 'protein_to_cds')
	human_target_proteins_fasta = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'human_target_proteins.fasta')

	aa_msa_fasta_folder = osp.join(script_dir, '..', 'data', 'dnds', 'aa_msas')
	dnds_folder = osp.join(script_dir, '..', 'data', 'dnds', hyphy_folder)


	output_folder = osp.join(script_dir, '..', 'data', 'dnds', output_f)
	if not osp.exists(output_folder):
		os.mkdir(output_folder)

	#Load target proteins
	human_target_proteins = parse_fasta_file_into_dict(human_target_proteins_fasta, change_header=False, split_header_by='|', ind_to_keep = 1)

	# Get dictionary of tblastn alignments for human (this includes ensembl and refseq tblastn results)	
	dict_of_alignments = parse_human_tblastn_results(human_blast_results_folder, blast_results_folder, human_target_proteins)

	# Calculate and map dN/dS values onto the target proteins
	dict_of_targ_prot_dnds_values = map_dnds_to_targ_prot_residue_positions(dict_of_alignments, aa_msa_fasta_folder, dnds_folder)

	#Output
	output_dnds_results_on_targ_prot(dict_of_targ_prot_dnds_values, output_folder)
if __name__ == '__main__':
	main()
