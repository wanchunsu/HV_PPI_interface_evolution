'''
Get RBH orthologs from reciprocal BLAST alignments between human and selected species proteomes.
Calculate sequence mismatches between human target proteins and their orthologs.
================================================================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp
from Bio import SeqIO
from fasta_tools import parse_fasta_file_into_dict
import json
import argparse

def get_list_of_target_proteins(infile):
	list_of_target_proteins = []
	with open(infile) as infi:
		next(infi)
		for line in infi:
			line = line.strip("\n")
			targ_prot = line.split("\t")[0]
			list_of_target_proteins.append(targ_prot)
	return list_of_target_proteins


def make_dict_of_blast_best_hits(blast_results):
	""" Parse blast results into dictionary where hits are sorted by descending bit score, then by ascending e-value
	Args:
		blast_results: file containing blast results
	Return:
		dict_of_blast_best_hits (dict)
	"""
	dict_of_blast_best_hits = {}

	with open(blast_results) as b:
		for line in b:
			line = line.strip("\n")
			info = line.split("\t")
			if '|' in info[0]:
				query_id = info[0].split('|')[1]
			else:
				query_id = info[0]
			subject_id = info[1].split('|')[1]
			query_range = (int(info[6]), int(info[7]))
			subject_range = (int(info[8]), int(info[9]))
			e_val = float(info[10])
			score = float(info[11])
			qseq = info[12]
			sseq = info[13]

			if query_id not in dict_of_blast_best_hits:
				dict_of_blast_best_hits[query_id] = []

			dict_of_blast_best_hits[query_id].append([subject_id, query_range, subject_range, e_val, score, qseq, sseq])
	
	for query_id in dict_of_blast_best_hits:
		#sort the hits by bit score first then eval
		sorted_by_bit_score_then_eval = sorted(dict_of_blast_best_hits[query_id], key = lambda x:(-x[4], x[3]))
		dict_of_blast_best_hits[query_id] = sorted_by_bit_score_then_eval
	
	return dict_of_blast_best_hits


def get_reciprocal_best_hit(infile, human_vs_close_species_blast_results, close_species_vs_human_blast_results, out_rbh_file): 
	""" Get reciprocal best hit from human-close_species and close_species-human blast results for the target proteins of interest

	Args:
		infile: file containing all the exogenous data (to extract the list of target proteins of interest from)
		human_vs_close_species_blast_results: dictionary of blast results for human against close_species (sorted by descending bit score, then by ascending e-value)
		close_species_vs_human_blast_results: dictionary of blast results for close_species against human (sorted by descending bit score, then by ascending e-value)

	Return: 
		rbh_for_targ_prots: dictionary of reciprocal best hits for human target proteins (these should be the orthologs)
	"""
	rbh_for_targ_prots = {}
	human_vs_close_species_dict = make_dict_of_blast_best_hits(human_vs_close_species_blast_results)
	close_species_vs_human_dict = make_dict_of_blast_best_hits(close_species_vs_human_blast_results)

	list_of_target_proteins = get_list_of_target_proteins(infile)
	
	for t in list_of_target_proteins:
		
		best_hit_human_close_species =  human_vs_close_species_dict[t][0][0]
		
		best_hit_close_species_human = close_species_vs_human_dict[best_hit_human_close_species][0][0]
		
		#Case 1: Human target's top hit in close_species also has human target as top hit (RBH found)
		if t ==best_hit_close_species_human:
			
			rbh_for_targ_prots[t] = human_vs_close_species_dict[t][0]
			#print(f"{t}, Case 1: {rbh_for_targ_prots[t]}")
			continue

		#check for other best hits in human for best_hit_human_close_species
		other_hits_with_same_best_scores = [k for k in range(1, len(close_species_vs_human_dict[best_hit_human_close_species])) if close_species_vs_human_dict[best_hit_human_close_species][k][3] ==close_species_vs_human_dict[best_hit_human_close_species][0][3] and close_species_vs_human_dict[best_hit_human_close_species][k][4] ==close_species_vs_human_dict[best_hit_human_close_species][0][4]]
		

		if len(other_hits_with_same_best_scores)>0:
			rbh_found = False
			for k in other_hits_with_same_best_scores:
				# Case 2: Human target's top hit in close_species also has human target as top hit (just not listed at the top) (RBH found)
				if close_species_vs_human_dict[best_hit_human_close_species][k][0] == t:
					rbh_for_targ_prots[t] = human_vs_close_species_dict[t][0] 
					#print(f"{t}, Case 2: {rbh_for_targ_prots[t]}")
					rbh_found = True
	
					break
			if rbh_found == True:
				continue

		#target protein's best hit in close_species is say prot M but prot M's best hit is not target protein
		#Check if there is more than one of the same scoring hit in human-vs-close_species
		#get list of these same scoring hits close_species and if any have top hit in human targ prot, then we assign this as rbh
		best_bit_score = human_vs_close_species_dict[t][0][4]
		best_eval = human_vs_close_species_dict[t][0][3]
		other_highest_score_and_e_val_indices = []
		
		for i in range(1, (len(human_vs_close_species_dict[t]))):
			
			if human_vs_close_species_dict[t][i][4] == best_bit_score and human_vs_close_species_dict[t][i][3] == best_eval:
				other_highest_score_and_e_val_indices.append(i)

		if len(other_highest_score_and_e_val_indices) > 0: #if more than one has the best score and e-val, check if those have best hit as targ protein
			
			
			found_rbh_ind = -1
			for i in other_highest_score_and_e_val_indices:
				other_best_hit = human_vs_close_species_dict[t][i][0]
				best_score_in_human = close_species_vs_human_dict[other_best_hit][0][4]
				best_eval_in_human = close_species_vs_human_dict[other_best_hit][0][3]

				if close_species_vs_human_dict[other_best_hit][0][0] == t: 
					found_rbh_ind = i
					break
				#check if this other best hit in close_species has more than one best hit in human
				other_indices_with_highest_score_in_human = [j for j in range(1, len(close_species_vs_human_dict[other_best_hit])) if close_species_vs_human_dict[other_best_hit][j][4]==best_score_in_human and close_species_vs_human_dict[other_best_hit][j][3]==best_eval_in_human]
				if len(other_highest_score_and_e_val_indices) > 0:
					
					for j in other_indices_with_highest_score_in_human:
						if close_species_vs_human_dict[other_best_hit][j][0] == t:
							
							found_rbh_ind = i
							break
				if found_rbh_ind != -1:
					break

			if found_rbh_ind != -1: #Case 3: Human target's other best hits has human target RBH. (RBH found)
				rbh_for_targ_prots[t] = human_vs_close_species_dict[t][found_rbh_ind]
				#print(f"{t}, Case 3: {rbh_for_targ_prots[t]}")

				continue
			else: # Case 4: no rbh found even after searching the other best hit(s) in close_species and their best hit(s) in human; just take the first hit
				#print(t)
				rbh_for_targ_prots[t] = human_vs_close_species_dict[t][0]
				#print(f"{t}, Case 4: {rbh_for_targ_prots[t]}")

				continue
				
			
		
		else: #Case 5: only one close_species hit with best score and e-val, but not rbh; will take this first hit(even tho not rbh)
			rbh_for_targ_prots[t] = human_vs_close_species_dict[t][0]
			#print(f"{t}, Case 5: {rbh_for_targ_prots[t]}")
			continue

	with open(out_rbh_file, "w") as outfile:
		json.dump(rbh_for_targ_prots, outfile)	

	return rbh_for_targ_prots


def parse_alignment_for_reciprocal_best_hit(rbh_for_targ_prots, alignment_matches_out_fi): 

	
	""" Parses the alignments between proteins and their reciprocal best hit ortholog (matches positions)
	Args:
		rbh_for_targ_prots: rbh orthologs for each protein
		alignment_matches_out_fi: output file for the corresponding alignment positions 
	Return:
		alignment_matches_btwn_human_and_close_species_orthologs: dictionary of alignments between target protein and its ortholog 
	"""

	###
	alignment_matches_btwn_human_and_close_species_orthologs = {}
	#[subject_id, query_range, subject_range, e_val, score, qseq, sseq]
	#go through the aligned qseq and make a dict of the matches/mismatches between the query and subject at each position
	for targ in rbh_for_targ_prots:
		alignment_matches_btwn_human_and_close_species_orthologs[targ] = {}

		subject_id = rbh_for_targ_prots[targ][0]
		query_start = rbh_for_targ_prots[targ][1][0]
		subject_range = rbh_for_targ_prots[targ][2][0]
		qseq = rbh_for_targ_prots[targ][5]
		sseq = rbh_for_targ_prots[targ][6]
		curr_pos = query_start
		for h, m in zip(qseq, sseq):
			if h != '-': #human residue is not a gap
				if h == m:
					alignment_matches_btwn_human_and_close_species_orthologs[targ][curr_pos] = 1
					curr_pos +=1
				else:
					if m != '-': #mismatch (substitution)
						alignment_matches_btwn_human_and_close_species_orthologs[targ][curr_pos] = 0
					else: # gap (targ prot residue maps to a gap in the ortholog)
						alignment_matches_btwn_human_and_close_species_orthologs[targ][curr_pos] = -1
					curr_pos += 1
			else: # case where target seq residue is a gap (nothing to store since there is no targ seq residue at this position)
				continue
	with open(alignment_matches_out_fi, "w") as outfile:
		json.dump(alignment_matches_btwn_human_and_close_species_orthologs, outfile)	
	return alignment_matches_btwn_human_and_close_species_orthologs

def output_reciprocal_best_hit(rbh_for_targ_prots, rbh_with_sseq_out_fi): 

	
	""" Outputs file with reciprocal best hit ortholog (and its position and alignment) for each targ prot
	Args:
		rbh_for_targ_prots: rbh orthologs for each protein
		rbh_with_sseq_out_fi: output file with each prots rbh and the rbh's aligned seq 
	
	"""

	###
	rbh_with_sseq = {}
	#[subject_id, query_range, subject_range, e_val, score, qseq, sseq]
	#go through the aligned qseq and make a dict of the matches/mismatches between the query and subject at each position
	for targ in rbh_for_targ_prots:
		rbh_with_sseq[targ] = []

		subject_id = rbh_for_targ_prots[targ][0]
		query_start = rbh_for_targ_prots[targ][1][0]
		subject_range = rbh_for_targ_prots[targ][2][0]
		qseq = rbh_for_targ_prots[targ][5]
		sseq = rbh_for_targ_prots[targ][6]
		curr_pos = query_start
		

		rbh_with_sseq[targ]= [subject_id, rbh_for_targ_prots[targ][2], sseq]
	with open(rbh_with_sseq_out_fi, "w") as outfile:
		json.dump(rbh_with_sseq, outfile)	
	

def get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, residue_results):
	""" Get the fraction divergence between the target protein and it ortholog (store in dict)
		
	Args:
		alignment_matches_btwn_human_and_close_species_orthologs: dictionary of alignments between target protein and its ortholog 
		residue_results: file of residues (either exo-specific, endo-specific, mimicry, surface, or buried)
	Returns:
		divergence_dict: fraction divergence per target protein

	"""
	divergence_dict = {}
	with open(residue_results) as int_res:
		next(int_res)
		for line in int_res:
			line = line.strip("\n")
			targ_prot = line.split("\t")[0]
			if line.split("\t")[1] == '': #this could happen for the surface residues (some proteins don't have surface residues, so will be an empty string '')
				continue
			int_res = line.split("\t")[1].split(",")
			int_res = [int(r) for r in int_res]
			total_num_int_res = len(int_res)
			total_mismatches = 0
			for i in int_res:

				if i in alignment_matches_btwn_human_and_close_species_orthologs[targ_prot]:
					if alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == 0: #mismatch
						total_mismatches += 1
					elif alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == -1: #maps to a gap
						total_mismatches += 1
				else: # targ prot res not in alignment (i.e. doesn't map to anything on the ortholog -- treat as a gap too)
					total_mismatches += 1

			fraction_divergence = total_mismatches / total_num_int_res 
			divergence_dict[targ_prot] = fraction_divergence

	return divergence_dict

def get_divergence_all(alignment_matches_btwn_human_and_close_species_orthologs, residue_results):
	""" Get the fraction divergence between the target protein and its rbh ortholog for the given type of residue
	Args:
		alignment_matches_btwn_human_and_close_species_orthologs: dictionary of alignments between target protein and its ortholog 
		residue_results: file of residues (either exo-specific, endo-specific, mimicry, surface, or buried)
	Returns:
		fraction of divergence across all the residues in the specified category (rounded to 5 decimals)

	"""
	divergence_list = []
	total_num_int_res = 0
	total_mismatches = 0

	# separating mismatches into gaps and differences 
	total_gaps = 0
	total_differences = 0
	list_of_gaps_and_non_gaps = []
	list_of_subs_and_non_subs = []

	with open(residue_results) as int_res:
		next(int_res)
		for line in int_res:
			line = line.strip("\n")

			targ_prot = line.split("\t")[0]
			if line.split("\t")[1] == '': #this could happen for the surface residues (some proteins don't have surface residues, so will be an empty string '')
				continue
			int_res = line.split("\t")[1].split(",")
			int_res = [int(r) for r in int_res]
			total_num_int_res += len(int_res)
			
			for i in int_res:

				if i in alignment_matches_btwn_human_and_close_species_orthologs[targ_prot]:
					if alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == 0: #maps to a mismatch
						total_mismatches += 1
						total_differences += 1
					elif alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == -1: #maps to a gap
						total_mismatches += 1
						total_gaps += 1
				else: # targ prot res not in alignment (i.e. doesn't map to anaything on the ortholog -- treat as a gap too)
					total_mismatches += 1
					total_gaps += 1

	fraction_divergence = total_mismatches/total_num_int_res

	fraction_div_that_are_gaps = total_gaps/total_mismatches
	fraction_div_that_are_differences = total_differences/total_mismatches
	
	return round(fraction_divergence*100, 3), round(fraction_div_that_are_gaps*100, 1), round(fraction_div_that_are_differences*100, 1), total_num_int_res #return percent divergence, percent gaps, percent differences (note: gaps + differences = divergence)



def get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, residue_results):
	""" Get list of 0 and 1s representing substitutions for residues in a residue category, or representing gaps for residues in a residue category
		Used for boostrap method for calculating standard error for average divergence values for each residue type
	"""
	total_num_int_res = 0
	
	list_of_gaps_and_non_gaps = [] #store 1 if gap else store 0
	list_of_subs_and_non_subs = [] #store 1 if substituion else store 0
	list_of_div_and_non_divs = [] #store 1 if is a difference (either gap/substitution) else store 0

	with open(residue_results) as int_res:
		next(int_res)
		for line in int_res:
			line = line.strip("\n")

			targ_prot = line.split("\t")[0]
			if line.split("\t")[1] == '': #this could happen for the surface residues (some proteins don't have surface residues, so will be an empty string '')
				continue
			int_res = line.split("\t")[1].split(",")
			int_res = [int(r) for r in int_res]
			total_num_int_res += len(int_res)
			
			for i in int_res:

				if i in alignment_matches_btwn_human_and_close_species_orthologs[targ_prot]:
					if alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == 0: #maps to a mismatch
						
						list_of_subs_and_non_subs.append(1) #substitution so add a 1 to substitution list
						list_of_gaps_and_non_gaps.append(0) # not a gap so add a 0 to gap list
						list_of_div_and_non_divs.append(1) #is a difference so add 1 to divergence list

					elif alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == -1: #maps to a gap
						list_of_subs_and_non_subs.append(0) # not substitution so add a 0 to substitution list
						list_of_gaps_and_non_gaps.append(1) # gap so add a 1 to gap list
						list_of_div_and_non_divs.append(1) #is a difference so add 1 to divergence list
					
					else: #match, neither gap nor substitution, so add 0 to both lists
						list_of_subs_and_non_subs.append(0) # not a substitution so add a 0 to substitution list
						list_of_gaps_and_non_gaps.append(0) # not gap so add a 0 to gap list
						list_of_div_and_non_divs.append(0) # not a difference so add 0 to divergence list

				else: # targ prot res not in alignment (i.e. doesn't map to anaything on the ortholog -- treat as a gap too)
					list_of_subs_and_non_subs.append(0) # not substitution so add a 0 to substitution list
					list_of_gaps_and_non_gaps.append(1) # gap so add a 1 to gap list
					list_of_div_and_non_divs.append(1) #is a difference so add 1 to divergence list

	
	return  list_of_div_and_non_divs, list_of_subs_and_non_subs, list_of_gaps_and_non_gaps


def get_substitutions_and_gaps_all(alignment_matches_btwn_human_and_close_species_orthologs, residue_results):
	""" Get the fraction divergence (separating actual mismatches --substitutions and gaps) between the target protein and its rbh ortholog for the given type of residue
	Args:
		alignment_matches_btwn_human_and_close_species_orthologs: dictionary of alignments between target protein and its ortholog 
		residue_results: file of residues (either exo-specific, endo-specific, mimicry, surface, or buried)
	Returns:
		fraction of divergence for substitutions and gaps

	"""
	divergence_list = []
	total_num_int_res = 0
	total_mismatches = 0

	# separating mismatches into gaps and substitutions
	
	total_substitutions = 0
	total_gaps = 0

	with open(residue_results) as int_res:
		next(int_res)
		for line in int_res:
			line = line.strip("\n")

			targ_prot = line.split("\t")[0]
			if line.split("\t")[1] == '': #this could happen for the surface residues (some proteins don't have surface residues, so will be an empty string '')
				continue
			int_res = line.split("\t")[1].split(",")
			int_res = [int(r) for r in int_res]
			total_num_int_res += len(int_res)
			
			for i in int_res:

				if i in alignment_matches_btwn_human_and_close_species_orthologs[targ_prot]:
					if alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == 0: #real mismatch (substitution)
						total_mismatches += 1
						total_substitutions += 1
					elif alignment_matches_btwn_human_and_close_species_orthologs[targ_prot][i] == -1: #maps to a gap
						total_mismatches += 1
						total_gaps += 1	
				else: # targ prot res not in alignment (i.e. doesn't map to anaything on the ortholog -- treat as a gap too)
					total_mismatches += 1
					total_gaps += 1

	#fraction_divergence = total_mismatches/total_num_int_res

	
	fraction_substitutions = total_substitutions/total_num_int_res
	fraction_gaps = total_gaps/total_num_int_res

	return round(fraction_substitutions*100, 3), round(fraction_gaps*100, 3)


def get_avg_divergence(alignment_matches_btwn_human_and_close_species_orthologs, residue_results):
	""" Get the avg fraction divergence between the target protein and it ortholog (averaging across the target proteins)
		(i.e. we have one fraction divergence value per target protein; take avg across these values)
	Args:
		alignment_matches_btwn_human_and_close_species_orthologs: dictionary of alignments between target protein and its ortholog 
		residue_results: file of residues (either exo-specific, endo-specific, mimicry, surface, buried)
	
	Returns: 
		average divergence across target proteins (rounded to 5 decimals)
	"""
	divergence_dict = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, residue_results)

	sum_of_divergence = sum(divergence_dict.values())

	num_vals = len(divergence_dict)
	
	avg_divergence = sum_of_divergence/num_vals

	return round(avg_divergence*100, 3) #return percent divergence





def main():

	script_dir = osp.dirname(__file__)
	
	# Get command line arguments for specific folder and file
	
	parser = argparse.ArgumentParser()

	parser.add_argument('-f', '--folder')
	
	args = parser.parse_args()
	folder_name = args.folder #'v_target_h' 
	
	all_exo = osp.join(script_dir, '..', 'data', 'categorized_interfacial_residues', folder_name, 'exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv')
	
	# Files with exo-specific, endo-specific, and mimicry interfacial residues as well as surface and buried residues
	exo_specific_file =  osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_specific_interfacial_residues.tsv")
	endo_specific_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_specific_interfacial_residues.tsv")
	mimicry_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "mimicked_interfacial_residues.tsv")
	surface_res_file = osp.join(script_dir, '..', 'data', "surface_residues", folder_name, "valid_surface_residues_on_uniprot_seq.tsv")
	buried_res_file = osp.join(script_dir, '..', 'data', "buried_residues", folder_name, "valid_buried_residues_on_uniprot_seq.tsv")

	# Files for all exo and all endo
	all_exo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")
	all_endo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")

	output_fraction_divergence_json = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', folder_name, 'fraction_divergence_with_close_species.json')
	
	if not osp.exists(osp.dirname(output_fraction_divergence_json)):
		os.mkdir(osp.dirname(output_fraction_divergence_json))

	list_of_divergence_values = {}
	per_protein_divergence = {}

	print(f"###   Divergence values for {folder_name}: ###  \n")
	for blast_res in os.listdir(osp.join(script_dir, '..', 'data', 'blast', 'results', 'reciprocal_best_hit')):
		
		if blast_res.endswith("against_human.blast.out"): 
			
			species_name = blast_res.split("_")[0]
			print(f"### Calculating divergence of interfacial and non-interfacial residues ({species_name} comparison) ###")

			# human-to-close_species and close_species-to-human blast results
			human_vs_close_species_blast_results = osp.join(script_dir, '..', 'data', 'blast', 'results', 'reciprocal_best_hit', 'human_against_' + species_name +'.blast.out')
			close_species_vs_human_blast_results = osp.join(script_dir, '..', 'data', 'blast', 'results', 'reciprocal_best_hit', blast_res)

			orthologs_folder = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', folder_name, 'rbh', species_name)

			if not osp.exists(orthologs_folder):
				os.makedirs(orthologs_folder)


			#output files for verification
			out_rbh_file = osp.join(orthologs_folder, "rbh_orthologs.json")
			alignment_matches_out_fi = osp.join(orthologs_folder, "rbh_ortholog_alignment_matches.json")
			rbh_with_sseq_out_fi = osp.join(orthologs_folder, "rbh_orthologs_with_aligned_pos_and_seq.json")

			# get reciprocal best hit orthologs for each target protein
			rbh_for_targ_prots = get_reciprocal_best_hit(all_exo, human_vs_close_species_blast_results, close_species_vs_human_blast_results, out_rbh_file)

			#Get alignment between the human and close_species RBH orthologs
			alignment_matches_btwn_human_and_close_species_orthologs = parse_alignment_for_reciprocal_best_hit(rbh_for_targ_prots, alignment_matches_out_fi)

			#output the reciprocal best hit ortholog for each targ prot (as well as the rbh's aligned pos'n range and its aligned seq)
			output_reciprocal_best_hit(rbh_for_targ_prots, rbh_with_sseq_out_fi)

			############## Calculating divergence for each protein individually (each protein will have its own divergence value) ##############

			#Output files for fraction of divergence of each target protein's exo-specific, endo-specific, and mimicry-specific interfacial residues and surface residues & buried residues
			
			per_prot_divergence_out = osp.join(orthologs_folder, "frac_div_per_protein.json")

			exo_per_prot_divergence = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, exo_specific_file)
			endo_per_prot_divergence = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, endo_specific_file)
			mimicry_per_prot_divergence = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, mimicry_file)
			surface_res_per_prot_divergence = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, surface_res_file)
			buried_res_per_prot_divergence = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, buried_res_file)
			
			all_exo_per_prot_divergence = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, all_exo_file)
			all_endo_per_prot_divergence = get_divergence(alignment_matches_btwn_human_and_close_species_orthologs, all_endo_file)

			per_protein_divergence = {'exo': exo_per_prot_divergence, 'endo': endo_per_prot_divergence, 'mimicry': mimicry_per_prot_divergence,
			 'surface_res': surface_res_per_prot_divergence, 'buried_res': buried_res_per_prot_divergence, 
			 'all_exo': all_exo_per_prot_divergence, 'all_endo': all_endo_per_prot_divergence}
			
			with open(per_prot_divergence_out, 'w') as out_frac_div_per_prot:
				json.dump(per_protein_divergence, out_frac_div_per_prot)
			#Output the divergence on each target proteins exo-specific, endo-specific, and mimicry interfacial residues (and surface residues & buried residues)
			
			#initialize with entries we want to store
			list_of_divergence_values[species_name] = {'avg':{}, 'all':{},'all_percent_gaps': {}, 'total_num_res': {}}

			############## Calculating divergence based on average of divergence across proteins ##############

			# Get the avg fraction of divergence of interfacial and non-interfacial when comparing with close_species ortholog  (avg of divergence across proteins)
			exo_avg_divergence = get_avg_divergence(alignment_matches_btwn_human_and_close_species_orthologs, exo_specific_file)
			endo_avg_divergence = get_avg_divergence(alignment_matches_btwn_human_and_close_species_orthologs, endo_specific_file)
			mimicry_avg_divergence = get_avg_divergence(alignment_matches_btwn_human_and_close_species_orthologs, mimicry_file)
			surface_res_avg_divergence = get_avg_divergence(alignment_matches_btwn_human_and_close_species_orthologs, surface_res_file)
			buried_res_avg_divergence = get_avg_divergence(alignment_matches_btwn_human_and_close_species_orthologs, buried_res_file)

			all_exo_avg_divergence = get_avg_divergence(alignment_matches_btwn_human_and_close_species_orthologs, all_exo_file)
			all_endo_avg_divergence = get_avg_divergence(alignment_matches_btwn_human_and_close_species_orthologs, all_endo_file)

			# print(f"Avg divergence of exogenous specific interfacial residues: {exo_avg_divergence}")
			# print(f"Avg divergence of endogenous specific interfacial residues: {endo_avg_divergence}")
			# print(f"Avg divergence of mimicry specific interfacial residues: {mimicry_avg_divergence}")
			# print(f"Avg divergence of surface residues: {surface_res_avg_divergence}")
			# print(f"Avg divergence of buried residues: {buried_res_avg_divergence}\n")
			
			list_of_divergence_values[species_name]['avg'] = {"exo": exo_avg_divergence, "endo": endo_avg_divergence, "mimicry": mimicry_avg_divergence, 
			"surface_res": surface_res_avg_divergence, "buried_res": buried_res_avg_divergence, 
			"all_exo": all_exo_avg_divergence, "all_endo": all_endo_avg_divergence}
			
			

			############## Calculating divergence across all interfacial residues and surface residues & buried residues##############
			
			# Get the fraction of divergence of interfacial residues and non-interfacial residues when comparing with close_species ortholog  (divergence across all int res (regardless of which protein))
			# Output these divergence values along with the percentages of substitutions and gaps 
			exo_all_divergence, exo_all_gaps, exo_all_differences, exo_total_num_res = get_divergence_all(alignment_matches_btwn_human_and_close_species_orthologs, exo_specific_file)
			endo_all_divergence, endo_all_gaps, endo_all_differences, endo_total_num_res = get_divergence_all(alignment_matches_btwn_human_and_close_species_orthologs, endo_specific_file)
			mimicry_all_divergence, mimicry_all_gaps, mimicry_all_differences, mimicry_total_num_res = get_divergence_all(alignment_matches_btwn_human_and_close_species_orthologs, mimicry_file)
			surface_res_all_divergence, surface_res_all_gaps, surface_res_all_differences, surface_res_total_num_res = get_divergence_all(alignment_matches_btwn_human_and_close_species_orthologs, surface_res_file)
			buried_res_all_divergence, buried_res_all_gaps, buried_res_all_differences, buried_res_total_num_res = get_divergence_all(alignment_matches_btwn_human_and_close_species_orthologs, buried_res_file)

			all_exo_all_divergence, all_exo_all_gaps, all_exo_all_differences, all_exo_total_num_res = get_divergence_all(alignment_matches_btwn_human_and_close_species_orthologs, all_exo_file)
			all_endo_all_divergence, all_endo_all_gaps, all_endo_all_differences, all_endo_total_num_res = get_divergence_all(alignment_matches_btwn_human_and_close_species_orthologs, all_endo_file)

			print(f"Divergence across all exogenous specific interfacial residues: {exo_all_divergence} (G: {exo_all_gaps}%)")
			print(f"Divergence across all endogenous specific interfacial residues: {endo_all_divergence} (G: {endo_all_gaps}%)")
			print(f"Divergence across all mimicry specific interfacial residues: {mimicry_all_divergence} (G: {mimicry_all_gaps}%)")
			print(f"Divergence across all surface residues: {surface_res_all_divergence} (G: {surface_res_all_gaps}%)")
			print(f"Divergence across all buried residues: {buried_res_all_divergence} (G: {buried_res_all_gaps}%)\n")
			
			print(f"Divergence across all exogenous interfacial residues: {all_exo_all_divergence} (G: {all_exo_all_gaps}%)")
			print(f"Divergence across all endogenous interfacial residues: {all_endo_all_divergence} (G: {all_endo_all_gaps}%)\n")

			list_of_divergence_values[species_name]['all'] = {"exo": exo_all_divergence, "endo": endo_all_divergence, 
			"mimicry": mimicry_all_divergence, "surface_res": surface_res_all_divergence, "buried_res": buried_res_all_divergence,
			"all_exo": all_exo_all_divergence, "all_endo": all_endo_all_divergence}

			list_of_divergence_values[species_name]['all_percent_gaps'] = {"exo_percent_gaps_diff": (exo_all_gaps, exo_all_differences), "endo_percent_gaps_diff": (endo_all_gaps, endo_all_differences), 
			"mimicry_percent_gaps_diff": (mimicry_all_gaps, mimicry_all_differences), 
			"surface_res_percent_gaps_diff": (surface_res_all_gaps, surface_res_all_differences),
			"buried_res_percent_gaps_diff": (buried_res_all_gaps, buried_res_all_differences),
			"all_exo_percent_gaps_diff": (all_exo_all_gaps, all_exo_all_differences), 
			"all_endo_percent_gaps_diff": (all_endo_all_gaps, all_endo_all_differences)}

			
			#separating mismatches into substitutions and gaps (indels) and outputting divergence in terms of substitutions and gaps separately
			
			exo_substitutions, exo_gaps =  get_substitutions_and_gaps_all(alignment_matches_btwn_human_and_close_species_orthologs, exo_specific_file)
			endo_substitutions, endo_gaps =  get_substitutions_and_gaps_all(alignment_matches_btwn_human_and_close_species_orthologs, endo_specific_file)
			mimicry_substitutions, mimicry_gaps =  get_substitutions_and_gaps_all(alignment_matches_btwn_human_and_close_species_orthologs, mimicry_file)
			surface_res_substitutions, surface_res_gaps =  get_substitutions_and_gaps_all(alignment_matches_btwn_human_and_close_species_orthologs, surface_res_file)
			buried_res_substitutions, buried_res_gaps =  get_substitutions_and_gaps_all(alignment_matches_btwn_human_and_close_species_orthologs, buried_res_file)

			all_exo_substitutions, all_exo_gaps =  get_substitutions_and_gaps_all(alignment_matches_btwn_human_and_close_species_orthologs, all_exo_file)
			all_endo_substitutions, all_endo_gaps =  get_substitutions_and_gaps_all(alignment_matches_btwn_human_and_close_species_orthologs, all_endo_file)

			list_of_divergence_values[species_name]['all_substitutions'] = {
			"exo": (exo_substitutions), 
			"endo": (endo_substitutions), 
			"mimicry": (mimicry_substitutions), 
			"surface_res": (surface_res_substitutions),
			"buried_res": (buried_res_substitutions),
			"all_exo": (all_exo_substitutions), 
			"all_endo": (all_endo_substitutions)}

			list_of_divergence_values[species_name]['all_gaps'] = {
			"exo": (exo_gaps), 
			"endo": (endo_gaps), 
			"mimicry": (mimicry_gaps), 
			"surface_res": (surface_res_gaps),
			"buried_res": (buried_res_gaps),
			"all_exo": (all_exo_gaps), 
			"all_endo": (all_endo_gaps)}


			print(f"(substitutions, gaps) all exogenous specific interfacial residues: {(exo_substitutions, exo_gaps)})")
			print(f"(substitutions, gaps) all endogenous specific interfacial residues: {(endo_substitutions, endo_gaps)}")
			print(f"(substitutions, gaps) all mimicry specific interfacial residues: {(mimicry_substitutions, mimicry_gaps)}")
			print(f"(substitutions, gaps) all surface residues: {(surface_res_substitutions, surface_res_gaps)}")
			print(f"(substitutions, gaps) all buried residues: {(buried_res_substitutions, buried_res_gaps)}\n")

			print(f"(substitutions, gaps) all exogenous interfacial residues: {(all_exo_substitutions, all_exo_gaps)})")
			print(f"(substitutions, gaps) all endogenous interfacial residues: {(all_endo_substitutions, all_endo_gaps)}\n")


			# Get the total number of residues in each residue category (this will be used for plotting -- i.e. calculating standard err)
			list_of_divergence_values[species_name]['total_num_res'] = {"exo": exo_total_num_res, "endo": endo_total_num_res, 
			"mimicry": mimicry_total_num_res, "surface_res": surface_res_total_num_res, "buried_res": buried_res_total_num_res,
			"all_exo": all_exo_total_num_res, "all_endo": all_endo_total_num_res}


			################## Get list of 0/1s to store whether a residues in a category are non-subs/subs OR non-gap/gap #############
			# This will be for calculating standard errors using bootstrap
			# Get the list of binary 0/1s representing non-substitution/substitution or non-gap/gap at each residue in a residue category
			exo_div_list, exo_subs_list, exo_gaps_list = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, exo_specific_file)
			endo_div_list, endo_subs_list, endo_gaps_list = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, endo_specific_file)
			mimicry_div_list, mimicry_subs_list, mimicry_gaps_list = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, mimicry_file)
			surface_res_div_list, surface_res_subs_list, surface_res_gaps_list = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, surface_res_file)
			buried_res_div_list, buried_res_subs_list, buried_res_gaps_list = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, buried_res_file)

			all_exo_div_list, all_exo_subs_list, all_exo_gaps_list = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, all_exo_file)
			all_endo_div_list, all_endo_subs_list, all_endo_gaps_list = get_subs_and_gaps_binary_label_for_all_residues(alignment_matches_btwn_human_and_close_species_orthologs, all_endo_file)

			list_of_divergence_values[species_name]['div_list'] = {"exo": exo_div_list, "endo": endo_div_list, "mimicry": mimicry_div_list, "surface_res": surface_res_div_list, "buried_res": buried_res_div_list, "all_exo": all_exo_div_list, "all_endo": all_endo_div_list}

			list_of_divergence_values[species_name]['subs_list'] = {"exo": exo_subs_list, "endo": endo_subs_list, "mimicry": mimicry_subs_list, "surface_res": surface_res_subs_list, "buried_res": buried_res_subs_list, "all_exo": all_exo_subs_list, "all_endo": all_endo_subs_list}
	
			list_of_divergence_values[species_name]['gaps_list'] = {"exo": exo_gaps_list, "endo": endo_gaps_list, "mimicry": mimicry_gaps_list, "surface_res": surface_res_gaps_list, "buried_res": buried_res_gaps_list, "all_exo": all_exo_gaps_list, "all_endo": all_endo_gaps_list}
	

	#print(list_of_divergence_values)
	with open(output_fraction_divergence_json, "w") as outfile:
		json.dump(list_of_divergence_values, outfile)	
	

if __name__ == '__main__':
	main()

