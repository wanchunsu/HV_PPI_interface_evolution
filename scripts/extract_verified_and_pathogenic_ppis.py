'''

Extract verified and pathogenic exogenous PPIs into a new exogenous PPI file. 
Also makes a new endogenous file with only target proteins that exist in the 
new exogenous file.
============================================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp
import argparse
from get_interfacial_residues_using_sasa import get_stats_for_interface_residues_exogenous, get_stats_for_interface_residues_endogenous

def load_verification_info(verif_file):
	""" Load file with potential ppi verification annotations
	Args:
		verif_file (tsv file): file containing potential ppis with verificaiton annotations ('T' if this ppi is a verified pathogenic ppi, 'F' otherwise)

	Return:
		dict_of_verifications (dict): dictionary storing potential ppi pairs and their verification annotation
	"""
	dict_of_verifications = {}
	with open(verif_file) as v:
		for line in v:
			line = line.strip("\n")
			info = line.split("\t")

			targ = info[0]
			part = info[1]
			verif_anno = info[2]

			dict_of_verifications[(targ, part)] = verif_anno
	return dict_of_verifications

def check_and_output_known_verified_and_pathogenic_ppis(input_file, dict_of_verifications, output_file):
	''' Go through input_file of ppis and check that the ppi is pathogenic and either already known or verified (based on dict_of_verificaitons). If yes, output into new file, otherwise ignore
	Note: we do not want any immune PPIs, so can ignore those regardless of whether they are verified.
	Args:
		input_file: input file of exogenous ppis to check
		dict_of_verificaitons (dict): dictionary storing potential ppi pairs and their verification annotation
		output_file (tsv file): file to output the known/verified and pathogenic exo ppis

	'''
	with open(input_file) as inp, open(output_file, 'w') as outfi:
		header = next(inp)
		header = header.rsplit('\t', 1)[0] #get header (without the antibody annotation column)
		outfi.write(header + "\n")
		for line in inp:
			line = line.strip("\n")
			info_split_anno = line.rsplit('\t', 1) # this will split up all info from the anno(T or F) 
			# e.g. 'Q9BYF1	P59594	7xo8	F	A	531,534,574,575	F'' ==> ['Q9BYF1	P59594	7xo8	F	A	531,534,574,575', 'F']
			main_info = info_split_anno[0]
			antibody_anno = info_split_anno[1]
			
			if antibody_anno == 'T': # we don't want immune PPIs
				continue
			targ = main_info.split("\t")[0]
			part = main_info.split("\t")[1]
			targ_part_pair = (targ, part)
			if targ_part_pair not in dict_of_verifications: #this is already a known pathogenic PPI
				outfi.write(main_info + "\n") #write the original line (without the anitbody anno)
				#print(targ_part_pair)
			else:
				if dict_of_verifications[targ_part_pair] == 'T':
					outfi.write(main_info + "\n")
				else:# 'F' not a verified pathogenic PPI, so skip
					#print(targ_part_pair)
					continue

def get_only_endo_with_verified_pathogenic_ppis(exo_file, endo_input, endo_output):
	""" Go through endogenous PPI file and remove PPIs where the target protein does not have any pathogenic, verified (or known) exo ppis
	Args:
		exo_file: file containing only pathogenic and verified/known exogenous ppis 
		endo_input: file containing endogenous PPIs (to check whether targ prots have verified and pathogenic ppis)
		endo_output: output file with only endogenous ppis where the target proteins participate in pathogenic and verified exo ppis

	"""
	list_of_exo_targs = []
	list_of_endo_targs_to_delete = []
	with open(exo_file) as ex:
		next(ex)
		for line in ex:
			targ_prot = line.split("\t")[0]
			if targ_prot not in list_of_exo_targs:
				list_of_exo_targs.append(targ_prot)

	with open(endo_input) as endin, open(endo_output, 'w') as endout:
		header = next(endin)
		endout.write(header)
		for line in endin:
			targ_p = line.split("\t")[0]
			if targ_p not in list_of_exo_targs:
				if targ_p not in list_of_endo_targs_to_delete:
					list_of_endo_targs_to_delete.append(targ_p)
			else:
				endout.write(line)

	print(list_of_endo_targs_to_delete)

def main():
	script_dir = osp.dirname(__file__)

	# file with both known and potential PPIs
	input_file = osp.join(script_dir, '..', 'data', 'annotated_combined_homology_and_pdb_anno_human_virus_interface_residues.tsv') 

	# file with pathogenic potential ppis annotated with 'T' (verified) or 'F' (not verified)
	verif_file = osp.join('..', 'data', 'pdb_annotations', 'verify_ppis', 'pathogenic_exogenous_ppis_to_verify.tsv')

	# output file of only known & verified pathogenic PPIs (unverified removed)
	exo_output_file = osp.join(script_dir, '..', 'data', 'exo_verified_pathogenic_int_res.tsv') 

	# input endo file (to check if target prots have pathogenic and verified PPIs)
	endo_input = osp.join(script_dir, '..', 'data', 'endo_structural_annotation', 'interface_residues', 'endo_human_interface_residues.tsv')

	#output endo file with only targ prots that participate in pathogenic and verified exo ppis
	endo_output = osp.join(script_dir, '..', 'data', 'endo_with_verified_exo_int_res.tsv')

	# Load potential ppis annotated with 'T' (verified) or 'F' (not verified)
	dict_of_verifications = load_verification_info(verif_file)
	#print(dict_of_verifications)

	# Get only pathogenic verified (known and verified potential ppis) PPIs and output into a new file
	check_and_output_known_verified_and_pathogenic_ppis(input_file, dict_of_verifications, exo_output_file)

	#Print stats for exo output file
	print('###   Stats for exogenous ppis   ###')
	get_stats_for_interface_residues_exogenous(exo_output_file)

	#Get only endo where targ prots have pathogenic and verified exo ppis
	get_only_endo_with_verified_pathogenic_ppis(exo_output_file, endo_input, endo_output)

	#Print stats for endo output file
	print('\n###   Stats for endogenous ppis   ###')
	get_stats_for_interface_residues_endogenous(endo_output)

if __name__ == '__main__':
	main()