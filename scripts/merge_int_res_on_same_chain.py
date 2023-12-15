'''

Merge interfacial residues on the same chains and 
output relevant fasta files for mapping interfacial
residues on PDB chains back onto human target proteins.
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
import os
from Bio import SeqIO

def format_into_dict(results):
	""" Formats the results into a dictionary to store the interfacial residues on the target protein 
	This will merge together results where the interfacial residues are on the same chain


	Args:
		results (tsv file): file with results (exo or endo)

	Returns:
		dictionary of results formatted as follows:
			{target_prot: 
				{pdbid: 
					{chainid: [], ...
					} 
				...
				}
			...
			}

	"""
	results_dict = {}
	
	with open(results) as res:
		next(res)
		for line in res:
			line = line.strip("\n")
			info = line.split("\t")
			targ_prot = info[0]
			pdb = info[2]
			chainid = info[3]
			pdb_chain = pdb+'_'+ chainid
			int_res = info[5].split(",") #convert to list

			if targ_prot not in results_dict:
				results_dict[targ_prot] = {}
			if pdb_chain not in results_dict[targ_prot]:
				results_dict[targ_prot][pdb_chain] = []
					
			results_dict[targ_prot][pdb_chain].extend(int_res) #merge together with existing residues
			results_dict[targ_prot][pdb_chain] = list(set(results_dict[targ_prot][pdb_chain])) #to delete any repeated residues
	return results_dict

def write_results(results_dict, out_fi, exo_or_endo):
	""" Write results in the results_dict into an output file with the following header
	"""
	with open(out_fi, 'w') as o:
		if exo_or_endo =="exo":
			o.write("target_prot\tpdb_chain\texo_interfacial_residues\n")
		else:
			o.write("target_prot\tpdb_chain\tendo_interfacial_residues\n")
		for targ_prot in results_dict:
			for pdb_chain in results_dict[targ_prot]:
				
				int_res = results_dict[targ_prot][pdb_chain]
				o.write(targ_prot +"\t" + pdb_chain  +"\t" +  ",".join(int_res) + "\n")

def get_lists_of_targ_seqs_and_pdb_ids(exo_dict, endo_dict):
	""" Get lists of target sequences and corresponding pdb ids in exogenous and endogenous dicts
	"""
	list_of_targ_prots = []
	list_of_pdb_chains = []

	for targ_prot in exo_dict:
		if targ_prot not in list_of_targ_prots:
			list_of_targ_prots.append(targ_prot)
		for pdb_chain in exo_dict[targ_prot]:

			if pdb_chain not in list_of_pdb_chains:
				list_of_pdb_chains.append(pdb_chain)
	
	for targ_prot in endo_dict:
		if targ_prot not in list_of_targ_prots:
			list_of_targ_prots.append(targ_prot)
		for pdb_chain in endo_dict[targ_prot]:
			if pdb_chain not in list_of_pdb_chains:
				list_of_pdb_chains.append(pdb_chain)

	return list_of_targ_prots, list_of_pdb_chains

def format_exo_and_endo_into_one_dict(exo_dict, endo_dict):
	"""
	{targ_prot1: [pdb_chain1, pdb_chain2, ...],
	...}
	"""
	exo_and_endo_dict = {}

	for targ_prot in exo_dict:
		if targ_prot not in exo_and_endo_dict:
			exo_and_endo_dict[targ_prot] = []
		for pdb_chain in exo_dict[targ_prot]:
			if pdb_chain not in exo_and_endo_dict[targ_prot]:
				exo_and_endo_dict[targ_prot].append(pdb_chain)

	for targ_prot in endo_dict:
		if targ_prot not in exo_and_endo_dict:
			exo_and_endo_dict[targ_prot] = []
		for pdb_chain in endo_dict[targ_prot]:
			if pdb_chain not in exo_and_endo_dict[targ_prot]:
				exo_and_endo_dict[targ_prot].append(pdb_chain)
	return exo_and_endo_dict


def get_relevant_dict_of_fasta(fasta_iterable, list_to_check):
	dict_of_relevant_fasta = {}
	for fasta in fasta_iterable:
		name, sequence = fasta.id, str(fasta.seq)
		if name in list_to_check and name not in dict_of_relevant_fasta:
			dict_of_relevant_fasta[name] = sequence
	
	#this will also weed out erroneous annotations where the human target protein is not derived from human (e.g. in pdb 1jma (P57083) and 5fgy(P03274)
	return dict_of_relevant_fasta



def get_fasta_seqs_as_individual_files(exo_dict, endo_dict, pdb_fasta, human_fasta, pdb_out_folder, human_out_folder):
	"""
		Get the fasta sequences for each of the elements in the lists of target sequences and corresponding pdb ids and save as individual files.
		These individual fasta files  will be used for mapping the pdb chain residues to its corresponding residue in the uniprot seq .
		The uniprot seq is the gold standard, so we need to find the corresponding interfacial residues that are present on the mapped pdb chains.

	"""
	list_of_targ_prots, list_of_pdb_chains = get_lists_of_targ_seqs_and_pdb_ids(exo_dict, endo_dict) 
	#print(len(list_of_targ_prots))
	#print(len(list_of_pdb_chains))

	dict_of_pdbs_for_each_targ_prot = format_exo_and_endo_into_one_dict(exo_dict, endo_dict) # targ_prot -> list of pdb_chains associated with it
	#print(len(dict_of_pdbs_for_each_targ_prot))

	all_pdb_fasta_sequences = SeqIO.parse(open(pdb_fasta),'fasta') #load all the pdb fasta sequences
	all_human_fasta_sequences = SeqIO.parse(open(human_fasta),'fasta') #load all the pdb fasta sequences
	
	pdb_fasta = get_relevant_dict_of_fasta(all_pdb_fasta_sequences, list_of_pdb_chains) #get list of fasta seqs for the pdb_chains of interest
	target_human_fasta = get_relevant_dict_of_fasta(all_human_fasta_sequences, list_of_targ_prots) #get list of fasta seqs for the target human proteins


	for human_seq in target_human_fasta: #go through each target protein
		name = human_seq
		sequence = target_human_fasta[name]

		output_fi_name = name + '.fasta'
		pdb_out_file = name + '.pdbs.fasta'

		with open(osp.join(human_out_folder, output_fi_name), 'w') as of: #for each target protein, output its fasta seq into a separate file (i.e. uniprotid.fasta)
			of.write(">"+ name + "\n")
			of.write(sequence + "\n")

			with open(osp.join(pdb_out_folder,pdb_out_file), "w") as pdb_out:
				#for each uniprot protein: output fasta seqs for corresponding chains into a file  (i.e. uniprotid.pdbs.fasta)
				for pdb_chain in dict_of_pdbs_for_each_targ_prot[name]:
					pdb_chain_seq = pdb_fasta[pdb_chain]
					pdb_out.write(">" + pdb_chain +"\n")
					pdb_out.write(pdb_chain_seq + "\n")


def main():
	script_dir = osp.dirname(__file__)

	endo_results = osp.join(script_dir, '..', 'data', 'endo_with_verified_exo_int_res.tsv') #only known endo PPIs (no pdb annotations), target prots all participate in pathogenic and verified exogenous PPIs
	exo_results = osp.join(script_dir, '..', 'data', 'exo_verified_pathogenic_int_res.tsv') #only pathogenic and verified exogenous PPIs (we will ignore any immune and/or non-verified PPIs)

	out_endo = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'endogenous_interfacial_residues.tsv')
	out_exo = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'exogenous_interfacial_residues.tsv')
	
	pdb_fasta = osp.join(script_dir, '..', 'data', 'PDB', 'pdb_seqres.txt') #all pdb reviewed fasta seqs
	human_fasta = osp.join(script_dir, '..', 'data', 'pdb_annotations', 'fasta_files', 'human_reviewed_uniprot.formatted_header.fasta') #all human reviewed fasta seqs


	pdb_out_folder = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'pdb_fasta_folder')
	if not osp.exists(pdb_out_folder):
		os.makedirs(pdb_out_folder)

	human_out_folder = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'human_target_fasta_folder')
	if not osp.exists(human_out_folder):
		os.mkdir(human_out_folder)
	
	
	#make dictionary of exogenous and endogenous results	
	endo_dict = format_into_dict(endo_results)
	exo_dict = format_into_dict(exo_results)


	#write results into new output files (these are different from prev. files in that we group together all the results on the same chain)
	write_results(endo_dict, out_endo, exo_or_endo = "endo")
	write_results(exo_dict, out_exo, exo_or_endo = "exo")


	
	get_fasta_seqs_as_individual_files(exo_dict, endo_dict, pdb_fasta, human_fasta, pdb_out_folder, human_out_folder)

if __name__ == '__main__':
	main()

