'''

Retrieving human-derived and virus-derived chains from PDB
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

from Bio import SeqIO
from fasta_tools import output_fasta_fi_from_list_of_fasta
from pdb_annotations_tools import make_pdb_chain_to_taxid_mapping_dict, get_list_of_taxids_for_viruses_affecting_humans
import os.path as osp
import os



def get_chains_with_certain_taxids(pdb_to_taxid_mapping, list_of_taxids):
	""" Get a list of pdb chains whose taxid anntotaions exist in our taxids of interest

	Args:
		pdb_to_taxid_mapping: file with pdb chain-to-taxid mappings 
		list_of_taxids: taxids of interest
	Returns:
		the list of pdb chains corresponding to the list of taxids of interest

	e.g. in terms of humans (taxid = 9606) we want all pdb chains that have 9606 annotated as their taxid
	"""
	pdb_to_taxid = make_pdb_chain_to_taxid_mapping_dict(pdb_to_taxid_mapping) #pdb -> pdb_chain -> corresponding taxid(s)
	

	list_of_pdb_chains = []
	for pdb in pdb_to_taxid:
		for chain in pdb_to_taxid[pdb]:
			if len(pdb_to_taxid[pdb][chain]) ==1: #only one taxid; not chimeric with proteins from two different organisms. We don't want anything with more than one taxid, these are chimeric proteins with more than one organism
				if pdb_to_taxid[pdb][chain][0] in list_of_taxids:
					pdb_chain = pdb+"_"+ chain
					if pdb_chain not in list_of_pdb_chains:
						list_of_pdb_chains.append(pdb_chain)

	return list_of_pdb_chains

def get_fasta(pdb_seqres, list_of_pdb_chains, outfi): 
	""" From fasta sequence file of all pdb chains, output the subset of fasta sequences corresponding to pdb chains that are in our list of interest
	Args:
		pdb_seqres: fasta file with sequences of all pdb chains in the pdb database
		list_of_pdb_chains: list of pdb chains that we're interested in
		outfi: fasta file with the sequences for the chains in the list_of_pdb_chains
	"""
	list_of_fasta_to_output = []
	fasta_sequences = SeqIO.parse(open(pdb_seqres), 'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq) #fasta.id will only take the first word in the header
		if name in list_of_pdb_chains:
			list_of_fasta_to_output.append(fasta)

	output_fasta_fi_from_list_of_fasta(list_of_fasta_to_output, outfi)
	

def check(fasta, list_of_pdb_chains):
	"""
	check which fasta seqs are missing (most likely because pdb seqres file is older and missing info on new chains)
	"""
	f = open(fasta, 'r').read()
	c = 0
	for p in list_of_pdb_chains:
		if p not in f:
			print(p)
			c+=1
	print(c)

def main():
	script_dir = osp.dirname(__file__)
	human_affecting_virus_taxids = os.path.join(script_dir, '../data/pdb_annotations/taxonomy_ids_for_viruses_infecting_humans.tab')
	pdb_to_taxid_mapping = osp.join(script_dir, '../data/pdb_annotations/pdb_chain_taxonomy_rel.tsv')
	pdb_seqres = osp.join(script_dir, '../data/PDB/pdb_seqres.txt')
	outfi_virus = osp.join(script_dir, '../data/PDB/virus_only_pdb_seqres.txt')
	outfi_human = osp.join(script_dir, '../data/PDB/human_only_pdb_seqres.txt')

	#get list of taxid(s) for the organism(s) of interest
	list_of_virus_taxids = get_list_of_taxids_for_viruses_affecting_humans(human_affecting_virus_taxids) #for human-affecting viruses
	list_of_human_taxids = ['9606'] # for humans

	#get list of pdb chains associated with the organism(s) of interest
	list_of_pdb_chains_virus = get_chains_with_certain_taxids(pdb_to_taxid_mapping, list_of_virus_taxids)
	list_of_pdb_chains_human = get_chains_with_certain_taxids(pdb_to_taxid_mapping, list_of_human_taxids)
	
	
	#get fasta files for the chains of interest
	get_fasta(pdb_seqres, list_of_pdb_chains_virus, outfi_virus)
	get_fasta(pdb_seqres, list_of_pdb_chains_human, outfi_human)

	#checking what's missing:
	print('checking virus')
	check(outfi_virus, list_of_pdb_chains_virus)
	print('checking human')
	check(outfi_human, list_of_pdb_chains_human)
	
if __name__ == '__main__':
	main()