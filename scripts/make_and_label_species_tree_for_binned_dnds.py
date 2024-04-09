import os
import os.path as osp
import argparse
from Bio import SeqIO
import re

''' Relabeling species tree (simply converting scientific name to common name)

'''

def load_species_tree(species_tree_file):
	""" Load species tree
	Args:
		species_tree_file: newick file for phylogenetic species tree
	Return:
		species_tree (loaded species tree as a string)

	"""
	species_tree = ''
	with open(species_tree_file) as s_tree:
		species_tree = s_tree.readline()
		#print(species_tree)
	return species_tree

def load_mapping_between_taxonomy_and_species_name(mapping_file):
	""" Load mapping between species common names and taxonomy names
	Args:
		mapping_file: file with species common names and corresponding taxonomy names (separated by tabs), one species per line

	Return:
		mapping_dict: dictionary mapping  common names to taxonomy names (for correspondence between our ortholog headers and the names used for making the species trees)
	"""
	mapping_dict = {}
	with open(mapping_file) as mf:
		for line in mf:
			line = line.strip("\n")
			name = line.split("\t")[0]
			taxonomy_name = line.split("\t")[1].replace(" ", '_') #add a '_' between spaces in taxonomy name to match the tree file (newick file)
			mapping_dict[name] = taxonomy_name
	return mapping_dict

def make_tree(msa, species_tree, mapping_dict, relabeled_species_tree_fi):
	""" Relabel tree by mapping from species taxonomy to protein sequence header
	Go through fasta file and save the headers, then rename the tree according to headers by using the mapping file

	Args:
		msa: multiple sequence alignment file for a specific residue type for seq header retrieva; (either of the 5 works b/c their headers are the same)
		species_tree: species tree showing the phylogenetic relationship between the species in the MSA (to be relabelled with the headers in msa)
		mapping_dict: mapping dictionary between common names and taxonomy names

		relabeled_species_tree_fi: file for storing relabeled species tree
	"""
	
	list_of_fasta_headers = []
	msa_fasta_fi = SeqIO.parse(open(msa), 'fasta')


	#print(targ_prot_name)
	for fasta in msa_fasta_fi:
		name = fasta.id
		list_of_fasta_headers.append(name)

	for header in list_of_fasta_headers:
		
		taxonomy_name = mapping_dict[header]
		species_tree = species_tree.replace(taxonomy_name, header)
	species_tree = re.sub('\'.*?\'', '', species_tree)
	
	with open(relabeled_species_tree_fi, 'w') as relabeled_tree_fi:
		relabeled_tree_fi.write(species_tree)

	

def main():
	script_dir = osp.dirname(__file__)

	species_tree_file = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'taxnames_19_orthologs.nwk') # species_tree_file (file with species tree where nodes are the species scientific names)
	
	mapping_file = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'taxnames_mapping.tsv') #mapping file that maps common species names to scientific names
	
	msa = osp.join(script_dir, '..', 'data', 'binned_dnds', 'exo.binned_codons.fasta') #could also be any of the other 6 residue types (we just want their header names, which are the same for all 7)
	relabeled_species_tree_folder = osp.join(script_dir, '..', 'data', 'binned_dnds')

	relabeled_species_tree_fi  = osp.join(relabeled_species_tree_folder, 'binned_dnds.species_tree.nwk')
	


	species_tree = load_species_tree(species_tree_file)

	mapping_dict = load_mapping_between_taxonomy_and_species_name(mapping_file)
	
	make_tree(msa, species_tree, mapping_dict, relabeled_species_tree_fi)

if __name__ == '__main__':
	main()
