'''

Make and relabel species trees for dN/dS analyses
=================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp
import argparse
from Bio import SeqIO
import re

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

def make_tree(msa, species_tree, mapping_dict, relabeled_species_tree_folder):
	""" Relabel tree by mapping from species taxonomy to protein sequence header
	Go through fasta file and save the headers, then rename the tree according to headers by using the mapping file

	Args:
		msa: multiple sequence alignment file for a specific target protein and its orthologs
		species_tree: species tree showing the phylogenetic relationship between the species in the MSA (to be relabelled with the headers in msa)
		mapping_dict: mapping dictionary between common names and taxonomy names

		relabeled_species_tree_folder: folder storing relabeled species tree
	"""
	
	list_of_fasta_headers = []
	msa_fasta_fi = SeqIO.parse(open(msa), 'fasta')

	targ_prot_name = (osp.basename(msa)).split('_')[0] #e.g. path/A0A075B6I9_msa.fasta => A0A075B6I9_msa.fasta => A0A075B6I9
	#print(targ_prot_name)
	for fasta in msa_fasta_fi:
		name = fasta.id
		list_of_fasta_headers.append(name)

	for header in list_of_fasta_headers:
		organism_name = ''
		if '_' not in header: #this is the target protein, its header is only the uniprot id (unlike the orthologous proteins which have the corresponding species in the header)
			organism_name = 'human'
		else:
			organism_name = header.split("_")[-1]
		
		taxonomy_name = mapping_dict[organism_name]
		species_tree = species_tree.replace(taxonomy_name, header)
	species_tree = re.sub('\'.*?\'', '', species_tree)
	relabeled_species_tree_fi  = osp.join(relabeled_species_tree_folder, targ_prot_name + '.species_tree.nwk')
	with open(relabeled_species_tree_fi, 'w') as relabeled_tree_fi:
		relabeled_tree_fi.write(species_tree)

	

def main():
	script_dir = osp.dirname(__file__)
	
	mapping_file = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'taxnames_mapping.tsv') #mapping file that maps common species names to scientific names
	parser = argparse.ArgumentParser()

	parser.add_argument('-m', '--msa_file')
	parser.add_argument('-s', '--species_tree_fi')
	parser.add_argument('-o', '--output_folder')

	args = parser.parse_args()
	
	msa_fasta_fi = args.msa_file #e.g. O00459_msa.fasta
	species_tree_file_name = args.species_tree_fi #e.g.  taxnames_19_orthologs.nwk
	output_folder_id = args.output_folder # e.g. output_folder_id=species_tree_out
	
	msa = osp.join(script_dir, '..', 'data', 'dnds', 'aa_msas', msa_fasta_fi) #e.g. ../data/dnds/aa_msas/O00459_msa.fasta

	relabeled_species_tree_folder = osp.join(script_dir, '..', 'data', 'dnds', output_folder_id)

	species_tree_file = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', species_tree_file_name)

	

	if not osp.exists(relabeled_species_tree_folder):
		os.mkdir(relabeled_species_tree_folder)

	species_tree = load_species_tree(species_tree_file)

	mapping_dict = load_mapping_between_taxonomy_and_species_name(mapping_file)
	
	make_tree(msa, species_tree, mapping_dict, relabeled_species_tree_folder)

if __name__ == '__main__':
	main()
