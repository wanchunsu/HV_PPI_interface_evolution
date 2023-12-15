'''

Retrieve fasta sequences for human target proteins.
====================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp
from Bio import SeqIO

def get_list_of_target_proteins(infile):
	list_of_target_proteins = []
	with open(infile) as infi:
		next(infi)
		for line in infi:
			line = line.strip("\n")
			targ_prot = line.split("\t")[0]
			list_of_target_proteins.append(targ_prot)
	return list_of_target_proteins

def make_fasta_with_all_human_targ_seqs(list_of_target_proteins, fasta_folder, out_fi):

	fasta_fis = [f for f in os.listdir(fasta_folder) if f.split(".")[0] in list_of_target_proteins]
	with open(out_fi, 'w') as o:
		
		for file in fasta_fis:
			fasta_file = SeqIO.parse(open(osp.join(fasta_folder, file)),'fasta') 
			for fasta in fasta_file:
				name, sequence = fasta.id, str(fasta.seq)
				o.write(">" + name + "\n" + sequence + "\n")



def main():
	script_dir = osp.dirname(__file__)
	infile = osp.join(script_dir, '..', 'data', 'categorized_interfacial_residues', 'virus_target_human', 'exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv')
	fasta_folder = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'human_target_fasta_folder')
	out_fi = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'human_target_proteins.fasta')
	
	if not osp.exists(osp.dirname(out_fi)):
		os.makedirs(osp.dirname(out_fi))

	list_of_target_proteins = get_list_of_target_proteins(infile)

	make_fasta_with_all_human_targ_seqs(list_of_target_proteins, fasta_folder, out_fi)


	
if __name__ == '__main__':
	main()