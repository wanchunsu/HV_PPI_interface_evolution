'''

Format alignments between UniProt sequences and their corresponding chains
Keep track of missing PDB chains (omitted because too short)
==========================================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
import os
from Bio import SeqIO
from fasta_tools import output_fasta_fi_from_list_of_fasta
def map_pdb_residues_to_uniprot(blast_result_folder):
	blast_results_dict = {} #to store {uniprot: {pdb_chain: pdb_res_id: corresponding uniprot_res_id}}
	for blast_result in os.listdir(blast_result_folder):
		
		#print(f"Currently mapping interfacial residues for: {blast_result}")
		with open(osp.join(blast_result_folder, blast_result)) as b:
			
			
			for line in b: #e.g. A0A0A0MS14	6nnj_D	50.943	106	44	1	20	117	1	106	2.16e-39	117	QMQLVQSGAEVKKTGSSVKVSCKASGYTFTYRYLHWVRQAPGQALEWMGWITPFNGNTNYAQKFQDRVTITRD--------RSMSTAYMELSSLRSEDTAMYYCAR	QGQLVQSGATTTKPGSSVKISCKTSGYRFNFYHINWIRQTAGRGPEWMGWISPYSGDKNLAPAFQDRVNMTTDTEVPVTSFTSTGAAYMEIRNLTSDDTGTYFCAK
				
				line = line.strip("\n")
				info = line.split("\t")
				uniprot_id = info[0]
				pdb_chain_id = info[1]
				uniprot_start = int(info[6])
				uniprot_end = int(info[7])
				pdb_chain_start = int(info[8])
				pdb_chain_end = int(info[9])
				uniprot_aligned = info[12]
				pdb_chain_aligned = info[13]

				uniprot_curr_idx = uniprot_start
				pdb_curr_idx = pdb_chain_start

				if uniprot_id not in blast_results_dict:
					blast_results_dict[uniprot_id] = {}
				blast_results_dict[uniprot_id][pdb_chain_id] = {}
				
				#print(uniprot_id, pdb_chain_id, str(uniprot_start), str(uniprot_end), str(pdb_chain_start), str(pdb_chain_end), uniprot_aligned, pdb_chain_aligned)
				for uniprot_aligned_res, pdb_chain_aligned_res in zip(uniprot_aligned, pdb_chain_aligned):
					#print(uniprot_aligned_res, pdb_chain_aligned_res)
					if uniprot_aligned_res == '-' and pdb_chain_aligned_res != '-': #pdb chain res maps to a gap on uniprot seq
						blast_results_dict[uniprot_id][pdb_chain_id][pdb_curr_idx] = -1 #no corresponding residue on uniprot sequence
						pdb_curr_idx +=1 #update to next index
					
					elif uniprot_aligned_res != '-' and pdb_chain_aligned_res == '-': #uniprot seq res maps to a gap on pdb chain seq, nothing to store in dict
						uniprot_curr_idx += 1 #update to next index

					elif uniprot_aligned_res == '-' and pdb_chain_aligned_res == '-': #shouldn't even see this case...
						print("Impossible, can't have gaps on both sequences!")
					
					elif uniprot_aligned_res != '-' and pdb_chain_aligned_res != '-':
						blast_results_dict[uniprot_id][pdb_chain_id][pdb_curr_idx] = uniprot_curr_idx
						pdb_curr_idx += 1 #update to next index
						uniprot_curr_idx += 1 #update to next index

	return(blast_results_dict)

def get_missing_pdb_chain_blast_results(blast_results_dict, pdb_fasta_folder):

	for uniprot in blast_results_dict:
		pdb_chains = list(blast_results_dict[uniprot].keys())

		
		pdbs_fasta = SeqIO.parse(open(osp.join(pdb_fasta_folder, uniprot + ".pdbs.fasta")),'fasta')
		list_of_missing_pdbs = []
		for fasta in pdbs_fasta:
			name = fasta.id
			if name not in pdb_chains: #this pdb_chain was not in blast results (most likely because too short)
				list_of_missing_pdbs.append(fasta)
		if list_of_missing_pdbs!= []:
			missing_pdbs_fasta = osp.join(pdb_fasta_folder, uniprot + ".missing_pdbs.fasta")
			output_fasta_fi_from_list_of_fasta(list_of_missing_pdbs, missing_pdbs_fasta)#(fasta_list, outfi)
		

def main():
	script_dir = osp.dirname(__file__)
	blast_result_folder = osp.join(script_dir, '..', 'data', 'blast', 'results', 'targ_prot_against_pdb_chains')

	targ_fasta_folder = osp.join(script_dir, '..', 'data',"exo_and_endo", "human_target_fasta_folder")
	pdb_fasta_folder= osp.join(script_dir, '..', 'data',"exo_and_endo", "pdb_fasta_folder")

	blast_results_dict = map_pdb_residues_to_uniprot(blast_result_folder)
	
	print("Done making pdb label res to uniprot res dictionary.")

	print("Getting missing pdbs that weren't in BLAST results . . . ")
	get_missing_pdb_chain_blast_results(blast_results_dict, pdb_fasta_folder)



	
if __name__ == '__main__':
	main()