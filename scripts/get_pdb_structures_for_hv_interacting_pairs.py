'''

Script that maps PDB structures to human-virus PPIs.
Downloads PDB files (mmcif).
====================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

from pathlib import Path
import os
import os.path as osp
import pickle
import Bio
from Bio.PDB import PDBList
from pdb_annotations_tools import download_pdb_files_from_file, load_pdb_resolution




def format_blast_results_into_dict(blast_results, pdb_resolution_dict):
	'''
	Formatting blast results into dict. 
	Filtering for only the alignment between the same protein and PDB chain that gives the lowest e-value (we will store only the lowest e-val).

	Args:
		blast_results(Path): path to blast results

	Return:

		formatted nested dict containing protein IDs as keys and value is a dict with aligned pdbid as keys for a nested dict containing the chains as keys and corresponding evals and coverage as value
		e.g.

			{protein1: 
				{mapped_pdbid: 
					{chain1: [best_eval_for_this_chain, coverage], 
					 chain2: [best_eval_for_this_chain, coverage],...}
				 mapped_pdbid: 
				 .
				 .
				 .
				 },
			 protein2: ...
			}
	'''


	blast_res_dict = {}

	with open(blast_results, 'r') as  b:
		for line in b:
			line = line.strip("\n")
			info = line.split("\t")
			uniprot_id = info[0]
			pdbid = info[1].split("_")[0]
			chainid = info[1].split("_")[1]
			coverage = info[2]
			e_val = float(info[10])

			#Check if this pdbid has a resolution else, skip
			if pdbid not in pdb_resolution_dict:
				print(f'PDBid does not have resolution: {pdbid}')
				continue
			#If resolution of pdbid > 6 or is equal to -1.0, then continue (don't need to store)
			# if pdb_resolution_dict[pdbid] > 6.0 or pdb_resolution_dict[pdbid] == -1.0:
			if pdb_resolution_dict[pdbid] > 6.0 or pdb_resolution_dict[pdbid] == -1.0:
				continue

			if uniprot_id not in blast_res_dict:
				blast_res_dict[uniprot_id] = {}

			if pdbid not in blast_res_dict[uniprot_id]:
				blast_res_dict[uniprot_id][pdbid] = {}

			if chainid not in blast_res_dict[uniprot_id][pdbid]:
				blast_res_dict[uniprot_id][pdbid][chainid] = [e_val, coverage]
			else: #both uniprot id and pdb and chainid are in dict already, so we update if the e-val is smaller than what's stored
				if e_val < blast_res_dict[uniprot_id][pdbid][chainid][0]:
					blast_res_dict[uniprot_id][pdbid][chainid] = [e_val, coverage]

	return blast_res_dict




def map_pdb_structures_to_human_virus_ppi(pdb_resolution_fi, interacting_hv_pair, human_blast_results, virus_blast_results, mapped_pdb_fi, list_of_pdb_structs_fi):
	''' Gets PDB structures and their corresponding chains that map to human-virus ppis and output as a tsv file
		Args:
			interacting_hv_pairs (Path): path to list of ppis
			human_blast_results (Path): path to human blast results (already filtered by e-val thresholds)
			virus_blast_results (Path): path to virus blast results (already filtered by e-val thresholds)

		Outputs:
			1. mapped_pdb_fi; tsv file with the following:
				humanid, virusid, pdbID, human protein chains, virus protein chains
			2. list of unique pdb structures in the file above
	'''
	print('### Loading pdb resolution dict ###')
	pdb_resolution_dict = load_pdb_resolution(pdb_resolution_fi)

	print('### Formatting blast results into dict ###')
	human_blast_res_dict = format_blast_results_into_dict(human_blast_results, pdb_resolution_dict)
	virus_blast_res_dict = format_blast_results_into_dict(virus_blast_results, pdb_resolution_dict)
	# print(len(human_blast_res_dict))
	# print(len(virus_blast_res_dict))

	print('### Mapping blast results to known human-virus PPIs to find representative structures ###')
	human_virus_mapped_pdb_dict = {} #dict storing the ppis and corresponding pdb chains

	list_of_unique_pdb_structs = []
	list_of_unique_ppis = []
	list_of_unique_human_prots = []
	list_of_unique_virus_prots = []

	with open(interacting_hv_pair, 'r') as hv_pairs:
		next(hv_pairs)
		with open(mapped_pdb_fi, 'wt') as outfi:
			outfi.write('human_prot' + "\t" +  'virus_prot' + "\t" + 'PDBid' + "\t" + 'human_chains' + "\t" + 'virus_chains' +"\n") #write header
			for line in hv_pairs:
				line = line.strip("\n")
				human_prot = line.split("\t")[0]
				virus_prot = line.split("\t")[1]

				human_virus_ppi = human_prot+"_"+virus_prot


				if human_prot in human_blast_res_dict and virus_prot in virus_blast_res_dict: #checking that both proteins mapped to at least one pdb structure
					for mapped_pdb in human_blast_res_dict[human_prot]: #go through each pdb structure that this human protein maps to
						if mapped_pdb in virus_blast_res_dict[virus_prot]: #human and virus pair share a pdb structure
							human_chains = list(human_blast_res_dict[human_prot][mapped_pdb].keys()) #get list of chains that mapped to the human protein
							virus_chains = list(virus_blast_res_dict[virus_prot][mapped_pdb].keys()) #get list of chains taht mapped to the virus protein
							#write to file
							outfi.write(human_prot + "\t" + virus_prot +"\t" + mapped_pdb + "\t" 
								+ ",".join(human_chains) +"\t" + ",".join(virus_chains) + "\n")
							if mapped_pdb not in list_of_unique_pdb_structs:
								list_of_unique_pdb_structs.append(mapped_pdb)

							#for stat purposes (want to see how many unique ppis, human prots and virus prots are mapped to pdb structs via blast)
							if human_virus_ppi not in list_of_unique_ppis:
								list_of_unique_ppis.append(human_virus_ppi)

							if human_prot not in list_of_unique_human_prots:
								list_of_unique_human_prots.append(human_prot)

							if virus_prot not in list_of_unique_virus_prots:
								list_of_unique_virus_prots.append(virus_prot)

	print("Number of unique pdb structures: ", str(len(list_of_unique_pdb_structs)))
	print("Number of unique ppis: ", str(len(list_of_unique_ppis)))
	print("Number of unique human proteins: ", str(len(list_of_unique_human_prots)))
	print("Number of unique virus proteins: ", str(len(list_of_unique_virus_prots)))


	with open(list_of_pdb_structs_fi, 'wt') as list_of_pdb_out:
		for pdb_struc in list_of_unique_pdb_structs:
			list_of_pdb_out.write(pdb_struc +"\n")






def main():
	script_dir = osp.dirname(__file__)

	pdb_resolution_fi = osp.join(script_dir, '../data/PDB/resolu.idx') #resolutions file

	path_to_blast_results = osp.join(script_dir, '../data/blast/results') #blast results folder

	hv_pairs_fi =  osp.join(script_dir, '../data/IntAct_human_virus/intact_exogenous_human_ppis_with_valid_ids.tsv')
	
	

    ## For human e-val threshold: 1e-10 and virus e-val threshold: 1e-5
	
	human_filtered_eval_fi = osp.join(path_to_blast_results,'human_against_pdb.blast_out.1e-10.txt') #human blast results with e-values filtered
	virus_filtered_eval_fi = osp.join(path_to_blast_results,'virus_against_pdb.blast_out.1e-5.txt') #virus blast results with e-values filtered
	

	mapped_pdb_fi = osp.join(script_dir,'../data/exo_structural_annotation/human_virus_ppis_mapped_to_pdb_chains.tsv')
	list_of_pdb_structs_fi = osp.join(script_dir,'../data/exo_structural_annotation/unique_mapped_pdb_structs.txt')
	pdb_fi_dir = osp.join(script_dir, '../data/pdb_files')

	if not osp.exists(osp.dirname(mapped_pdb_fi)):
		os.mkdir(osp.dirname(mapped_pdb_fi))

	if not osp.exists(pdb_fi_dir):
		os.mkdir(pdb_fi_dir)

	
	map_pdb_structures_to_human_virus_ppi(pdb_resolution_fi, hv_pairs_fi, human_filtered_eval_fi, virus_filtered_eval_fi, mapped_pdb_fi, list_of_pdb_structs_fi)

	print('Downloading PDB files . . .')

	download_pdb_files_from_file(list_of_pdb_structs_fi, pdb_fi_dir)
	
	print('Done!')

	
			

if __name__ == "__main__":
	main()

