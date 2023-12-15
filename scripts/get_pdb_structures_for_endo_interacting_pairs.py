'''

Script that maps PDB structures to human-human PPIs.
Downloads PDB files (mmcif).
====================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
import os
from pdb_annotations_tools import download_pdb_files_from_list, load_pdb_resolution

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

			#2023 Update: Check if this pdbid has a resolution else, skip
			if pdbid not in pdb_resolution_dict:
				print(f'PDBid does not have resolution: {pdbid}')
				continue
			#2023 update: If resolution of pdbid > 6 or is equal to -1.0, then continue (don't need to store)
			# if pdb_resolution_dict[pdbid] > 3.0 or pdb_resolution_dict[pdbid] == -1.0:
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


def map_pdb_structures_to_human_human_ppis(pdb_resolution_fi, interacting_pairs, blast_results, mapped_pdb_fi):
	""" Gets PDB structures and their corresponding chains that map to human-human ppis (where at least one human protein is a target protein)
	
	Args:
		interacting_pairs: list of interacting human protein pairs
		blast_results: path to human protein and PDB blast results (already filtered by e-value threshold)
		mapped_pdb_fi: output tsv file with the following: target_human_id, human_partner_od, pdbID, target_pdb_chains, partner_pdb_chains
		

	"""
	print('### Loading pdb resolution dict ###')
	pdb_resolution_dict = load_pdb_resolution(pdb_resolution_fi)



	list_of_unique_pdb_structs = []
	list_of_unique_ppis = []
	list_of_unique_target_prots = []
	list_of_unique_partner_prots = []

	print('### Formatting blast results into dict ###')
	blast_res_dict = format_blast_results_into_dict(blast_results, pdb_resolution_dict)

	print('### Mapping blast results to known human-human PPIs to find representative structures ###')
	human_ppis_mapped_pdb_dict = {} # dict for storing the ppis and corresponding pdb chains

	
	with open(interacting_pairs, 'r') as pairs:
		next(pairs)
		with open(mapped_pdb_fi, 'w') as outfi:
			outfi.write('target_prot' + "\t" +  'partner_prot' + "\t" + 'PDBid' + "\t" + 'target_chains' + "\t" + 'partner_chains' +"\n") #write header
			for line in pairs:
				line = line.strip("\n")
				target_prot = line.split("\t")[0]
				partner_prot = line.split("\t")[1]

				endogenous_ppi_pair = (target_prot, partner_prot)

				if target_prot in blast_res_dict and partner_prot in blast_res_dict: #checking that both proteins mapped to at least one pdb structure
					for mapped_pdb in blast_res_dict[target_prot]: #go through each pdb structure that this target protein maps to
						if mapped_pdb in blast_res_dict[partner_prot]: #both target and partner protein map to the same pdb structure)

							human_target_chains = list(blast_res_dict[target_prot][mapped_pdb].keys())
							partner_target_chains = list(blast_res_dict[partner_prot][mapped_pdb].keys())

							outfi.write(target_prot + "\t" + partner_prot + "\t" + mapped_pdb + "\t" 
								+ ",".join(human_target_chains) + "\t" + ",".join(partner_target_chains) + "\n")

							if mapped_pdb not in list_of_unique_pdb_structs:
								list_of_unique_pdb_structs.append(mapped_pdb)

							#check whether a pair of interacting proteins (order doesn't matter) are already in list of unique ppis
							if (endogenous_ppi_pair[0], endogenous_ppi_pair[1]) not in list_of_unique_ppis and (endogenous_ppi_pair[1], endogenous_ppi_pair[0]) not in list_of_unique_ppis:
								list_of_unique_ppis.append(endogenous_ppi_pair)

							if target_prot not in list_of_unique_target_prots:
								list_of_unique_target_prots.append(target_prot)

							if partner_prot not in list_of_unique_partner_prots:
								list_of_unique_partner_prots.append(partner_prot)
	list_of_repeated = [] # interacting pairs where both proteins are target proteins, checking just for reference
	for pair in list_of_unique_ppis:
		p1 = pair[0]
		p2 = pair[1]
		if (p2, p1) in list_of_unique_ppis:
			if (p1,p2) not in list_of_repeated and (p2, p1) not in list_of_repeated:
				list_of_repeated.append((p1,p2))
	for i in list_of_repeated:
		print(i)
	print("Number of unique pdb structures: ", str(len(list_of_unique_pdb_structs)))
	print("Number of unique ppis: ", str(len(list_of_unique_ppis)))
	print("Number of unique target proteins: ", str(len(list_of_unique_target_prots)))
	print("Number of unique partner proteins: ", str(len(list_of_unique_partner_prots)))



def download_pdb_files(pdb_human_human_mapped_fi, pdb_dir):
	""" Downloads the pdb files for the pdb structures present in our annotation file
	Input: 
		pdb_hv_annotations_fi: Info of pdb chains that have human and virus annotations

	"""
	list_of_pdbs_to_download = []

	with open(pdb_human_human_mapped_fi, 'r') as pdb_info:
		next(pdb_info)
		for line in pdb_info:
			line = line.strip("\n")
			pdb_id = line.split("\t")[2]
			if pdb_id not in list_of_pdbs_to_download:
				list_of_pdbs_to_download.append(pdb_id)

	download_pdb_files_from_list(list_of_pdbs_to_download, pdb_dir) #download the pdb structures



def main():
	script_dir = osp.dirname(__file__)

	pdb_resolution_fi = osp.join(script_dir, '../data/PDB/resolu.idx') #resolutions file

	human_blast_filtered_eval_fi = osp.join(script_dir, '..', 'data', 'blast', 'results', 'human_endo_prots_against_pdb.blast_out.1e-10.txt')

	interacting_pairs = osp.join(script_dir, '..', 'data', 'IntAct_human_human', "intact_endogenous_human_ppis_with_valid_ids.tsv")

	
	output_dir = osp.join(script_dir, '..', 'data', 'endo_structural_annotation')
	pdb_dir = osp.join(script_dir, '../data/pdb_files')
	if not osp.exists(output_dir):
		os.mkdir(output_dir)

	mapped_pdb_fi = osp.join(output_dir, 'human_endogenous_ppis_mapped_to_pdb_chains.tsv')
	
	map_pdb_structures_to_human_human_ppis(pdb_resolution_fi, interacting_pairs, human_blast_filtered_eval_fi, mapped_pdb_fi)
	
	#download pdb files
	print("Downloading required pdb files . . .")
	download_pdb_files(mapped_pdb_fi, pdb_dir)

	print("Done!")

if __name__ == '__main__':
	main()