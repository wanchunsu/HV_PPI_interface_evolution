"""
Script for finding whether residues in a protein which 
interact with a partner protein are interfacial residues
=========================================================
Adapted from Leah Pollet's script build_protein_models.py from the following repo:
https://github.com/LeahPollet/interface_structure_evolution

For more info see: Pollet et al., 2022. https://doi.org/10.1016/j.jmb.2022.167750


Author: Wan-Chun Su (wan.su@mail.mcgill.ca)
"""


import argparse
import os.path as osp
import os
from calc_sasa import sasa_scan, load_structure
import json


def get_residue_info(structure, chainID, partnerChain, minDeltaRSAThreshold, tmpPath, dssp_res_folder):
	'''Finding whether residues in the protein modeled as chainID, 
	which interacts with the protein modeled as partnerChain in the structure pdbID, are interfacial residues
	- For all residues in the protein we compute 
		rsaMonomer (float)			: 	The solvent accessibility of the residue when the protein is in monomer state
		rsaComplex (float)			: 	The solvent accessibility of the residue when the protein is co-complexed with partnerChain
		drsa (float)	  			:	The change in solvent accessibility upon complex formation (rsaMonomer-rsaComplex)
		interface (bool)  			:	Whether or not the residue is interfacial (has a change in solvent accessibility upon complex formation i.e. drsa !=0)

	Args:
		structure (gemmi structure)	:	Pre-loaded structure	
		chainID (str)				:	Chain ID for the protein of interest to the calculations
		partnerChain (str)			:	Chain ID for the partner protein in the PPI ('' if we are modeling a single protein)
		minDeltaRSAThreshold		: 	Threshold for a residue to be interfacial (i.e. residues with delta RSA >= deltaRSAThreshold are considered interfacial) 
		tmpPath (path)				: 	Path to file to store newly created pdb files (will overrride this each time)
		dssp_res_folder (path)		: 	Folder to store the dssp results
	Returns:
		sasa (OrderedDict): Format: { residueChainID : {'aa':'M','rsaMonomer': float, 'rsaComplex': float, 'drsa': float, 'interface': bool}, ...}
		list_of_interfacial_residues (list): list of interfacial residue ids
	--------------------------
	Adapted from calculate_structural_properties function in build_protein_models.py written by Leah Pollet
	'''
	
	# (1) Compute solvent accessibilities (rsaMonomer,rsaComplex,drsa) for all residues in the protein
	sasa = sasa_scan(structure, chainID, chainID + ',' + partnerChain, True, tmpPath, dssp_res_folder)[0]
	if sasa == {}: #if sasa for this chain pair can't be calculated (i.e. at least one chain is a ligand -- only has HETATM entries)
		return {}, []

	# (2) Label interfacial residues
	for k, v in sasa.items():
		if sasa[k]['drsa'] >= minDeltaRSAThreshold: 
			sasa[k]['interface'] = True
		else:
			sasa[k]['interface'] = False

	# (3) Get list of residues that are interfacial residues
	list_of_interfacial_residues = [k for k in sasa if sasa[k]['interface'] == True] 
	

	return sasa, list_of_interfacial_residues


def make_dict_of_unique_pdb_chain_pairs(ppis_mapped_to_pdb_fi):
	'''
	Args:
		ppis_mapped_to_pdb_fi (tsv file): file containing PPIs mapped to pdb chains
		
			Sample ppis_mapped_to_pdb_fi file:
			human_prot	virus_prot	PDBid	human_chains	virus_chains
			P52292	P03427	4uaf	B	E
			P52292	P03427	4uae	A	F
			P52292	P03427	4uad	A	E
			P52292	P03427	2jdq	B,A	E,D
			O00505	P03427	4uae	A	F
		
	Return: 
		dict_of_pdb_chain_pairs (dict): dictionary storing the unique pdb chain pairs to check for interfacial residues

			Format of dict_of_pdb_chain_pairs:
			{	pdbid1: 	[	(chainid1, partnerchain1), 
								(chainid2, partnerchain2), ...], 
				pdbid2: 	[...], 
				...
			}
		list_of_pdb_chain_pairs (list): list storing the unique pdb chain pairs to check for interfacial residues

			Format of list_of_pdb_chain_pairs:
			[pdbid1_chainid1_partnerchainid1, pdbid2_chainid2_partnerchainid2, ... ]
			e.g. 
			['1egh_D_C', '1egh_A_B', '2awl_A,C', '2awl_C_A', '7o7y_AB_BA'...]
	'''
	dict_of_pdb_chain_pairs = {}

	with open(ppis_mapped_to_pdb_fi) as m:
		next(m) #skip header
		for line in m:
			line = line.strip("\n")
			info = line.split("\t")
			pdbid = info[2]
			target_chains = info[3].split(",")
			partner_chains = info[4].split(",")
			
			for tc in target_chains:
				for pc in partner_chains:
					if tc == pc: # don't check pairs with the same chain
						continue
					if pdbid not in dict_of_pdb_chain_pairs:
						dict_of_pdb_chain_pairs[pdbid] = []
					if (tc, pc) not in dict_of_pdb_chain_pairs[pdbid]:
						dict_of_pdb_chain_pairs[pdbid].append((tc, pc))

	return dict_of_pdb_chain_pairs
			


def get_interfacial_residues_for_unique_chain_pairs(int_res_memoize_fi, dict_of_pdb_chain_pairs, minDeltaRSAThreshold, tmpPath, dssp_res_folder):
	"""
	Go through all unique chain pairs get their interfacial residues if not already in int_res_memoize_fi

	Args:
		int_res_memoize_fi (json file)	: 	json file storing dictionary of memoized interfaicla residues (will be loaded as int_res_memoize_dict)
		dict_of_pdb_chain_pairs (dict)	: 	dicitonary containing unique pdb chain pairs to check for interfacial residues
		minDeltaRSAThreshold (float)	: 	Threshold for a residue to be interfacial (i.e. residues with delta RSA >= deltaRSAThreshold are considered interfacial) 
		tmpPath (path)					: 	Path to file to store newly created pdb files (will overrride this each time)
		dssp_res_folder (path)			: 	Folder to store the dssp results

	Return:
		int_res_memoize_dict (dict)		: 	dictionary storing the memoized interfacial residues 

	Steps:
	1. First: load int_res_memoize fi as int_res_memoize_dict (if exists)
	2. Run through dict_of_pdb_chain_pairs
	3. Check if each chain pair exists in int_res_memoize_dict:
		a. if yes: continue
		b. o/w: run dssp and store results in int_res_memoize_dict
	4. Update the int_res_memoize_fi with the new int_res_memoize_dict 
	5. Return int_res_memoize_dict (we will parse this dict in the next function to output interfacial residues for our pdb chain pairs of interest)
	
	Format of int_res_memoize_dict:
		int_res_memoize_dict = 
		{	pdbid1_chainid1_partnerchainid1: [intres, intres, ...],
			pdbid2_chainid2_partnerchainid2: [intres, intres, ...],
			...
		}


	"""
	
	int_res_memoize_dict = {}
	if osp.exists(int_res_memoize_fi):
		with open(int_res_memoize_fi) as mem:
			int_res_memoize_dict = json.load(mem)

	for pdb in dict_of_pdb_chain_pairs:
		structure = load_structure(osp.join('..', 'data', 'pdb_files', pdb + '.cif')) #load structure 
		for chain_pair_tuple in dict_of_pdb_chain_pairs[pdb]:
			chainID = chain_pair_tuple[0] 
			partnerChain = chain_pair_tuple[1]

			pdb_chain_pair = pdb + '_' + chainID + '_' + partnerChain  #to check/store in int_res_memoize_dict

			if pdb_chain_pair in int_res_memoize_dict:
				print(f'Found {pdb_chain_pair} in memoized dict, no need to re-run dssp')
				continue # already exists in int_res_memoize_dict, so no need to run dssp again
			else: # We haven't seen this pdb chain pair, so run dssp and store in pdb_chain_pair
				# Run get_residue_info to get list of interfacial residues (we probably don't need the sasa dict -- this is just for verification purposes)
				print(f'{pdb_chain_pair} not in memoized dict, running dssp and storing. . .')
				sasa, list_of_interfacial_residues = get_residue_info(structure, chainID, partnerChain, minDeltaRSAThreshold, tmpPath, dssp_res_folder)
				
				if sasa == {}: #if sasa for this chain pair can't be calculated (i.e. at least one chain is a ligand -- only has HETATM entries)
					continue 
					
				int_res_memoize_dict[pdb_chain_pair] = list_of_interfacial_residues
	
	# Update the int_res_memoize_fi with the new int_res_memoize_dict
	with open(int_res_memoize_fi, 'w') as mem_fi_write:
		json.dump(int_res_memoize_dict, mem_fi_write)

	return int_res_memoize_dict


def parse_int_res_results_into_fis(ppis_mapped_to_pdb_fi, int_res_memoize_dict, int_res_results_fi):
	""" Output int res results for each targ_prot-partner-prot pdb chain pair 

	Args:
		ppis_mapped_to_pdb_fi (tsv file): 	file containing PPIs mapped to pdb chains
		int_res_memoize_dict (dict)		: 	dictionary storing the memoized interfacial residues 
		int_res_results_fi (tsv file)	: 	output file storing the interfacial residues for each targ_prot-partner-prot pdb chain pair 


	"""
	with open(ppis_mapped_to_pdb_fi) as mapped_fi, open(int_res_results_fi, 'w') as outfi:
		outfi.write("target_prot\tpartner_prot\tPDBid\ttarget_chainid\tpartner_chainid\ttarget_interfacial_residues\n")
		
		next(mapped_fi)
		for line in mapped_fi:

			line = line.strip("\n")
			info = line.split("\t")
			target_prot = info[0]
			partner_prot = info[1]
			PDBid = info[2]
			target_chains = info[3].split(",") #list of chains associated with this human prot
			partner_chains = info[4].split(",") #list of chains associated with this virus prot
			
			for t in target_chains:
				for p in partner_chains:
					if t == p: # don't check pairs with the same chain
						continue

					#Get interfacial residues for the pdb chain pair from int_res_memoize_dict
					pdb_chain_pair_to_find = PDBid + '_' + t + '_' + p
					if pdb_chain_pair_to_find not in int_res_memoize_dict: #for cases where one chain or more is a ligand (so can't calculate dssp and hence no interfacial residues)
						continue
					int_res_list = int_res_memoize_dict[pdb_chain_pair_to_find]
					if int_res_list == []:
						continue #no need to store if no interfacial residues for this pdb chain pair
					else:
						#Write (tab separated line) into int_res_results_fi:
						#target_prot, partner_prot, PDBid, target_chainid, partner_chainid, and interfacial_residues
						int_res_list_str = [str(r) for r in int_res_list]
						outfi.write(target_prot + "\t" + partner_prot + "\t" + PDBid + "\t" + t + "\t" + p + "\t" + ",".join(int_res_list_str) + "\n")

def get_stats_for_interface_residues_exogenous(exo_interface_residues_fi):
    """
    Prints out stats related to the exogenous interface residues

    Args:
        interface_residues_fi (path): tsv fi with the interacting proteins, their corresponding pdb structures, chains, and interfacial residues on the target chain
    """
    list_of_unique_ppis = []
    list_of_unique_pdb_structs = []
    list_of_unique_human_prots = []
    list_of_unique_virus_prots = []

    with open(exo_interface_residues_fi, 'r') as int_res:
        next(int_res)
        for line in int_res:
            info  = line.split("\t")
            human_prot = info[0]
            virus_prot = info[1]
            if human_prot not in list_of_unique_human_prots:
                list_of_unique_human_prots.append(human_prot)

            if virus_prot not in list_of_unique_virus_prots:
                list_of_unique_virus_prots.append(virus_prot)

            human_virus_prot_ids = human_prot+"_"+virus_prot
            if human_virus_prot_ids not in list_of_unique_ppis:
                list_of_unique_ppis.append(human_virus_prot_ids)

            PDBid = info[2]
            if PDBid not in list_of_unique_pdb_structs:
                list_of_unique_pdb_structs.append(PDBid)

    print("Number of unique ppis: ", str(len(list_of_unique_ppis)))

    print("Number of unique pdb structures: ", str(len(list_of_unique_pdb_structs)))
    
    print("Number of unique human proteins: ", str(len(list_of_unique_human_prots)))
    print("Number of unique virus proteins: ", str(len(list_of_unique_virus_prots)))


def get_stats_for_interface_residues_endogenous(endo_interface_residues_fi):
    """
    Prints out stats related to the interface residues (endogenous)

    Args:
        interface_residues_fi (path): tsv fi with the interacting proteins, their corresponding pdb structures, chains, and interfacial residues on the target chain
    """
    list_of_unique_ppis = []
    list_of_unique_pdb_structs = []
    list_of_unique_target_prots = []
    list_of_unique_partner_prots = []

    with open(endo_interface_residues_fi, 'r') as int_res:
        next(int_res)
        for line in int_res:
            info  = line.split("\t")
            target_prot = info[0]
            partner_prot = info[1]
            if target_prot not in list_of_unique_target_prots:
                list_of_unique_target_prots.append(target_prot)

            if partner_prot not in list_of_unique_partner_prots:
                list_of_unique_partner_prots.append(partner_prot)

            

            target_partner_prot_ids = target_prot+"_"+partner_prot
            partner_target_prot_ids = partner_prot + "_" + target_prot

            if target_partner_prot_ids not in list_of_unique_ppis and partner_target_prot_ids not in list_of_unique_ppis:
                list_of_unique_ppis.append(target_partner_prot_ids)

            PDBid = info[2]
            if PDBid not in list_of_unique_pdb_structs:
                list_of_unique_pdb_structs.append(PDBid)

    print("Number of unique ppis: ", str(len(list_of_unique_ppis)))

    print("Number of unique pdb structures: ", str(len(list_of_unique_pdb_structs)))

    print("Number of unique target proteins: ", str(len(list_of_unique_target_prots)))
    
    print("Number of unique partner proteins: ", str(len(list_of_unique_partner_prots)))


def main():
	script_dir = osp.dirname(__file__)
	parser = argparse.ArgumentParser()

	#read in pdb_chain_pairs file and output file as arguments! 
	#So that we can use this script for other input files too!
	parser.add_argument('-i', '--input_ppi_to_pdb_file')
	parser.add_argument('-o', '--output_interfacial_residues_file')
	parser.add_argument('-t', '--type_of_interaction')
	args  = parser.parse_args()

	ppis_mapped_to_pdb_fi = args.input_ppi_to_pdb_file
	int_res_results_fi = args.output_interfacial_residues_file
	endo_or_exo = args.type_of_interaction

	# make folder to store interfacial residues (if doesn't already exist)
	if not osp.exists(osp.dirname(int_res_results_fi)):
		os.mkdir(osp.dirname(int_res_results_fi))

	# file to store memoized pdb chain pairs and their interfacial residues
	# we can just load this in each time, so that we don't have to re-compute the interfacial residues
	# for pdb chain pairs that we've already seen
	int_res_memoize_fi = osp.join('..', 'data', 'memoized', 'int_res_memoize.json')
	if not osp.exists(osp.dirname(int_res_memoize_fi)):
		os.mkdir(osp.dirname(int_res_memoize_fi))

	# temporary path to store the newly-made pdb files required for computing dssp 
	# Note: this file gets overwritten each time with create a new pdb file
	tmpPath = osp.join('..', 'data', 'tmp_pdb_fi', 'tmp.pdb')
	if not osp.exists(osp.dirname(tmpPath)):
		os.mkdir(osp.dirname(tmpPath))

	# folder for storing all dssp results
	dssp_res_folder = osp.join('..', 'data', 'dssp_results')
	if not osp.exists(dssp_res_folder):
		os.mkdir(dssp_res_folder)
	
	#set delta RSA threshold (if residues have delta RSA >= minDeltaRSAThreshold, they are considered interfacial residues)
	# delta RSA is the difference in RSA of the residue in the unbound and bound (in complex with its partner) state
	minDeltaRSAThreshold = 0.001 # As 


	########################## Calling functions ##########################
	
	#e.g. For H-V known PPIs to pdb: 
	#python get_interfacial_residues_using_sasa.py -i ../data/exo_structural_annotation/human_virus_ppis_mapped_to_pdb_chains.tsv -o ../data/exo_structural_annotation/interface_residues/exo_hv_interface_residues.tsv -t exo
	
	# Step 1: First make dict of unique pdb chain pairs that we would like to find interfacial residues for:
	dict_of_pdb_chain_pairs = make_dict_of_unique_pdb_chain_pairs(ppis_mapped_to_pdb_fi)

	"""
	Step 2:
		Next, load in the interfacial residues memoize dictionary from file (if exists) and run through dict_of_pdb_chain_pairs
		If a pdb chain pair has not been stored yet, we run dssp and store in int_res_memoize_dict
		Then update the int_res_memoize_fi with the new int_res_memoize_dict (for future use)
		Return new int_res_memoize_dict
	"""
	print("##### Running interfaical residues for chains pairs if not already in int_res_memoize_dict #####\n")
	int_res_memoize_dict = get_interfacial_residues_for_unique_chain_pairs(int_res_memoize_fi, dict_of_pdb_chain_pairs, minDeltaRSAThreshold, tmpPath, dssp_res_folder)

	# Step 3: Parse new int_res_memoize_dict & output interfacial residues for each targ_prot-partner-prot pdb chain pair in our ppis_mapped_to_pdb_fi
	print("\n##### Parsing interfacial residue results into file #####\n")
	parse_int_res_results_into_fis(ppis_mapped_to_pdb_fi, int_res_memoize_dict, int_res_results_fi)

	print("\nDone!")

	# Print out stats for interfacial residues
	print(f"##### Stats for {int_res_results_fi} #####\n")
	if endo_or_exo == 'exo':
		get_stats_for_interface_residues_exogenous(int_res_results_fi)
	elif endo_or_exo == 'endo': #this will ensure that we don't overcount duplicate PPIs (i.e. A-B and B-A are the same PPI)
		get_stats_for_interface_residues_endogenous(int_res_results_fi)


if __name__ == '__main__':
	main()