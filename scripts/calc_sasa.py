"""
Wrapper for doing SASA and deltaSASA calculations 

=================================================
Adapted from Eric Franzosa's work (eric.franzosa@gmail.com)
(Franzosa and Xia, 2009. https://doi.org/10.1093/molbev/msp146)

With the following modifications:
	1) Load in mmCIF files using gemmi in load_structure(previously: parsed manually and only pdb files allowed, but mmCIF is the new default)
	2) Write chains of interest into new pdb files using gemmi (previously, wrote only 'ATOM' lines, but this input format is invalid for current version of DSSP)
	3) Changed DSSP to most recent version (mkdssp); see https://github.com/PDB-REDO/dssp. 
	4) Added memoization: Save DSSP outputs as json dicts (so that we don't need to re-run DSSP for chains we've already seen)
	5) chainsInComplex now a str separated by commas (for cases where chainids are not a single letter -- e.g. chains AA and AB --> 'AA,AB)
	6) Changed max SASA values to theoretical values recommended by Tien et al., 2013 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3836772/)


	Author: Wan-Chun Su (wan.su@mail.mcgill.ca)
"""


import os
import os.path as osp
import sys
import re
import subprocess
import gzip
from collections import OrderedDict
import gemmi
import json

def load_structure(path): #load structure using gemmi 

	structure = gemmi.read_structure(path)
	return structure

def write_chains(structure, chainsToWrite, chainids_shortened, numerical_chains, outPath): #write selected chains into a new pdb file
	""" Writes the specified chains into a new pdb file
	Args:
		structure: gemmi structure
		chainsToWrite: chains that we want to extract
		chainids_shortened (boolean): whether at least one of the chains in chainsToWrite is more than one char long -- will need to shorten the chain(s)
		numerical_chains (boolean): whether at least one of the chains in chainsToWrite is numerical
		outPath: where to write new pdb file
	Return:
		shortened_chain_dict: dict storing mappings between original chainid and its shortened id
		numerical_chain_dict: dict storing mappings between original chainid and its modfiied id
	"""

	shortened_chain_dict = {}
	numerical_chain_dict = {}
	numerical_chain_dict_rev = {}
	sel = gemmi.Selection('//' + chainsToWrite)
	st_se = sel.copy_structure_selection(structure)

	st_se.remove_ligands_and_waters()

	if numerical_chains: # if at least one chain is numerical
		chains_in_order = [c.name for c in st_se[0]]
		for c_orig in chains_in_order: #add a duplicate chain with a new id
			st_se[0].add_chain(st_se[0][c_orig], unique_name=True)
		for c_orig in chains_in_order: #remove the duplicated chain
			st_se[0].remove_chain(c_orig)
		
		chains_after_editing = [ce.name for ce in st_se[0]]
		numerical_chain_dict = {chains_after_editing[j]: chains_in_order[j] for j in range(len(chains_in_order))}
		print(numerical_chain_dict)

	if chainids_shortened: #if at least one of the chains in chainsToWrite is more than 1 letter long (need to shorten)
		# if one chain in chainsToWrite will rename to 'A'
		# two chains: will rename one chain to 'A' and another to 'B' (depending on order at which it appears in the structure)
		chains_in_order = [c.name for c in st_se[0]]
		
		st_se.shorten_chain_names() 

		chains_shortened = [cs.name for cs in st_se[0]]

		shortened_chain_dict = {chains_shortened[i]: chains_in_order[i] for i in range(len(chains_in_order))}
		print(shortened_chain_dict)

	st_se.write_minimal_pdb(outPath) #just write minimal_pdb (and add header manually -- required for mkdssp)
	
	#Write the header line
	try:
		h = structure.make_pdb_headers().split('\n')[0] #this is the header line
	except: #if for some reason header is not accessible, just put the pdb structure name after "HEADER" as the header line
		h = 'HEADER    '+ st_se.name 

	with open(outPath, 'r') as original: 
		data = original.read()
		#print(data)
	with open(outPath, 'w') as modified: 
		modified.write(h + "\n" + data)
	
	

	return shortened_chain_dict, numerical_chain_dict


def run_dssp(path, structure_name, chainIDs, chainids_shortened, shortened_chain_dict, numerical_chains, numerical_chain_dict, dssp_file):
	"""Run dssp on a PDB file and parse the output.

	Args:
		path (str): location of PDB file.
		structure_name (str): pdb structure name (for naming dssp output file)
		chainIDs (str): chains to return data for (separated by commas)
		chainids_shortened (boolean): whether the chainids were shortened
		shortened_chain_dict (dict): dictionary storing shortened chainids as keys and their original chainids as values
		dssp_file (path): file to store dssp results
	Returns:
		{str: list(str)/list(int)}: single letter amino acid code
									and solvent accessible area of
									residue for each chain.
	Example of how to run mkdssp on command line (output file is optional): mkdssp --output-format dssp ../data/pdb_files/1e9h.cif out_1e9h.dssp
	# Note: we've installed dssp in /usr/local/bin/ (accessible from anywhere), so don't need path
	"""
	
	if not os.path.exists(path):
		raise ValueError('Input file missing for dssp: ' + path)
	#strCommand = dsspPath + ' -i ' + path
	# update 2023: using recent dssp (mkdssp) -- see https://github.com/PDB-REDO/dssp for more info
	strCommand = 'mkdssp --output-format dssp ' + path # assuming mkdssp is added to /usr/local/bin/ (accessible from anywhere)
	cmd = subprocess.run(strCommand, shell=True, stdout=subprocess.PIPE)
	data = {'aminos': {}, 'acc': {}, 'id': {}}
	started = False
	results = cmd.stdout.decode('utf-8')
	#print(results)
	for line in results.splitlines():
		# Skip the header
		if not started:
			regex = re.compile('^  #  ')
			if regex.search(line):
				started = True
				continue
		else:
			chainID = line[11]
			if chainids_shortened: #if the chain ids were shortened, need to map them back to original chainid
				if chainID in shortened_chain_dict:
					chainID = shortened_chain_dict[chainID] #get the original chainid 
				else:
					continue
			if numerical_chains:
				if chainID in numerical_chain_dict:
					chainID = numerical_chain_dict[chainID]
				else:
					continue
			if chainID in chainIDs.split(','): #if chainID in chainIDs:
				acc = int(line[34:38])
				data['acc'].setdefault(chainID, []).append(acc)
				amino = line[13]
				data['aminos'].setdefault(chainID, []).append(amino)
				residueID = int(line[6:10])
				data['id'].setdefault(chainID, []).append(residueID)
	
	
	with open(dssp_file, 'w') as o:
		json.dump(data, o)
	return data


def rescale(values, sequence, maxSasaPerResidue):
	"""Scale solvent accessable surface area to relative solvent accessability.

	Args:
		values: solvent accessibility of each protein.
		sequence (list(str)): amino acid sequence of protein.
		maxSasaPerResidue ({str: float}): single letter amino acid code
										  --> maximum solvent accessibility
	"""
	return [min(1, values[i] / maxSasaPerResidue[sequence[i]]) for i in range(len(values))]


def sasa_scan(structure, chainIDs, chainsInComplex, normalize, tmpPath, dssp_res_folder): 
	"""

	Args:
		structure: dict(str: list(str)):
								Chain ID -> atom lines from parsed PDB file.
		chainIDs (str): chains of interest (separated by commas).
		chainsInComplex (str): chains in complex of interest (separated by commas)
		normalize (bool): normalize to get relative solvent accessability.
		tmpPath (str): path to write out temporary PDB structure.
		

	Returns:
		list(OrderedDict(int: dict)):
			PDB residue ID to dssp results of amino acid, solvent accessability
			and difference in solvent accessability between the complexed and
			uncomplexed states
	"""
	#Update 2023: Changed maxSASA values to those recommended by Tien et al., 2013 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3836772/)
	maxSasa = {'A': 129.,
			   'R': 274.,
			   'N': 195.,
			   'D': 193.,
			   'C': 167.,
			   'Q': 225.,
			   'E': 223.,
			   'G': 104.,
			   'H': 224.,
			   'I': 197.,
			   'L': 201.,
			   'K': 236.,
			   'M': 224.,
			   'F': 240.,
			   'P': 159.,
			   'S': 155.,
			   'T': 172.,
			   'W': 285.,
			   'Y': 263.,
			   'V': 174.}
	""" Previous:
	maxSasa = {'A': 115.,
			   'R': 175.,
			   'N': 154.,
			   'D': 154.,
			   'C': 118.,
			   'Q': 180.,
			   'E': 177.,
			   'G': 91.,
			   'H': 179.,
			   'I': 147.,
			   'L': 159.,
			   'K': 206.,
			   'M': 198.,
			   'F': 188.,
			   'P': 145.,
			   'S': 128.,
			   'T': 140.,
			   'W': 181.,
			   'Y': 178.,
			   'V': 135.}
	"""
	structure_name = structure.name.lower() #get pdb id

	chainids_shortened = False
	numerical_chains = False

	
	#if at least one chainid is a number, we need to rename the chain ids (this will also take care of cases where we need to shorten the chainids)
	if any(i.isdigit() for i in chainsInComplex):
		numerical_chains = True
	else: # if at least one chainid is more than one char long, we need to shorten the chain(s)
		# we don't need to check for this if numerical_chains is also True, because we will also be modifying the chain ids, so that will take care of any chains that may have more than 1 char
		for c in chainsInComplex.split(','): 
			if len(c) > 1:
				chainids_shortened = True
				break	


	
	#Windows is case-insenstive (for example: "Af" and "AF" will be seen as the same)
	# so need to check for lowercase in chainids 
	# if yes, need to add a '~' after the lowercase so that it can be differentiated from the uppercase letter:
	chains_in_complex_labelled = chainsInComplex
	if not chainsInComplex.isupper():
		
		chains_in_complex_labelled = re.sub('([a-z]{1})', r'\1~', chainsInComplex) # "Ak,AB" ==> "Ak~,AB"



	### Compute sasa for each chain in the lone state
	dsspOutput = {}
	for chainID in chainIDs.split(','):

		chainID_labelled = chainID
		if not chainID.isupper():
			chainID_labelled = re.sub('([a-z]{1})', r'\1~', chainID) ## "Af" ==> "Af~"

		dssp_fi_name = osp.join(dssp_res_folder, structure_name +'_'+ chainID_labelled +'.dssp.json')

		# Check if dssp has already been run on this pdbid and chain:
		# If yes, load the dssp results
		if osp.exists(dssp_fi_name):
			print(f'{dssp_fi_name} already exists, loading dssp results . . .')
			with open(dssp_fi_name) as df:
				dsspOutput[chainID] = json.load(df)
		else: # Otherwise, write new pdb file and run dssp on it
			print(f'Computing dssp and saving results in {dssp_fi_name} . . .')
			shortened_chain_dict1, numerical_chain_dict1 = write_chains(structure, chainID, chainids_shortened, numerical_chains, tmpPath)
			dsspOutput[chainID] = run_dssp(tmpPath, structure_name, chainID, chainids_shortened, shortened_chain_dict1, numerical_chains, numerical_chain_dict1, dssp_fi_name) # run_dssp(tmpPath, chainID, dsspPath)

	### Compute sasa for each chain in the complex
	chainsInComplex_list = chains_in_complex_labelled.split(',')
	# Both of the following give the same dssp results, so check if either exists 
	dssp_fi_name1 = osp.join(dssp_res_folder, structure_name + '_' + chainsInComplex_list[0] + '-' + chainsInComplex_list[1] +'.dssp.json') #e.g. 1e9h_C-D.dssp.json
	dssp_fi_name2 = osp.join(dssp_res_folder, structure_name + '_' + chainsInComplex_list[1] + '-' + chainsInComplex_list[0] +'.dssp.json') #e.g. 1e9h_D-C.dssp.json
	
	# Check if dssp has already been run on this pdbid and chain pair:
	# If yes, load the dssp results
	if osp.exists(dssp_fi_name1):
		print(f'{dssp_fi_name1} already exists, loading dssp results . . .')

		with open(dssp_fi_name1) as df1:
			dsspOutput[chainsInComplex] = json.load(df1)

	elif osp.exists(dssp_fi_name2):
		print(f'{dssp_fi_name2} already exists, loading dssp results . . .')

		with open(dssp_fi_name2) as df2:
			dsspOutput[chainsInComplex] = json.load(df2)

	else: # otherwise, write new pdb file and run dssp on it
		print(f'Computing dssp and saving results in {dssp_fi_name1} . . .')

		shortened_chain_dict2, numerical_chain_dict2 = write_chains(structure, chainsInComplex, chainids_shortened, numerical_chains, tmpPath)

		dsspOutput[chainsInComplex] = run_dssp(tmpPath, structure_name, chainsInComplex, chainids_shortened, shortened_chain_dict2, numerical_chains, numerical_chain_dict2, dssp_fi_name1) #run_dssp(tmpPath, chainsInComplex, dsspPath)


	### Calculate delta SASA for each residue by comparing SASA for the monomer and the complex
	results = []
	for chainID in chainIDs.split(','):
		if chainID not in dsspOutput[chainID]['id']: #if this chainID is a ligand (only has HETATM rows) -- won't have dssp results so just append an empty dict {}
			results.append({})
			print(f"Note: {structure_name}, chain {chainID} is a ligand, so skipped")
			continue
		residueIDs = dsspOutput[chainID]['id'][chainID]
		
		seq = dsspOutput[chainID]['aminos'][chainID]
		sasaMonomer = dsspOutput[chainID]['acc'][chainID]
		sasaComplex = dsspOutput[chainsInComplex]['acc'][chainID]
		if len(sasaMonomer) != len(sasaComplex):
			raise UserWarning('Something went wrong in SASA code.')
		if normalize:
			sasaMonomer = rescale(sasaMonomer, seq, maxSasa)
			sasaComplex = rescale(sasaComplex, seq, maxSasa)
		delta = [chainRSA - complexRSA for chainRSA, complexRSA in zip(sasaMonomer, sasaComplex)]
		resInfo = list(map(lambda a, b, c, d: {'aa': a, 'rsaMonomer': b, 'rsaComplex': c, 'drsa': d}, seq, sasaMonomer, sasaComplex, delta))
		
		
		results.append(OrderedDict(zip(residueIDs, resInfo)))
	
	return results