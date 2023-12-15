'''

Tools for working with and annotating PDB structures.
=====================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import urllib.parse
import urllib.request
import pandas as pd
from io import StringIO
import os
from Bio.PDB import PDBList

def load_pdb_resolution(pdb_resolution_fi):
	''' Load resolutions for each PDB structure (note some structures will have more than one resolution entry; pick the smallest)
	Args:
		pdb_resolution_fi: file containing resolutions for each PDB structure 

	Return:
		pdb_resolution_dict: dictionary with resolutions for each PDB structure (key: PDBid, value: resolution) processed from resolution_fi
	''' 
	pdb_resolution_dict = {}
	with open(pdb_resolution_fi) as res_fi:
		for line in res_fi.readlines()[6:]:
			PDBid = line.split(';')[0].lower().strip()
			resol = line.split(';')[1].strip()
			if resol == '' or PDBid == '': #either resolution or PDBid field is empty
				#print(line)
				continue
			resol = float(resol)
			if PDBid in pdb_resolution_dict: #if there is already a resolution value stored for this id, check if current resol is smaller
				if resol < pdb_resolution_dict[PDBid] and resol != -1.0: #if current resolution is smaller than what's already stored (given that it's not -1.0), then store current
					pdb_resolution_dict[PDBid] = resol
			else:
				pdb_resolution_dict[PDBid] = resol
	
	return pdb_resolution_dict


def get_list_of_taxids_for_viruses_affecting_humans(infile):
	""" Get a list of taxids for viruses that affect humans from uniprot human-host viruses taxids data
	"""
	list_of_human_affecting_virus_taxids = []
	with open(infile, 'r') as inf:
		next(inf)
		for line in inf:
			taxid = line.split("\t")[0]
			if taxid not in list_of_human_affecting_virus_taxids:
				list_of_human_affecting_virus_taxids.append(taxid)
	return list_of_human_affecting_virus_taxids

def make_pdb_chain_to_uniprot_mapping_dict(pdb_to_uniprot_mapping):
	""" 
	Input: 
		pdb_to_uniprot_mapping(tsv file): file containing information on pdb chains and their corresponding uniprot ids

	Output:
		dictionary of pdb ids, their corresponding chains, and uniprot ids as follows; a chain can correspond to more than one uniprot id if it is chimeric:
		{pdbid1: {chain1:[uniprotid_chain1_1, uniprotid_chain1_2],  
				  chain2: [uniprotid], 
				...},
	     pdbid2: ...

	    }

	"""
	mapping_dict = {}
	with open(pdb_to_uniprot_mapping, 'r') as mapping_fi:
		next(mapping_fi)
		next(mapping_fi)
		for line in mapping_fi:
			line = line.strip("\n")
			info = line.split("\t")
			pdb = info[0]
			chain = info[1]
			uniprot = info[2]

			if pdb not in mapping_dict:
				mapping_dict[pdb] = {}
				mapping_dict[pdb][chain] = [uniprot]
			else:
				if chain in mapping_dict[pdb]: #more than one uniprot id associated with a pdb chain
					if uniprot not in mapping_dict[pdb][chain]:
						mapping_dict[pdb][chain].append(uniprot)
						#print(pdb, chain, mapping_dict[pdb][chain])
				else:
					mapping_dict[pdb][chain] = [uniprot]

	return mapping_dict


def make_pdb_chain_to_uniprot_blast_mapping_dict(pdb_to_uniprot_blast_mapping):
	""" 
	Input: 
		pdb_to_uniprot_blast_mapping(tsv file): file containing best blast results for pdb chains and their corresponding uniprot ids

	Output:
		dictionary of pdb ids, their corresponding chains, and uniprot ids as follows; a chain can correspond to more than one uniprot id if it is chimeric:
		{pdbid1: {chain1:[uniprotid_chain1_1, uniprotid_chain1_2],  
				  chain2: [uniprotid], 
				...},
	     pdbid2: ...

	    }

	"""
	mapping_dict = {}
	with open(pdb_to_uniprot_blast_mapping, 'r') as mapping_fi:
		next(mapping_fi)
		for line in mapping_fi:
			line = line.strip("\n")
			info = line.split("\t")
			pdb = info[0]
			chain = info[1]
			uniprot = info[2]

			if pdb not in mapping_dict:
				mapping_dict[pdb] = {}
			mapping_dict[pdb][chain] = uniprot
					
	return mapping_dict

def make_pdb_chain_to_taxid_mapping_dict(pdb_to_taxid_mapping):
	"""
	Input: 
		pdb_to_taxid_mapping(tsv file): file containing information on pdb chains and their corresponding taxonomy ids

	Output:
		dictionary of pdb ids, their corresponding chains, and taxonomy ids as follows; a chain can correspond to more than one taxid if it is chimeric with proteins from different organisms:
		{pdbid1: {chain1: [taxid_chain1,...], 
				  chain2: [taxid_chain2,...], 
				...},
	     pdbid2: ...

	    }

	"""
	
	mapping_dict = {}
	with open(pdb_to_taxid_mapping, 'r') as mapping_fi:
		next(mapping_fi)
		next(mapping_fi)
		for line in mapping_fi:
			line = line.strip("\n")
			info = line.split("\t")
			pdb = info[0]
			chain = info[1]
			taxid = info[2]

			if pdb not in mapping_dict:
				mapping_dict[pdb] = {}
				mapping_dict[pdb][chain] = [taxid]
			else:
				if chain in mapping_dict[pdb]: #more than one uniprot id associated with a pdb chain
					if taxid not in mapping_dict[pdb][chain]:
						mapping_dict[pdb][chain].append(taxid)
						#print(pdb, chain, mapping_dict[pdb][chain])
				else:
					mapping_dict[pdb][chain] = [taxid]
	return mapping_dict


def download_pdb_files_from_file(list_of_pdb_structs_fi, pdb_dir):
	""" Downloads pdb files (cif) format if it doesn't already exist in the pdb directory

	Input:
		list_of_pdb_structs_fi: file of pdb structures (one pdb id per line)
		pdb_dir: where to store the pdb files that are downloaded
	"""
	if not os.path.isdir(pdb_dir):
		os.mkdir(pdb_dir)
	list_of_pdb_structs = []
	with open(list_of_pdb_structs_fi, 'r') as list_of_pdb:
		for line in list_of_pdb:
			line = line.strip("\n")
			list_of_pdb_structs.append(line)

	'''Selecting structures from PDB'''
	pdbl = PDBList()
	
	for p in list_of_pdb_structs:
		if not os.path.isfile(os.path.join(pdb_dir, p + '.cif')):
			pdbl.retrieve_pdb_file(p, pdir=pdb_dir)


def download_pdb_files_from_list(list_of_pdb_structs, pdb_dir):
	""" Downloads pdb files (cif) format if it doesn't already exist in the pdb directory

	Input:
		list_of_pdb_structs (list): list of pdb structures (each element in lsit is one pdb id)
		pdb_dir: where to store the pdb files that are downloaded
	"""
	if not os.path.isdir(pdb_dir):
		os.mkdir(pdb_dir)

	'''Selecting structures from PDB'''
	pdbl = PDBList()
	
	for p in list_of_pdb_structs:
		if not os.path.isfile(os.path.join(pdb_dir, p + '.cif')):
			pdbl.retrieve_pdb_file(p, pdir=pdb_dir)


def map_locations_of_protein_seq_on_chain(pdb_chain_uniprot_file, pdb_id, chain_id, uniprot_id):
	""" Gives the location of a protein sequence on a specific pdb chain (in terms of the PDB file positions)
		This is for identifying the locations of the chimeric proteins
	Args:
			pdb_chain_uniprot_file (tsv file): file containing uniprot to pdb chains mappings (as well as positions)
			pdb_id (str): the pdb id of interest
			chain_id (char): the chain id of interest
			uniprot_id (str): the uniprot id of interest
	Returns:
			(start, end) tuple: The start and end positions of the specific protein sequence (defined by uniprot_id) on the chain of the pdb structure of interest
	"""
	mapping_dict = {}
	with open(pdb_chain_uniprot_file, 'r') as mapping_fi:
		next(mapping_fi)
		next(mapping_fi)
		for line in mapping_fi:
			for line in mapping_fi:
				info = line.split("\t")
				pdb = info[0]
				chain = info[1]
				uniprot = info[2]
				start_pdb = info[5]
				end_pdb = info[6]

				if pdb not in mapping_dict:
					mapping_dict[pdb] = {}
					mapping_dict[pdb][chain] = {uniprot: (start_pdb, end_pdb)}
					
				else:
					if chain in mapping_dict[pdb]: #more than one uniprot id associated with a pdb chain
						mapping_dict[pdb][chain][uniprot] = (start_pdb, end_pdb)
					else:
						mapping_dict[pdb][chain] = {uniprot: (start_pdb, end_pdb)}

	return mapping_dict[pdb_id][chain_id][uniprot_id] #returns (pdb_start, pdb_end)


def get_taxid(uniprot_ids):
	"""
	Retrieve taxids based on a list of uniprot sequence identifier.

	this function is based on
	https://www.biostars.org/p/67822/ and https://www.biostars.org/p/94422/

	Parameters:
	uniprot_ids: List, list of uniprot identifier

	Returns:
		pd.DataFrame, pandas dataframe with uniprot id column and sequence
	"""

	url = 'https://www.uniprot.org/uploadlists/'  # This is the webserver to retrieve the Uniprot data
	#https://www.uniprot.org/uniprot/?query=accession:P00750&format=tab&columns=id,lineage-id
	params = {
		'from': "ACC",
		'to': 'ACC',
		'format': 'tab',
		'query': " ".join(uniprot_ids),
		'columns': 'id,lineage-id'}

	data = urllib.parse.urlencode(params)
	data = data.encode('ascii')
	request = urllib.request.Request(url, data)
	with urllib.request.urlopen(request) as response:
		res = response.read()
	df = pd.read_csv(StringIO(res.decode("utf-8")), sep="\t")
	df.columns = ["uniprotid", "taxid", "query"]
	df['taxid'] = df['taxid'].astype(str) #convert taxid to string
	return df[["uniprotid", "taxid"]]