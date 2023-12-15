'''

Tools for retrieving and processing UniProt data
=================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import urllib.request
import urllib.parse
import pandas as pd
from io import StringIO
import os
import requests

def get_uniprot_sequence(uniprot_id): 
	#programmatically downloading the sequence and status for a uniprot protein (id as query)
	#e.g. https://rest.uniprot.org/uniprotkb/search?query=P00750&format=tsv&fields=accession,reviewed,sequence
	url = ('https://rest.uniprot.org/uniprotkb/search?query='+ uniprot_id + '&format=tsv&fields=accession,reviewed,sequence')
	with urllib.request.urlopen(url) as response:
		res = response.read()
	df = pd.read_csv(StringIO(res.decode("utf-8")), sep="\t")
	
	return df


def obsolete_or_alt_ids(mapping_df): 
	""" Get list of uniprotids who's sequence and status are "NaN" -- likely obsolete or alternate
	-- will delete these entries when formating the data
	"""
	df_to_delete = mapping_df[mapping_df.isnull().any(axis=1)]["uniprotid"].to_list()

	return df_to_delete

def get_unreviewed(mapping_df):
	""" Get list of unreviewed uniprot ids, we will be deleting any ppis where at least one of the protein partners is unreviewed
	"""
	list_of_unreviewed_prots = []
	for index, row in mapping_df.iterrows():
		
		if row['status'] == 'unreviewed':
			if row['query'] not in list_of_unreviewed_prots: #add the unreviewed uniprot id in if not alredy there
				list_of_unreviewed_prots.append(row['query'])
	return list_of_unreviewed_prots

def get_uniprot_name_lineage(uniprot_id): 
	# programmatically downloading the taxonimic lineage of a uniprot protein (id as query)
	
	url = ('https://rest.uniprot.org/uniprotkb/search?query='+ uniprot_id + '&format=tsv&fields=accession,organism_name,protein_name,lineage')
	with urllib.request.urlopen(url) as response:
		res = response.read()
	df = pd.read_csv(StringIO(res.decode("utf-8")), sep="\t")
	
	return df

##### Functions for using downloaded uniprot data #####
def load_uniprot_info(uniprot_tsv):
	mapping_df = pd.read_csv(uniprot_tsv, sep='\t', header=0)
	return mapping_df

def obsolete_or_alt_ids2(protein_list, mapping_df):
	""" 
	Get list of uniprotids who's "From" and "Entry" column values in mapping_df are different -- alternate ids (were merged with another protein into a new id, or demerged into separate ids)
	and uniprotids that don't exist in mapping_df (not in Uniprot anymore) -- obsolete/deleted ids
	-- will delete these entries when formating the data

	Args:
		protein_list: list of proteins that we want to retrieve uniprot info from
		mapping_df: pandas df with info downloaded from uniprot (using protein_list file as input)
		--- Note: not all proteins in protein_list will be in mapping_df -- i.e. obsolete ids will not appear. Also, some proteins will also have alternate ids.
	
	Return:
		list_to_delete: list containing proteins that are obsolete or have alternate ids on uniprot
	"""
	list_to_delete = []

	# get obsolete/deleted ids
	for p in protein_list:
		if p not in mapping_df['From'].values: #this uniprot entry doesn't exist (not even with an alternate id)
			list_to_delete.append(p)

	#get uniprotids that now have alternate ids
	for index, row in mapping_df.iterrows():
		if (row['From'] != row['Entry']):
			if row['From'] not in list_to_delete:
				list_to_delete.append(row['From'])
	return list_to_delete

	
def get_unreviewed2(mapping_df):
	''' Get list of uniprot ids in mapping_df that are unreviewed
	Args:
		mapping_df: pandas df with info downloaded from uniprot 
	Return:
		unreviewed_prots: list of unreviewed proteins
	'''
	unreviewed_prots = []
	for index, row in mapping_df.iterrows():
		if row['Reviewed'] == 'unreviewed':
			if row['From'] not in unreviewed_prots: #add the unreviewed uniprot id in if not alredy there
				unreviewed_prots.append(row['From'])
	return unreviewed_prots


def write_uniprot_seqs_to_fasta(uniprotid_seq_mapping, out_fa):
	""" Write contents of pandas dataframe of uniprotid to sequence mapping into a fasta file
	Args:
		uniprotid_seq_mapping (pandas dataframe): pandas dataframe with two columns, "uniprotid" and "sequence"
		out_fa (file): fasta file to output

	"""
	with open(out_fa, 'w') as of:
		for index, row in uniprotid_seq_mapping.iterrows():
			of.write('>' + row['From'] + "\n" + row['Sequence'] + "\n")



# Return jobid for specified id mapping using uniprot API
def submit_id_mapping(from_db, to_db, ids):
	""" get jobid for id mapping job

	Args:
		from_db: database to map from (e.g. UniProtKB_AC-ID: original uniprot accession id)
		to_db: database to map to (e.g. RefSeq_Protein: ref seq protein accession)
		ids: list of protein ids that we want to map (must be in from_bd format)


	"""
	API_URL = "https://rest.uniprot.org"

	request = requests.post(
		f"{API_URL}/idmapping/run",
		data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
	)

	request.raise_for_status()
	return request.json()["jobId"]

def get_refseq_prot_id(uniprotid):
	""" Retrieve the refseq protein id from uniprot accession using uniprot api (only do this if can't find uniprot accession in our downloaded ref_seq to uniprot mapping file)
	
	Args:
		uniprotid: uniprot accession id that we want to map to refseq
	"""
	from_db="UniProtKB_AC-ID"
	to_db="RefSeq_Protein"
	ids=[uniprotid] 

	jobid = submit_id_mapping(from_db, to_db, ids)
	url = 'https://rest.uniprot.org/idmapping/stream/' + jobid + '?&format=tsv'
	#https://rest.uniprot.org/idmapping/stream/ea3489307a51b332d489b720a7bcc85f609b7909?&format=tsv

	with urllib.request.urlopen(url) as response:
		res = response.read()
	df = pd.read_csv(StringIO(res.decode("utf-8")), sep="\t")
	''' df will look something like this (column names "From" and "To") with corresponding elements in first row
	    From           To
		0  P62861  NP_001988.1

	'''
	print(df)
	if df.empty: #this uniprot id doesn't correspond to a refseq protein (e.g. it's an immunoglobuilin protein --only corresponds to genomic NG refseq regions)
		return ''

	mapped_id = df['To'][0].split(".")[0] # get the mapped id only -- first element of 'To' column (e.g. NP_001988)
	return mapped_id



	

	