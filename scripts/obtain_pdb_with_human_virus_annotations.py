'''

Obtain PDB chains with human and virus annotations using BLAST results and PDB-to-UniProt annotations. 
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import urllib.parse
import urllib.request
import pandas as pd
from io import StringIO
import os
import os.path as osp
from Bio import SeqIO
from pdb_annotations_tools import download_pdb_files_from_list, make_pdb_chain_to_uniprot_mapping_dict, make_pdb_chain_to_taxid_mapping_dict, make_pdb_chain_to_uniprot_blast_mapping_dict, map_locations_of_protein_seq_on_chain, get_taxid, get_list_of_taxids_for_viruses_affecting_humans, load_pdb_resolution

##read in list of virus pdb structures, check if taxid id is amongst the list of virus ids, then add into dict


def verify_correct_pdb_human_annotation(all_human_reviewed_list, uniprot):
	in_list = uniprot in all_human_reviewed_list
	return in_list

def get_human_virus_pdbs(human_affecting_virus_taxids, pdb_to_uniprot_mapping_reviewed, pdb_to_taxid_mapping, human_pdb_to_uniprot_blast, virus_pdb_to_uniprot_blast, outfi, all_human_reviewed_fasta_fi, pdb_resolution_fi):
	""" Gets list of pdb structure ids that have both human and virus annotations
	The pdb_to_uniprot_mapping_reviewed file contains entries where pdb chains correspond to reviewed proteins.
	The human_pdb_to_uniprot_blast and virus_pdb_to_uniprot_blast mapping files will be used if the pdb chain does not have a corresponding uniprot id in pdb_to_uniprot_mapping_reviewed file


	Args:

		human_affecting_virus_taxids (txt file): list of taxids for viruses affecting humans
		pdb_to_uniprot_mapping_reviewed (tsv file): pdb-to-uniprot mapping annotations for all pdb files
		pdb_to_taxid_mapping (tsv file): pdb-to-taxid mapping annotations for all pdb files
		human_pdb_to_uniprot_blast (tsv file): pdb-to-uniprot mappings from blasting human pdb and reviewed human uniprot ids 
		virus_pdb_to_uniprot_blast (tsv file): pdb-to-uniprot mappings from blasting virus pdb and reviewed virus uniprot ids 
	Output:
		outfi (txt file): list of pdb ids that have both human and virus annotations

	"""
	print("Creating initial mapping dictionaries . . .")
	#get list of human-affecting virus taxonomy ids
	list_of_virus_ids = get_list_of_taxids_for_viruses_affecting_humans(human_affecting_virus_taxids)

	pdb_to_uniprot = make_pdb_chain_to_uniprot_mapping_dict(pdb_to_uniprot_mapping_reviewed) #pdb -> pdb_chain -> corresponding uniprot id(s)
	pdb_to_uniprot_human_blast = make_pdb_chain_to_uniprot_blast_mapping_dict(human_pdb_to_uniprot_blast)
	pdb_to_uniprot_virus_blast = make_pdb_chain_to_uniprot_blast_mapping_dict(virus_pdb_to_uniprot_blast)

	pdb_to_taxid = make_pdb_chain_to_taxid_mapping_dict(pdb_to_taxid_mapping) #pdb -> pdb_chain -> corresponding taxid(s)

	#make list of all human reviewed protein uniprot ids (this is for verifying PDB annotations, because some PDB annotations could be incorrect; i.e. claims that chain is derived from human but incorrect)
	#e.g. P57083(1jma) and P03274(5fgy)
	all_human_reviewed_fasta = SeqIO.parse(open(all_human_reviewed_fasta_fi),'fasta')
	all_human_reviewed_list = []
	for fasta in all_human_reviewed_fasta:
		name = fasta.id
		if name not in all_human_reviewed_list:
			all_human_reviewed_list.append(name)

	#now we have pdb_to_taxid mapping, so go through each and only take pdb structures with human and virus proteins 
	print("Getting human-virus pdb structures . . .")
	dict_of_human_virus = {}

	for pdb in pdb_to_taxid:
		dict_of_chains = {}
		for chain in pdb_to_taxid[pdb]:
			if len(pdb_to_taxid[pdb][chain]) ==1: #only one taxid; not chimeric with proteins from two different organisms
				if pdb_to_taxid[pdb][chain][0] =='9606' or pdb_to_taxid[pdb][chain][0] in list_of_virus_ids:
					dict_of_chains[chain] = pdb_to_taxid[pdb][chain][0]

		#go through dict_of_chains and check that there are human and human-affecting virus proteins present:
		if len(set(dict_of_chains.values())) > 1: #there is at least one human protein ('9606') and at least one virus protein in this pdb structure
			dict_of_human_virus[pdb] = dict_of_chains


	dict_of_human_virus_with_uniprot_ids = {}
	
	for pdb in dict_of_human_virus:
		dict_of_human_virus_with_uniprot_ids[pdb] = {}
		for chain in dict_of_human_virus[pdb]:
			taxid = dict_of_human_virus[pdb][chain]
			#print(pdb, chain)
			found = False
			if pdb in pdb_to_uniprot:
				if chain in pdb_to_uniprot[pdb]: #some of these pdb ids and pdb chains do not have uniprot annotations, so check blast results
					uniprot_ids_associated_with_this_chain = pdb_to_uniprot[pdb][chain] # might have more than one uniprot id associated with a taxid

					if len(uniprot_ids_associated_with_this_chain) > 1: #these are the cases where a chain has more than one protein annotation
						print(f"chimeric: {pdb}, {chain}") #(stored in pdb_annos_with_more_than_one_protein_on_a_chain.txt just to see what it is) 
						continue #ignore; these are chimeric proteins (not many cases)
					
					dict_of_human_virus_with_uniprot_ids[pdb][chain] = []
					for uniprot_id in uniprot_ids_associated_with_this_chain:
						if taxid == '9606':
							#update 2023:need to verify that the annotated uniprot id is indeed a human protein
							if verify_correct_pdb_human_annotation(all_human_reviewed_list, uniprot_id): #verify that this uniprot id is indeed a human protein (PDB could have incorrect annotations)
								dict_of_human_virus_with_uniprot_ids[pdb][chain].append((uniprot_id, taxid)) #keep track of whether this is human or virus (will be easier to search for interface residues)
						else:
							dict_of_human_virus_with_uniprot_ids[pdb][chain].append((uniprot_id, taxid)) #keep track of whether this is human or virus (will be easier to search for interface residues)
						found = True
				#Now check blast results if this pdb chain is not in the pdb to uniprot annoation file
				elif taxid == '9606': #check the blast results for human chain to uniprot mappings
					if pdb in pdb_to_uniprot_human_blast:
						if chain in pdb_to_uniprot_human_blast[pdb]:
							uniprot_mapped = pdb_to_uniprot_human_blast[pdb][chain]
							dict_of_human_virus_with_uniprot_ids[pdb][chain] = [(uniprot_mapped, taxid)]
							#print(pdb, chain, uniprot_mapped, taxid)
							found = True
				elif taxid != '9606': #check the blast results for virus chain to uniprot mappings
					if pdb in pdb_to_uniprot_virus_blast:
						if chain in pdb_to_uniprot_virus_blast[pdb]:
							uniprot_mapped = pdb_to_uniprot_virus_blast[pdb][chain]
							dict_of_human_virus_with_uniprot_ids[pdb][chain] = [(uniprot_mapped, taxid)]
							#print(pdb, chain, uniprot_mapped, taxid)
							found = True
			#Now check blast results if this pdb is not in the pdb to uniprot annoation file
			elif taxid == '9606': #check the blast results for human chain to uniprot mappings
				if pdb in pdb_to_uniprot_human_blast:
					if chain in pdb_to_uniprot_human_blast[pdb]:
						uniprot_mapped = pdb_to_uniprot_human_blast[pdb][chain]
						dict_of_human_virus_with_uniprot_ids[pdb][chain] = [(uniprot_mapped, taxid)]
						#print(pdb, chain, uniprot_mapped, taxid)
						found = True
			elif taxid != '9606': #check the blast results for virus chain to uniprot mappings
				if pdb in pdb_to_uniprot_virus_blast:
					if chain in pdb_to_uniprot_virus_blast[pdb]:
						uniprot_mapped = pdb_to_uniprot_virus_blast[pdb][chain]
						dict_of_human_virus_with_uniprot_ids[pdb][chain] = [(uniprot_mapped, taxid)]
						#print(pdb, chain, uniprot_mapped, taxid)
						found = True
			
			# if found == False: #just to check which pdb chains were unable to be mapped by neither the annotation file nor the blast mapping files
			# 	print(pdb, chain)

	
	dict_of_human_virus_final = {} #to store the chains corresponding to each uniprot id
	#e.g. pdb_id: {(uniprotid1, taxid1): ['A', 'B', 'C'], (uniprotid2, taxid)}

	for pdb in dict_of_human_virus_with_uniprot_ids:
		dict_of_human_virus_final[pdb] = {}
		for chain in dict_of_human_virus_with_uniprot_ids[pdb]:
			for uniprot_taxid_pair in dict_of_human_virus_with_uniprot_ids[pdb][chain]:

				if uniprot_taxid_pair not in dict_of_human_virus_final[pdb]:
					dict_of_human_virus_final[pdb][uniprot_taxid_pair] = [chain]
				else:
					dict_of_human_virus_final[pdb][uniprot_taxid_pair].append(chain)
	#print(dict_of_human_virus_final)
	
	pdb_resolution_dict = load_pdb_resolution(pdb_resolution_fi) #2023 update

	print("Writing output ...")
	with open(outfi, 'wt') as o:
		o.write("human_prot\tvirus_prot\tPDBid\thuman_chains\tvirus_chains\n")
		for pdb in dict_of_human_virus_final:

			#2023 update: check if this pdbid has a resolution else, skip
			if pdb not in pdb_resolution_dict:
				print(f'PDBid does not have resolution: {pdb}')
				continue

			#2023 update: Don't add any entries with pdb structures that have resolution >3 or ==-1.0
			# if pdb_resolution_dict[pdb] > 3.0 or pdb_resolution_dict[pdb] == -1.0:
			if pdb_resolution_dict[pdb] > 6.0 or pdb_resolution_dict[pdb] == -1.0:
				continue

			human_prots = [uniprot for uniprot in dict_of_human_virus_final[pdb] if uniprot[1] == '9606']
			virus_prots = [uniprot for uniprot in dict_of_human_virus_final[pdb] if uniprot[1] != '9606']
			if human_prots != [] and virus_prots != []: #some of these might have annotations that don't have uniprot ids, which we skipped, so need to verify that we still have existing uniprot ids from human and virus
				for human_p in human_prots:
					human_p_id = human_p[0]
					chains_h = ','.join(dict_of_human_virus_final[pdb][human_p])
					for virus_p in virus_prots:
						virus_p_id = virus_p[0]
						chains_v = ','.join(dict_of_human_virus_final[pdb][virus_p])
						o.write(human_p_id + "\t" + virus_p_id + "\t" + pdb + "\t" + chains_h + "\t" + chains_v +"\n")


	print("Done!")

def download_pdb_files_for_annotations(pdb_hv_annotations_fi, pdb_dir):
	""" Downloads the pdb files for the pdb structures present in our annotation file
	Input: 
		pdb_hv_annotations_fi: Info of pdb chains that have human and virus annotations

	"""
	list_of_pdbs_to_download = []

	with open(pdb_hv_annotations_fi, 'r') as pdb_info:
		next(pdb_info)
		for line in pdb_info:
			line = line.strip("\n")
			pdb_id = line.split("\t")[2]
			if pdb_id not in list_of_pdbs_to_download:
				list_of_pdbs_to_download.append(pdb_id)

	download_pdb_files_from_list(list_of_pdbs_to_download, pdb_dir) #download the pdb structures

def main():
	script_dir = os.path.dirname(__file__)
	list_of_virus_taxids = os.path.join(script_dir, '..', 'data', 'pdb_annotations', 'taxonomy_ids_for_viruses_infecting_humans.tab')
	
	pdb_to_uniprot_mapping = os.path.join(script_dir, '..', 'data', 'pdb_annotations', 'pdb_chain_uniprot_reviewed_only.tsv')
	
	pdb_to_taxid_mapping = os.path.join(script_dir, '..', 'data', 'pdb_annotations', 'pdb_chain_taxonomy_rel.tsv')
	human_pdb_to_uniprot_blast = osp.join(script_dir, '..', 'data', 'pdb_annotations', 'human_pdb_chain_uniprot.blast_out.tsv')
	virus_pdb_to_uniprot_blast = osp.join(script_dir, '..', 'data', 'pdb_annotations', 'virus_pdb_chain_uniprot.blast_out.tsv')
	all_human_reviewed_fasta_fi = osp.join(script_dir, '..', 'data', 'pdb_annotations', 'fasta_files', 'human_reviewed_uniprot.formatted_header.fasta') #for verifying that PDB annotations are correct
	pdb_resolution_fi = osp.join(script_dir, '..', 'data', 'PDB', 'resolu.idx')
	outfi = os.path.join(script_dir, '../data/pdb_annotations/pdbs_with_human_virus_annotations.txt')
	pdb_dir = osp.join(script_dir, '../data/pdb_files')
	get_human_virus_pdbs(list_of_virus_taxids, pdb_to_uniprot_mapping, pdb_to_taxid_mapping, human_pdb_to_uniprot_blast, virus_pdb_to_uniprot_blast, outfi, all_human_reviewed_fasta_fi, pdb_resolution_fi)
	
	print("Downloading required pdb files . . .")
	download_pdb_files_for_annotations(outfi, pdb_dir)
if __name__ == '__main__':
	main()