'''

Map between UniProt and RefSeq IDs
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

from uniprot_tools import get_refseq_prot_id
from fasta_tools import parse_fasta_file_into_dict
from get_divergence_of_orthologs_in_closely_related_species import get_list_of_target_proteins
import os
import os.path as osp
from Bio import Entrez
import json 

def load_refseq_to_uniprot_mapping_dict(refseq_to_uniprot_file, list_of_proteins):
	dict_of_mappings = {}
	with open(refseq_to_uniprot_file) as rs_to_uni:
		next(rs_to_uni)
		for line in rs_to_uni:
			line = line.strip("\n")
			rs_id = line.split("\t")[0]
			uniprot_id = line.split("\t")[1]
			if uniprot_id not in list_of_proteins: #only store refseq ids that map to the uniprot ids corresponding to our target proteins of interest
				continue
			
			if uniprot_id not in dict_of_mappings:
				dict_of_mappings[uniprot_id] = [rs_id]
			else:
				dict_of_mappings[uniprot_id].append(rs_id)
	return dict_of_mappings
	

def get_refseq_prot_ids(list_of_proteins, dict_of_mappings):
	dict_of_refseq_mappings = {}
	for prot in list_of_proteins:
		#print(prot)
		if prot not in dict_of_mappings: #if not in refseq_to_uniprot_mapping file, we use the uniprot_to_refseq idmapping API tool on uniprot
			try:
				uniprot_to_rf_mapping  = get_refseq_prot_id(prot)
				if uniprot_to_rf_mapping != '':
					dict_of_refseq_mappings[prot] =  uniprot_to_rf_mapping
					
				else:
					print(f'refseq protein not found for {prot}')
			except:
				print(f'Something went wrong when using uniprot API mapping for: {prot}')
			
		else:
			mapped_refseqs = dict_of_mappings[prot]
			refseq_select = ''
			for mapped_rs in mapped_refseqs:
				handle = Entrez.efetch(db="protein", id=mapped_rs, rettype="gp")
				for line in handle.readlines():
					if line.startswith('KEYWORDS'):
						if 'Select' in line:
							refseq_select = mapped_rs
						break
			if refseq_select == '':
				refseq_select = mapped_refseqs[0] # if no refseq annotated as the 'Select', then pick first one
			
			dict_of_refseq_mappings[prot] = refseq_select

	print(f'Number of uniprot protein ids successfully mapped to refseq: {len(dict_of_refseq_mappings)}')
	# print(dict_of_refseq_mappings)

	return dict_of_refseq_mappings

def map_refseq_prot_to_refseq_nuc(refseq_prot_id):
	"""Map the refseq protein id to the refseq nucleotide (transcript) id
	
	Args:
		refseq_prot_id: refseq protein id whose nucleotide transcript we want
	"""
	try:
		nucleotide_fasta = Entrez.efetch(db="protein", id=refseq_prot_id, rettype="fasta_cds_na")
		fasta_txt = nucleotide_fasta.read()
	except:
		fasta_txt = ''

	return fasta_txt

def make_uniprot_to_refseq_nuc_mapping(dict_of_refseq_mappings):
	""" Make a dictionary where keys are uniprot protein ids and values are lists of refseq info: [refseq_prot_id, refseq_nucleotide_id, refseq_nucleotide_seq]

	"""
	dict_of_uniprot_to_refseq_prot_and_nuc_mappings = {}
	for prot in dict_of_refseq_mappings:
		refseq_prot_id = dict_of_refseq_mappings[prot]
		nuc_fasta = map_refseq_prot_to_refseq_nuc(refseq_prot_id)
		if nuc_fasta == '':
			print(f"No refseq nucleotide for the following refseq prot: {refseq_prot_id}")
			continue

		nuc_header = nuc_fasta.split()[0].split('|')[1].split('.')[0] #get just the refseq cds transcript id (e.g. NM_016937)
		nuc_fasta_seq = nuc_fasta.split('\n', 1)[1].replace("\n", '') #get the fasta sequence (without any newlines)
		dict_of_uniprot_to_refseq_prot_and_nuc_mappings[prot] = [refseq_prot_id, nuc_header, nuc_fasta_seq]

	return dict_of_uniprot_to_refseq_prot_and_nuc_mappings



def make_human_protein_and_cds_fasta_dict(dict_of_uniprot_to_refseq_prot_and_nuc_mappings, targ_prots_fasta_fi, targ_prot_out_folder, targ_prot_cds_out_folder):
	targ_prots_fasta_dict = parse_fasta_file_into_dict(targ_prots_fasta_fi, change_header=False, split_header_by='|', ind_to_keep = 1)

	for targ_prot in dict_of_uniprot_to_refseq_prot_and_nuc_mappings:
		with open(osp.join(targ_prot_out_folder, targ_prot + '.aa.fasta'), 'w') as out_targ_prot_fa, open(osp.join(targ_prot_cds_out_folder, targ_prot + '.cds.fasta'), 'w') as out_targ_prot_cds_fa:
			targ_prot_sequence = targ_prots_fasta_dict[targ_prot]
			targ_prot_cds_seq = dict_of_uniprot_to_refseq_prot_and_nuc_mappings[targ_prot][2]

			out_targ_prot_fa.write('>' + targ_prot + "\n" + targ_prot_sequence + "\n")
			out_targ_prot_cds_fa.write('>' + targ_prot + "\n" + targ_prot_cds_seq + "\n")




def main():
	script_dir = osp.dirname(__file__)
	exogenous_int_res = osp.join(script_dir, '..', 'data', 'categorized_interfacial_residues', 'v_target_h', 'exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv')
	refseq_to_uniprot_file = osp.join(script_dir, '..', 'data', 'dnds', 'gene_refseq_uniprotkb_collab')
	rbh_folder = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'v_target_h', 'rbh')
	proteomes_folder = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'proteomes')
	targ_prots_fasta_fi = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'human_target_proteins.fasta')

	#Output folders for protein and cds sequences (mkdir if doesn't exist yet)
	targ_prot_out_folder = osp.join(script_dir, '..', 'data', 'dnds', 'human_targ_prot_fasta')
	targ_prot_cds_out_folder = osp.join(script_dir, '..', 'data', 'dnds', 'human_targ_prot_cds_fasta_refseq')
	if not osp.exists(targ_prot_out_folder):
		os.mkdir(targ_prot_out_folder)
	if not osp.exists(targ_prot_cds_out_folder):
		os.mkdir(targ_prot_cds_out_folder)

	# Load the list of human target proteins (we're only focusing on virus_target_human targets)	
	list_of_human_target_proteins = get_list_of_target_proteins(exogenous_int_res)

	# Load mappings between all refseq and uniprot ids into dictionary
	print("Building uniprot-to-refseq mapping dictionary using refseq_to_uniprot file . . .")
	dict_of_mappings = load_refseq_to_uniprot_mapping_dict(refseq_to_uniprot_file, list_of_human_target_proteins)

	#Get mappings between refseq and specified uniprot ids (note: we will skip immunoglobulin proteins b/c they don't correspond to a specific cds seq)
	print("Mapping ids for: human . . .")
	dict_of_refseq_mappings_human = get_refseq_prot_ids(list_of_human_target_proteins, dict_of_mappings)

	
	#Map protein to corresponding nucleotide and output nuclotide fasta for human
	print("Mapping proteins to cds sequences . . . ")
	dict_of_uniprot_to_refseq_prot_and_nuc_mappings = make_uniprot_to_refseq_nuc_mapping(dict_of_refseq_mappings_human)


	#Output fasta files for human target proteins and their cds sequences (one seq per file)
	print("Outputting proteins and their cds sequences . . .")
	make_human_protein_and_cds_fasta_dict(dict_of_uniprot_to_refseq_prot_and_nuc_mappings, targ_prots_fasta_fi, targ_prot_out_folder, targ_prot_cds_out_folder)

	

if __name__ == '__main__':
	main()