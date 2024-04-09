import os
import os.path as osp
import json
from Bio import SeqIO
from Bio import AlignIO
from process_dn_ds import parse_ensembl_tblastn_results, parse_human_tblastn_results, get_translated_aa_msa_alignment_fasta
from fasta_tools import parse_fasta_file_into_dict
'''
1) Go through amino acids alignments (following along with codon alignments (take 3 at a time)) and keep track of indices on the actual protein itself 
-- map from msa aa (keeping track of codon for each aa) --> targ prot aa --> residue type (if maps to targ prot) --> place into

Go through aa alignment, if aa residue in aa alignment is part of aa sequence, store (aa residue num -- on targ prot), (aa residue), (corresponding codon)
targ_prot1: {targ_prot1:{res_num_on_targ_prot: human:[], }

Then go through the indices in the above dictionary place into dict for exo, endo, mim, surface, buried, all_exo, all_endo
exo:{targ_prot1:[codon1_in_exo....], baboon_ortholog_for_targ_prot1:[codon1_in_exo]}


'''
def annotate_residues_with_category(exo_residues, endo_residues, mim_residues, surface_residues, buried_residues, all_exo_residues, all_endo_residues, dna_msas):
	''' For every target protein, we annotate all of its residues as exo, endo, mimicry, surface, buried. We will build a dict storing the target proteins, their residues (those that map to PDB), and residue annotation
	e.g. {targprot1: {res1: 'surface', res2: 'surface', res3: 'exo', ...}, targprot2: {...}, ...}

	
	'''

	list_of_dna_msas = os.listdir(dna_msas)

	list_of_targ_prots_for_dnds = [targ_prot_msa_fi.split('_')[0] for targ_prot_msa_fi in list_of_dna_msas] #get list of target prots
	
	dict_of_targprot_residue_annos = {}

	dict_of_labels = {exo_residues: 'exo', endo_residues: 'endo', mim_residues: 'mimicry', surface_residues: 'surface', buried_residues: 'buried', all_exo_residues: 'all_exo', all_endo_residues: 'all_endo'}


	for residue_cat in [exo_residues, endo_residues, mim_residues, surface_residues, buried_residues, all_exo_residues, all_endo_residues]:

		with open(residue_cat) as res:
			next(res)
			for line in res:
				line = line.strip("\n")
				targ_prot = line.split("\t")[0]

				if targ_prot not in dict_of_targprot_residue_annos:
					dict_of_targprot_residue_annos[targ_prot] = {}

				if targ_prot not in list_of_targ_prots_for_dnds: #case where protein is an immunoglobulin (so doesn't have cds sequence -- thus ignored for dnds analysis)
					continue
				if line.split("\t")[1] == '':
					continue
				residues = line.split("\t")[1].split(",") #split up residues separated by "," into a list of residues
				residues = [int(r) for r in residues]

				if residues == []:
					continue
				
				for r in residues:
					if r not in dict_of_targprot_residue_annos[targ_prot]:

						dict_of_targprot_residue_annos[targ_prot][r] = []
					dict_of_targprot_residue_annos[targ_prot][r].append(dict_of_labels[residue_cat])

	# for t in dict_of_targprot_residue_annos:
	# 	for r in dict_of_targprot_residue_annos[t]:
	# 		if len(dict_of_targprot_residue_annos[t][r])>1:
	# 			print(dict_of_targprot_residue_annos[t][r])

	#print(dict_of_targprot_residue_annos)
	return dict_of_targprot_residue_annos

def get_codon_mappings(aa_msas, dna_msas):
	""" Map aa msa positions and residues to codon msa. Parse through both aa alignment and codon alignment to match. 

	Args:
		aa_msas: folder containing aa msas (one for each target protein)
		dna_msas: folder containing codon msas (one for each target protein)

	Return:
		codon_mappings_dict: dictionary mapping aa msa positions to their corresponding aa residue and codon


	"""
	codon_mappings_dict = {} #targprot:{msa_aa_res_pos: {targprot:(aa_res, codon_for_res_pos), baboon:(aa_res, codon_for_res_pos), ...}}
	#print(os.listdir(aa_msas))
	for targ_prot_aa_msa_fi in os.listdir(aa_msas) :
		if targ_prot_aa_msa_fi.endswith('_msa.fasta'):
			targ_prot_id = targ_prot_aa_msa_fi.split('_')[0]
			targ_prot_dna_msa_fi = targ_prot_id + '_dna_msa.fasta'
			aa_alignment = AlignIO.read(osp.join(aa_msas, targ_prot_aa_msa_fi), "fasta")
			dna_alignment = AlignIO.read(osp.join(dna_msas, targ_prot_dna_msa_fi), 'fasta')
			codon_mappings_dict[targ_prot_id] = {} #initialize inner dict for current targ prot
			aa_pos = 1
			headers = [dna_alignment[i].id for i in range(len(aa_alignment[:,0]))]
			headers = [h.split('_')[-1] for h in headers]
			
			
			while aa_pos <= len(aa_alignment[0]):
				dna_pos = aa_pos * 3 - 2 #update codon start position in msa
				
				list_of_aa_residues_at_pos = aa_alignment[:,aa_pos-1]
				list_of_codons_at_pos = dna_alignment[:, dna_pos-1:dna_pos-1+3]
				
				
				codon_mappings_dict[targ_prot_id][aa_pos] = {}
			
				for org_ind in range(len(headers)):
					# for each targ prot aa_pos in msa for each organism, store the corresponding aa residue at aa_pos and corresponding codon as a tuple
					codon_mappings_dict[targ_prot_id][aa_pos][headers[org_ind]] = (list_of_aa_residues_at_pos[org_ind], str(list_of_codons_at_pos[org_ind].seq)) #store 

				aa_pos += 1 #update aa msa pos
				
			#print(codon_mappings_dict[targ_prot_id])
	#print(len(codon_mappings_dict))
	return codon_mappings_dict
			


def map_translated_aa_residues_to_msa_aligned_residues(translated_aa_tblastn_alignment_with_orig_prot, codon_mappings_dict, targ_prot, aa_msa_file):
	""" map aa residues and codons to translated cds alignment positions 

	Args: 
		translated_aa_tblastn_alignment_with_orig_prot: translated cds seq aligned with orig prot -- need this to map residues and dnds vals to orig protein
		codon_mappings_dict: dictionary containing positions in the msa and their corresponding dN/dS values
		aa_msa_file: file containing aa msa for the current protein (this will be for finding special cases for alignment positions -- dS==0)
	Return:
		dict_of_translated_cds_pos_and_aa_codon: dictionary storing amino acid positions on the translated cds and their corresponding aa and codon residues for each organism (note: human is denoted by targprotid)
		
		i.e.
		index_on_tblastn_translated_alignment (sseq in tblastn results): {targprot:(aa_res, codon_for_res_pos), baboon:(aa_res, codon_for_res_pos), ...}
	"""

	#print(translated_aa_tblastn_alignment_with_orig_prot)
	dict_of_translated_cds_pos_and_aa_codon = {}
	ind_translated_aa = 0


	#Run through each residue position in the msa alignment --corresponds to positions in the dnds results
	# Match the residues to the original translated aa aligned seq from tblastn and store their indices and dnds results
	# Note there may be gaps in the msa that don't exist in the original translated aa aligned seq, so we skip those since they won't help with mapping back to the original targ prot seq
	for msa_alignment_pos in codon_mappings_dict[targ_prot]:
		if ind_translated_aa < len(translated_aa_tblastn_alignment_with_orig_prot): #if we're past the last residue in the tblastn alignment, then that means everything else is a gap in the msa
			(aa_msa_residue, dna_msa_codon) = codon_mappings_dict[targ_prot][msa_alignment_pos][targ_prot] #get the aa and corresponding codon at the given msa position
			if translated_aa_tblastn_alignment_with_orig_prot[ind_translated_aa] == aa_msa_residue: #residue in translated aa seq alignment corresponds to same residue in msa alignment pos
				
				# store the aa and codons at the given position for each organism
				dict_of_translated_cds_pos_and_aa_codon[ind_translated_aa+1] = codon_mappings_dict[targ_prot][msa_alignment_pos] # store {targprot:(aa_res, codon_for_res_pos), baboon:(aa_res, codon_for_res_pos), ...}
				
				# if msa_alignment_residue == '-': #these are the positions that were already gaps on the tblastn alignment (need to keep these b/c correspond to residues on the alignment of the orig targ protein)
				# 	print((ind_translated_aa+ 1, dnds_val, msa_alignment_residue, msa_alignment_pos))
				ind_translated_aa += 1
				
			else: #this is a gap in the msa alignment(doesn't correspond to the human prot)
				continue

	return dict_of_translated_cds_pos_and_aa_codon



def map_dnds_to_targ_prot_residue_positions(dict_of_alignments, aa_msas, codon_mappings_dict):
	''' Map codons (originally in terms of msa positions) back onto original target seq (so that we now have everything in terms of positions on the target prot)
	Will need to do mappings from msa positions to tblastn subject positions (i.e. translated cds aa seq) then to the positions on the original target protein
	For each targ protein (outer key) store the aligned position indices on the original targ prot seq (key) and aa and codon as a tuple (value)

	Args:
		dict_of_alignments: dictionary storing tblastn alignments
		aa_msas: folder storing amino acid msas
		codon_mappings_dict: dictionary containing msa positions and aa and codon mappings

	Return: 
		dict_of_targ_prot_aa_codons: Dictionary where each targ protein (outer key) stores the aligned position indices on the original targ prot seq (next key) 
		and organism (targprot for human or organism name) as most inner key and aa and codon (tuple) (value)
		e.g. {targ_prot1: {pos1: {human_targprot_id:(aa_res, codon_for_res_pos), baboon:(aa_res, codon_for_res_pos), ...}, pos2: {..}, ...}, targ_prot2: ...}
	'''

	
	print('Mapping msa aas and codons onto the original target protein seq . . .')
	dict_of_targ_prot_aa_codons = {}

	for targ_prot in dict_of_alignments:

		dict_of_targ_prot_aa_codons[targ_prot] = {}


		#we want to keep indices in terms of targ seq
		
		pos = dict_of_alignments[targ_prot][0]

		targ_prot_aligned  = dict_of_alignments[targ_prot][4]
		translated_aa_aligned = dict_of_alignments[targ_prot][5]

		# parse and get necessary mapping dicts btwn translated aa aligned seq and msa and dn, ds positions and values
		aa_msa_file = osp.join(aa_msas, targ_prot + '_msa.fasta')

		msa_aligned_aa_seq = get_translated_aa_msa_alignment_fasta(aa_msa_file) #get human msa alignment seq (this is the seq from msa alignment btwn the translated cds aa seq with its orthologs)
		
		dict_of_translated_cds_pos_and_aa_codon = map_translated_aa_residues_to_msa_aligned_residues(translated_aa_aligned, codon_mappings_dict, targ_prot, aa_msa_file) # dictionary storing amino acid positions on the translated cds and their correspinding aa and codons

		for targ_prot_res, translated_aa_res, ind_translated_aa in zip(targ_prot_aligned, translated_aa_aligned, dict_of_translated_cds_pos_and_aa_codon):
			
			(aa_msa_residue, dna_msa_codon) = dict_of_translated_cds_pos_and_aa_codon[ind_translated_aa][targ_prot] #for verification purposed
			if targ_prot_res == '-': #gap in targ prot alignment, so skip b/c no residue to store
				continue
			else:
				#targ prot mapped to a gap in translated aa (store '-') or mapped to an amino acid (store the amino acid letter)
				if aa_msa_residue != translated_aa_res: #just to verify that things worked
					print('incorrect')
				dict_of_targ_prot_aa_codons[targ_prot][pos] = dict_of_translated_cds_pos_and_aa_codon[ind_translated_aa]
				pos += 1
		# print(f'Done: {targ_prot}')
	
	return dict_of_targ_prot_aa_codons



def bin_residues(dict_of_targ_prot_aa_codons, dict_of_targprot_residue_annos):
	''' Save all residues that belong in a residue category

	go through dict_of_targ_prot_aa_codons and check each residue's category annotation using dict_of_targprot_residue_annos and store under appropriate dictionary key
	
	Args:
		dict_of_targ_prot_aa_codons: Dictionary where each targ protein (outer key) stores the aligned position indices on the original targ prot seq (key) and aa and codon (tuple) (value)
		dict_of_targprot_residue_annos: dict storing the target proteins, their residues (those that map to PDB), and residue annotation e.g. {targprot1: {res1: 'surface', res2: 'surface', res3: 'exo', ...}, targprot2: {...}, ...}
	Return:
	{'exo': {human:[codon1_info, codon2_info, ...], cow:[...], ...,} 'endo: {...}, ...'}, where codon*info is a triple with (targ_protid, organism codon, position_on_targ_prot)
	'''
	dict_of_all_residues_codons_in_each_category = {'exo': {}, 'endo': {}, 'mimicry': {}, 'surface': {}, 'buried': {}, 'all_exo': {}, 'all_endo': {}}
	for targ_prot in dict_of_targ_prot_aa_codons:
		for pos in dict_of_targ_prot_aa_codons[targ_prot]:
			#print(pos)
			#print(dict_of_targprot_residue_annos)
			if pos in dict_of_targprot_residue_annos[targ_prot]: #this position mapped to a pdb chain (so we know whether its one of the 3 interfacial residues or surface/buried)
				residue_annos = dict_of_targprot_residue_annos[targ_prot][pos]

				codons_to_store = dict_of_targ_prot_aa_codons[targ_prot][pos] #this consists of all organism codons at this position of the target protein
				
				for org in codons_to_store:
					
					organism_codon = codons_to_store[org][1] #second value in tuple is the codon

					if org == targ_prot:
						org_name = 'human'
					else:
						org_name = org

					for residue_anno in residue_annos:
						if org_name not in dict_of_all_residues_codons_in_each_category[residue_anno]:
							dict_of_all_residues_codons_in_each_category[residue_anno][org_name] = []
						dict_of_all_residues_codons_in_each_category[residue_anno][org_name].append((targ_prot, organism_codon, pos))

	return dict_of_all_residues_codons_in_each_category



def output_binned_codons(dict_of_all_residues_codons_in_each_category, output_folder):
	''' Output one fasta file per residue category, each fasta file with have 20 sequences (one per organsism)
	Made from concatenating the residues in that organism that belong to the specified category.

	Args:
		dict_of_all_residues_codons_in_each_category: dictionary storing codons from each organism for each residue category
		output_folder: folder to store the output fasta files

	'''
	
	for residue_anno in dict_of_all_residues_codons_in_each_category:
		with open(osp.join(output_folder, residue_anno + '.binned_codons.fasta'), 'w') as of:
			for org_name in dict_of_all_residues_codons_in_each_category[residue_anno]:
				
				header = org_name
				list_of_codons = [codon_info[1] for codon_info in dict_of_all_residues_codons_in_each_category[residue_anno][org_name]] #get just the codon
				concatenated_codon_seq = ''.join(list_of_codons)

				#replace any ambiguous chars (e.g. 'X') with a '?' (since only '?' is used to represent ambiguous chars in paml)
				concatenated_codon_seq = concatenated_codon_seq.replace('X', '?')

				of.write('>' + org_name + '\n')
				of.write(concatenated_codon_seq + '\n')

				print(f'\tBinned codons fasta written to file for {org_name}')
			
		

	#if targ_prot doesn't have a certain residue type, we skip that targ prot for the given residue type.
def main():
	script_dir = osp.dirname(__file__)

	aa_msas = osp.join(script_dir, '..', 'data', 'dnds', 'aa_msas')
	dna_msas = osp.join(script_dir, '..', 'data', 'dnds', 'dna_msas')

	
	# files needed for parsing tblastn alignments
	human_blast_results_folder = osp.join(script_dir, '..', 'data', 'blast', 'results', 'human_protein_to_cds')
	blast_results_folder = osp.join(script_dir, '..', 'data', 'blast', 'results', 'protein_to_cds')
	human_target_proteins_fasta = osp.join(script_dir, '..', 'data', 'homologs_and_conservation', 'orthologs', 'human_target_proteins.fasta')

	# Files with exo-specific, endo-specific, and mimicry interfacial residues as well as surface and buried residues
	folder_name = 'v_target_h'
	exo_specific_file =  osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_specific_interfacial_residues.tsv")
	endo_specific_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_specific_interfacial_residues.tsv")
	mimicry_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "mimicked_interfacial_residues.tsv")
	surface_res_file = osp.join(script_dir, '..', 'data', "surface_residues", folder_name, "valid_surface_residues_on_uniprot_seq.tsv")
	buried_res_file = osp.join(script_dir, '..', 'data', "buried_residues", folder_name, "valid_buried_residues_on_uniprot_seq.tsv")

	all_exo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "exogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")
	all_endo_file = osp.join(script_dir, '..', 'data', "categorized_interfacial_residues", folder_name, "endogenous_uniprot_interfacial_residues_grouped_by_uniprot_id.tsv")

	output_folder = osp.join(script_dir, '..', 'data', 'binned_dnds')

	if not osp.exists(output_folder):
		os.makedirs(output_folder)
	
	#Load target proteins
	human_target_proteins = parse_fasta_file_into_dict(human_target_proteins_fasta, change_header=False, split_header_by='|', ind_to_keep = 1)

	# Get dictionary of tblastn alignments for human (this includes ensembl and refseq tblastn results)	
	dict_of_alignments = parse_human_tblastn_results(human_blast_results_folder, blast_results_folder, human_target_proteins)

	# Get mapping dict between msa positions and corresponding aa residues and codons
	print('Making mapping dict between msa positions and corresponding aas and codons . . .')
	codon_mappings_dict = get_codon_mappings(aa_msas, dna_msas)


	# Get mapping between targ prot positions and corresponding msa aa and codons in each organism
	dict_of_targ_prot_aa_codons = map_dnds_to_targ_prot_residue_positions(dict_of_alignments, aa_msas, codon_mappings_dict)


	# Get residue category annotations for target protein residues
	print('Making dict to store residue category annotations . . .')
	dict_of_targprot_residue_annos = annotate_residues_with_category(exo_specific_file, endo_specific_file, mimicry_file, surface_res_file, buried_res_file, all_exo_file, all_endo_file, dna_msas)


	# Binning residues into the 5 residue categories (exo, endo, mimicry, surface, buried, all_exo, all_endo)
	print('Binning residues . . .')
	dict_of_all_residues_codons_in_each_category = bin_residues(dict_of_targ_prot_aa_codons, dict_of_targprot_residue_annos)


	# Output binned codons into fasta files (20 organsism codon sequences in each file; one fasta file per residue type)
	print('Outputing binned residues into their respective fasta files . . .')
	output_binned_codons(dict_of_all_residues_codons_in_each_category, output_folder)
if __name__ == '__main__':
	main()