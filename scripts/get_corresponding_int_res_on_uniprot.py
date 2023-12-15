'''

Get interfacial residues in terms of UniProt protein residue positions.
=======================================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
import os
from format_int_res_pdb_to_uniprot_blast_results import map_pdb_residues_to_uniprot
from get_interfacial_residues_based_on_pdb_label_res import get_label_auth_conversion_dict
from get_interfacial_residues_using_sasa import get_stats_for_interface_residues_exogenous, get_stats_for_interface_residues_endogenous



def get_corresponding_int_res_on_uniprot(interfacial_residues_with_pdb_label_residues, blast_results_dict, output_interfacial_residues_in_terms_of_uniprot):

	"""
	Go through the interfacial residues with pdb labels and for each of them, get the interfacial residues based on the uniprot sequence 
	if a certain interfacial residue is not located on the uniprot seq, or it maps to a gap on the uniprot seq, we ignore it

	Args:
		interfacial_residues_with_pdb_label_residues: file with interfacial residues in terms of the label residues (mappable to blast results)
		blast_results_dict: dictionary storing the mapping btwn pdb label residues and uniprot sequence residues
		output_interfacial_residues_in_terms_of_uniprot: output file with interfacial residues in terms of uniprot sequence
	"""

	with open(interfacial_residues_with_pdb_label_residues) as int_res_pdb:
		header = next(int_res_pdb)
		with open(output_interfacial_residues_in_terms_of_uniprot, 'w') as o:

			# write the headers (columns after interfacial residues on uniprot seq are mainly for verification purposes)
			if "exo" in header: #exogenous
				o.write("target_prot" +"\t" + "pdb_chain" + "\t" + "exo_interfacial_residues_on_uniprot_seq"+ "\t" + "mapped_pdb_label_int_res" + "\t" + "unmapped_pdb_label_int_res"+ "\t" + "pdb_label_int_res_mapped_to_gap" + "\n")
			else: #endogenous
				o.write("target_prot" +"\t" + "pdb_chain" + "\t" + "endo_interfacial_residues_on_uniprot_seq"+ "\t" + "mapped_pdb_label_int_res" + "\t" +"unmapped_pdb_label_int_res"+ "\t" + "pdb_label_int_res_mapped_to_gap" + "\n")
			
			for line in int_res_pdb:#P50750	4ec9_A	10,65,68,13,57,9,59,8,12,62,69,87,100,72,61,85,58,11,60,14
				line = line.strip("\n")

				info = line.split("\t")
				targ_prot = info[0]
				pdb_chain = info[1]
				int_res_pdb = info[2].split(",")
				int_res_pdb = [int(r) for r in int_res_pdb] #convert to int

				if targ_prot not in blast_results_dict:
					# only has one pdb chain and this pdb chain doesn't map to the uniprot significantly (so no blast results whatsoever)
					#  e.g. Q7Z6E9--6e5x_B
					print(f"   target prot: {targ_prot} not in blast_results dict")
					continue
				if pdb_chain not in blast_results_dict[targ_prot]: #this is the case where there is no significant blast result between uniprot seq and pdb chain 
					#(e.g. Q8WUM4--2zne_C, Q8WUM4--2zne_D, P04637--4fz3_B, )
					print(f"   For {targ_prot}, the pdb chain {pdb_chain} is not in the blast_results_dict")
					continue
				int_res_mapped = []
				int_res_on_uniprot = []
				int_res_not_mapped_to_uniprot = []
				int_res_mapped_to_gap = []
				for res in int_res_pdb:
					
					if res not in blast_results_dict[targ_prot][pdb_chain]: #this interfacial residue doesn't map to an aligned region on the uniprot seq
						int_res_not_mapped_to_uniprot.append(str(res))
						continue
					elif blast_results_dict[targ_prot][pdb_chain][res] == -1: #this interfacial residue maps to a gap on the uniprot seq
						int_res_mapped_to_gap.append(str(res))
					elif blast_results_dict[targ_prot][pdb_chain][res] != -1: #mapped to residue on the uniprot seq
						mapped_res_on_uniprot = blast_results_dict[targ_prot][pdb_chain][res]
						int_res_mapped.append(str(res))
						int_res_on_uniprot.append(str(mapped_res_on_uniprot))
				
				o.write(targ_prot + "\t" + pdb_chain +"\t" + ",".join(int_res_on_uniprot) + "\t" + ",".join(int_res_mapped) + "\t" 
					+ ",".join(int_res_not_mapped_to_uniprot) +"\t" + ",".join(int_res_mapped_to_gap) + "\n")


def group_interfacial_residues_by_uniprot_id(int_res_on_uniprot_seq, int_res_grouped_by_uniprot_seq):
	dict_of_int_res_on_uniprot = {}
	with open(int_res_on_uniprot_seq) as int_res:
		header = next(int_res)
		
		for line in int_res:
			line = line.strip("\n")
			info = line.split("\t")
			uniprot_id = info[0]
			#pdb_chain = info[1]
			if info[2] == '':
				print(line)
				continue
			int_res_on_uniprot = info[2].split(",")

			int_res_on_uniprot = [int(r) for r in int_res_on_uniprot]

			if uniprot_id not in dict_of_int_res_on_uniprot:
				dict_of_int_res_on_uniprot[uniprot_id] = []

			dict_of_int_res_on_uniprot[uniprot_id].extend(int_res_on_uniprot)
			sorted_and_unique = sorted(list(set(dict_of_int_res_on_uniprot[uniprot_id])))
			dict_of_int_res_on_uniprot[uniprot_id] = sorted_and_unique

	with open(int_res_grouped_by_uniprot_seq, 'w') as o:
		if "exo" in header:
			o.write("target_prot" +"\t" + "all_exo_interfacial_residues_on_uniprot_seq" + "\n")
		else:
			o.write("target_prot" +"\t" + "all_endo_interfacial_residues_on_uniprot_seq" + "\n")

		for uniprotid in dict_of_int_res_on_uniprot:
			int_res_converted_to_str = [str(i) for i in dict_of_int_res_on_uniprot[uniprotid]]

			o.write(uniprotid +"\t" + ",".join(int_res_converted_to_str) + "\n")


def remake_ppi_files_with_mapped_int_res_to_uniprot_seq(original_int_res, blast_results_dict, new_ppi_file_with_uniprot_int_res, pdb_files_folder, exo_int_res_grouped_by_uniprot_seq, endo_or_exo):
	"""
	"""
	#blast results_dict: {uniprot: {pdb_chain: pdb_res_id: corresponding uniprot_res_id}}
	#P50750	O60563	4ec9	A	B	12,9,11,8,7,57,58,64,59,60,61,56,10,13,86,99,71,67,84,68
	#P50750	4ec9_A	9,64,67,12,56,8,58,7,11,61,68,86,99,71,60,84,57,10,59,13
	dict_of_mapped_uniprot = {}
	with open(original_int_res) as orig:
		header = next(orig)
		with open(new_ppi_file_with_uniprot_int_res, 'w') as new_fi:
			if endo_or_exo == 'exo': #exogenous
				new_fi.write("target_prot" + "\t" + "partner_prot" + "\t" + "PDBid" + "\t" + 
					"target_chainid" + "\t" + "partner_chainid" + "\t" + "exo_interface_residues_on_target_uniprot_seq" +"\n")
			else: #endogenous
				new_fi.write("target_prot" + "\t" + "partner_prot" + "\t" + "PDBid" + "\t" + 
					"target_chainid" + "\t" + "partner_chainid" + "\t" + "endo_interface_residues_on_target_uniprot_seq" +"\n")

				#need to check if a target protein in endo has any exogenous int res that map to uniprot (otherwise pointless to keep it because our main focus is on exogenous ppis)
				list_of_exo_targ_prots_with_mapped_int_res= []
				with open(exo_int_res_grouped_by_uniprot_seq) as grouped_int_res:
					next(grouped_int_res)
					for line in grouped_int_res:
						targ_prot_exo = line.split("\t")[0]
						list_of_exo_targ_prots_with_mapped_int_res.append(targ_prot_exo)

			for line in orig:
				line = line.strip("\n")
				info = line.split("\t")
				targ_prot = info[0]

				#check only for endo:
				if endo_or_exo == 'endo': #endogenous file:
					if targ_prot not in list_of_exo_targ_prots_with_mapped_int_res: #no point in looking at a target protein that doesn't have any mapped exo interface residues on its uniprot seq
						print(f'targ prot in endo without exo results: {targ_prot}')
						continue
				partner = info[1]
				pdb = info[2]
				targ_chain = info[3]
				partner_chain = info[4]
				auth_int_res = info[5].split(",")
				targ_pdb_chain = pdb+'_'+targ_chain

				if targ_prot not in blast_results_dict: # only has one pdb chain and this pdb chain doesn't map to the uniprot significantly (so no blast results whatsoever)
					continue
				if targ_pdb_chain not in blast_results_dict[targ_prot]:#this is the case where there is no significant blast result between uniprot seq and pdb chain 
					continue

				auth_to_label_pdb_res_dict = get_label_auth_conversion_dict(pdb_files_folder, pdb) #{chainid: {auth_id: label_id, auth_id2: label_id2, ...}}
				auth_to_label_pdb_res_dict_chain = auth_to_label_pdb_res_dict[targ_chain] #mapping btwn auth_id and label_id for residues
				
				
				label_pdb_res_to_uniprot_mapping = blast_results_dict[targ_prot][targ_pdb_chain] #mapping between pdb_label_id and correspinding uniprot residue id

				#convert: auth_id --> label_id --> uniprot seq id
				converted_int_res_on_uniprot_seq = []
				for res in auth_int_res:
					label_id = auth_to_label_pdb_res_dict_chain[res] #get corresponding label id

					#get corresponding uniprot seq res id (if exists), else ignore this interfacial residue

					#mapping between pdb_label_id and uniprot residue id exists and is not a gap
					if int(label_id) in label_pdb_res_to_uniprot_mapping:
						uniprot_residue_id = label_pdb_res_to_uniprot_mapping[int(label_id)]
						if uniprot_residue_id != -1: 
							converted_int_res_on_uniprot_seq.append(uniprot_residue_id)
				if converted_int_res_on_uniprot_seq !=[]:
					converted_int_res_on_uniprot_seq = [str(i) for i in converted_int_res_on_uniprot_seq]
					new_fi.write(targ_prot + "\t" + partner + "\t" + pdb + "\t" + targ_chain + "\t" + partner_chain + "\t" + ",".join(converted_int_res_on_uniprot_seq) + "\n")


def get_ppi_stats(new_endo_ppi_file_with_uniprot_int_res, new_exo_ppi_file_with_uniprot_int_res):
	print("Endogenous:")
	get_stats_for_interface_residues_endogenous(new_endo_ppi_file_with_uniprot_int_res)
	print("______________________________")
	print("Exogenous:")
	get_stats_for_interface_residues_exogenous(new_exo_ppi_file_with_uniprot_int_res)

		
def main():
	script_dir = osp.dirname(__file__)
	pdb_files_folder = osp.join(script_dir, '..', 'data', 'pdb_files')

	
	blast_result_folder = osp.join(script_dir, '..', 'data', 'blast', 'results', 'targ_prot_against_pdb_chains')

	endo_int_res_with_pdb_label_residues = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'endogenous_interfacial_residues_with_pdb_label_residues.tsv')
	exo_int_res_with_pdb_label_residues = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'exogenous_interfacial_residues_with_pdb_label_residues.tsv')
	
	endo_int_res_on_uniprot_seq = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'endogenous_interfacial_residues_on_uniprot_seq.tsv')
	exo_int_res_on_uniprot_seq = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'exogenous_interfacial_residues_on_uniprot_seq.tsv')

	endo_int_res_on_uniprot_seq_grouped_by_uniprot_id = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'endogenous_interfacial_residues_on_uniprot_seq_grouped_by_uniprot_id.tsv')
	exo_int_res_on_uniprot_seq_grouped_by_uniprot_id = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'exogenous_interfacial_residues_on_uniprot_seq_grouped_by_uniprot_id.tsv')

	original_endo_results = osp.join(script_dir, '..', 'data', 'endo_with_verified_exo_int_res.tsv')
	original_exo_results = osp.join(script_dir, '..', 'data', 'exo_verified_pathogenic_int_res.tsv') 
	
	new_endo_ppi_file_with_uniprot_int_res = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'endogenous_ppis_with_int_res_on_uniprot_seq.tsv')
	new_exo_ppi_file_with_uniprot_int_res = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'exogenous_ppis_with_int_res_on_uniprot_seq.tsv')

	
	blast_results_dict = map_pdb_residues_to_uniprot(blast_result_folder)
	
	print("Done making pdb label res to uniprot res dictionary.")

	print("Getting corresponding endogenous interfacial residues on the uniprot sequence . . .")
	get_corresponding_int_res_on_uniprot(endo_int_res_with_pdb_label_residues, blast_results_dict, endo_int_res_on_uniprot_seq)
	print("Getting corresponding exogenous interfacial residues on the uniprot sequence . . .")
	get_corresponding_int_res_on_uniprot(exo_int_res_with_pdb_label_residues, blast_results_dict, exo_int_res_on_uniprot_seq)
	
	print("Grouping the interfacial residues by uniprot id and outputting together")
	group_interfacial_residues_by_uniprot_id(endo_int_res_on_uniprot_seq, endo_int_res_on_uniprot_seq_grouped_by_uniprot_id)
	group_interfacial_residues_by_uniprot_id(exo_int_res_on_uniprot_seq, exo_int_res_on_uniprot_seq_grouped_by_uniprot_id)

	print("Making new ppi files with uniprot seq interfacial residues . . .")
	remake_ppi_files_with_mapped_int_res_to_uniprot_seq(original_endo_results, blast_results_dict, new_endo_ppi_file_with_uniprot_int_res, pdb_files_folder,exo_int_res_on_uniprot_seq_grouped_by_uniprot_id, 'endo')
	remake_ppi_files_with_mapped_int_res_to_uniprot_seq(original_exo_results, blast_results_dict, new_exo_ppi_file_with_uniprot_int_res, pdb_files_folder, exo_int_res_on_uniprot_seq_grouped_by_uniprot_id, 'exo')
	

	#Printing stats 
	print("########################       Printing stats       ########################\n")
	
	print(f" Stats for ppis with interfacial residues that map onto the uniprot seq \n")
	get_ppi_stats(new_endo_ppi_file_with_uniprot_int_res, new_exo_ppi_file_with_uniprot_int_res)
	
	
if __name__ == '__main__':
	main()