'''

Annotate exogenous PPIs based on human target protein and PDB annotations
=========================================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp

def load_data(input_file_with_anno):
	dict_to_store_data = {}
	with open(input_file_with_anno) as annos:
		for line in annos:
			line = line.strip("\n")
			id_ = line.split("\t")[0]
			annotation = line.split("\t")[1]

			dict_to_store_data[id_] = annotation
	return dict_to_store_data

def annotate_exo_ppis(all_exo_ppis, annotated_human_targs_file, annotated_pdbs_file, output_annotated_ppis):
	""" Annotating exogenous ppis based on the annotations for the human target proteins and their pdbs

	Args:
		all_exo_ppis: file with all exogenous ppis
		annotated_human_targs_file: file with exo human targets and their antibody/mhc status (T = antibody/mhc protein, F = not antibody/mhc protein), separated by a tab char
		annotated_pdbs_file: file with exo pdbs and their antibody/mhc status separated by a tab char
			For the above file each pdb annotation represents the following:
			T = all ppis mapping to this pdb are human protein (antibody/mhc protein) targetting virus protein, 
			F = all ppis mapping to this pdb are virus targetting human protein (could potentially be an antibody/mhc protein being targetted) 
			D = contains more than one type of ppi (both T and F) mapping to this pdb; virus targetting non-antibody/mhc protein, antibody/mhc protein targetting virus
		output_annotated_ppis: exogenous ppis with annotations

	"""
	dict_of_human_targ_annos =  load_data(annotated_human_targs_file)
	dict_of_pdb_annos = load_data(annotated_pdbs_file)

	with open(all_exo_ppis) as exo_ppis, open(output_annotated_ppis, 'w') as o:
		header = next(exo_ppis)
		header = header.strip("\n")

		header = header + "\tantibody_or_mhc\n"
		o.write(header)

		for line in exo_ppis:
			line  = line.strip("\n")
			info = line.split("\t")

			human_prot = info[0]
			pdb_id = info[2]

			if dict_of_human_targ_annos[human_prot] == 'F': # human target protein is not an antibody/mhc protein
				o.write(line + "\t" + "F" + "\n")
			elif dict_of_human_targ_annos[human_prot] == 'T': #human target protein is an antibody/mhc protein
				if dict_of_pdb_annos[pdb_id] == 'T' or dict_of_pdb_annos[pdb_id] == 'D':
					o.write(line + "\t" + "T" + "\n")
				elif dict_of_pdb_annos[pdb_id] == 'F':
					o.write(line + "\t" + "F" + "\n")
					print(f"{pdb_id}: This antibody/mhc protein is targeted by viruses")
				else:
					print(f'invalid annotation for {pdb}')
			else:
				print(f'invalid annotation for {human_prot}')

def count_num_ppis_in_v_targ_h_or_h_targ_v(output_annotated_ppis):
	list_of_pathogenic = []

	list_of_immune = []
	with open(output_annotated_ppis) as fi:
		#target_prot	partner_prot	PDBid	target_chainid	partner_chainid	target_interfacial_residues	antibody_or_mhc
		next(fi)
		for line in fi:
			info = line.strip("\n").split("\t")
			pair = (info[0], info[1])
			anno = info[-1]
			if anno == 'F' and pair not in list_of_pathogenic:
				list_of_pathogenic.append(pair)
			if anno == 'T' and pair not in list_of_immune:
				list_of_immune.append(pair)

	targs_for_patho, parts_for_patho = count_unique_targs_and_partners(list_of_pathogenic)
	targs_for_immune, parts_for_immune = count_unique_targs_and_partners(list_of_immune)

	
	print(f'### V_target_H results: ###')
	print(f'Number of PPIs: {len(list_of_pathogenic)}')
	print(f'Number of target proteins: {len(targs_for_patho)}')
	print(f'Number of partner proteins: {len(parts_for_patho)}')

	print(f'\n### H_target_V results: ###')
	print(f'Number of PPIs: {len(list_of_immune)}')
	print(f'Number of target proteins: {len(targs_for_immune)}')
	print(f'Number of partner proteins: {len(parts_for_immune)}')

def count_unique_targs_and_partners(list_of_pairs):
	list_of_targ_prots = []
	list_of_part_prots = []
	for p in list_of_pairs:
		if p[0] not in list_of_targ_prots:
			list_of_targ_prots.append(p[0])
		if p[1] not in list_of_part_prots:
			list_of_part_prots.append(p[1])

	return list_of_targ_prots, list_of_part_prots
def main():
	script_dir = osp.dirname(__file__)
	all_exo_ppis = osp.join(script_dir, '..', 'data', 'combined_homology_and_pdb_anno_human_virus_interface_residues.tsv')
	annotated_human_targs_file = osp.join(script_dir, '..', 'data', 'antibody_mhc_complex_annotations',  'exogenous_human_targs_with_antibody_mhc_annotations.tsv')
	annotated_pdbs_file = osp.join(script_dir, '..', 'data', 'antibody_mhc_complex_annotations',  'exogenous_ppi_pdbs_with_antibody_mhc_annotations.tsv')
	output_annotated_ppis = osp.join(script_dir, '..', 'data', 'annotated_combined_homology_and_pdb_anno_human_virus_interface_residues.tsv')

	
	annotate_exo_ppis(all_exo_ppis, annotated_human_targs_file, annotated_pdbs_file, output_annotated_ppis)
	count_num_ppis_in_v_targ_h_or_h_targ_v(output_annotated_ppis)
if __name__ == '__main__':
	main()

