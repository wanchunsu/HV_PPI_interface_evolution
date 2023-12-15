'''

Convert interfacial residues ('auth') into 'label' residues
(this is for mapping between UniProt-to-PDB BLAST alignments)
============================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os.path as osp
import os


def get_label_auth_conversion_dict(pdb_files_folder, pdb_id):
	""" Making a dict to store the auth residue numbers and their corresponding label residue numbers for each chain in the PDB structure of interest
		Args:
			pdb_files_folder: folder with all pdb files
			pdb_id: pdb structure of interest
	"""
	pdb_file = osp.join(pdb_files_folder, pdb_id + '.cif')
	label_auth_mapping_dict = {} # {chainid: {auth_id: label_id, auth_id2: label_id2, ...}}

	list_of_atom_attr = []
	getting_attr = False
	got_attrib = False
	label_seq_index = -1
	auth_seq_index= -1
	auth_asym_index = ''
	res_num = 0
	chain_curr = ''
	with open(pdb_file) as pdb_f:
		for line in pdb_f:

			if line.startswith("_atom_site."):
				if getting_attr == False:
					getting_attr = True
				list_of_atom_attr.append(line.strip())
			elif not line.startswith("_atom_site.") and getting_attr==True: #attributes already stored
				
				label_seq_index = list_of_atom_attr.index('_atom_site.label_seq_id')
				auth_seq_index = list_of_atom_attr.index('_atom_site.auth_seq_id')
				auth_asym_index = list_of_atom_attr.index('_atom_site.auth_asym_id')
				getting_attr = False
				got_attrib = True
				
			if line.startswith("ATOM") and not line.startswith("ATOMS") and got_attrib == True:
				atom_info = line.split()
				if len(atom_info) != len(list_of_atom_attr):
					print(line)
				label_residue = atom_info[label_seq_index]
				if res_num != label_residue or (res_num == label_residue and atom_info[auth_asym_index] != chain_curr): #first atom of the residue (or in cases where prev. chain's last res has same res num as curr chain resnum, we check that the prev chain was different)

					auth_residue = atom_info[auth_seq_index]
					auth_chain_id = atom_info[auth_asym_index]

					if auth_chain_id not in label_auth_mapping_dict:
						label_auth_mapping_dict[auth_chain_id] = {}

					label_auth_mapping_dict[auth_chain_id][auth_residue] = label_residue

					res_num = label_residue
					chain_curr = auth_chain_id
			elif line.startswith('_pdbx_poly_seq_scheme.') or line.startswith('_atom_site_anisotrop.'): # we've gone through all the ATOMs
				break

	return label_auth_mapping_dict

def make_interfacial_residues_label_id_dict(interfacial_residues, int_res_with_pdb_label_residues, pdb_files_folder):
	""" Remake the interfacial residues file where the interfacial residues (auth ids) are now the label ids (this is for better mapping btwn pdb seqres and blast results)

	Args:
		interfacial residues: original interfacial residues file with the auth ids
		int_res_with_label_residues: new interfacial residues file (with label ids) to be outputted
		pdb_files_folder: folder storing all the pdb files

	"""
	with open(interfacial_residues) as int_res:
		header = next(int_res)
		with open(int_res_with_pdb_label_residues, 'w') as int_res_with_label_residues:
			int_res_with_label_residues.write(header)
			for line in int_res:
				line = line.strip("\n")
				info = line.split("\t")
				targ_prot =  info[0]
				pdb_chain = info[1]
				pdb_id = pdb_chain.split('_')[0]
				chain_id = pdb_chain.split('_')[1]
				int_res_auth = info[2].split(",")

				label_auth_mapping_dict = get_label_auth_conversion_dict(pdb_files_folder, pdb_id)
				int_res_label = [label_auth_mapping_dict[chain_id][i] for i in int_res_auth]
				int_res_with_label_residues.write(targ_prot + "\t" + pdb_chain +"\t" + ",".join(int_res_label)+"\n")

			
		


def main():
	script_dir = osp.dirname(__file__)
	pdb_files_folder = osp.join(script_dir, '..', 'data', 'pdb_files')

	endo_interfacial_residues = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'endogenous_interfacial_residues.tsv')
	exo_interfacial_residues= osp.join(script_dir, '..', 'data', 'exo_and_endo', 'exogenous_interfacial_residues.tsv')

	#since our current interfacial residues are 'auth' residues and pdb seqres (and blast results) are in 'label' residues, we need to convert 
	#make a new file with the corresponding 'label' residues so that we can map to blast results easily
	endo_int_res_with_pdb_label_residues = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'endogenous_interfacial_residues_with_pdb_label_residues.tsv')
	exo_int_res_with_pdb_label_residues = osp.join(script_dir, '..', 'data', 'exo_and_endo', 'exogenous_interfacial_residues_with_pdb_label_residues.tsv')
	print("Getting the label interfacial residues for endogenous ppis . . . ")
	make_interfacial_residues_label_id_dict(endo_interfacial_residues, endo_int_res_with_pdb_label_residues, pdb_files_folder)

	print("Getting the label interfacial residues for exogenous ppis . . . ")
	make_interfacial_residues_label_id_dict(exo_interfacial_residues, exo_int_res_with_pdb_label_residues, pdb_files_folder)



if __name__ == '__main__':
	main()