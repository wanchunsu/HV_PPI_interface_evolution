'''

Process human and virus pdb-to-uniprot blast results 
=========================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

import os
import os.path as osp

def get_only_reviewed_uniprot_to_pdb_mapping(pdb_to_uniprot_anno, all_reviewed_uniprot_ids, reviewed_pdb_uniprot_out_fi):
	""" Filter out non-reviewed uniprot ids from the pdb-to-uniprot annotation mapping file
	Args:
		pdb_to_uniprot_anno: original pdbchain to uniprot mapping file
		all_reviewed_uniprot_ids: list of all uniprot ids that are reviewed (swissprot)
		reviewed_pdb_uniprot_out_fi: out file with only pdbchain to uniprot mappings where the uniprot id is reviewed
	"""
	list_all_reviewed_uniprot_ids = []
	with open(all_reviewed_uniprot_ids) as rev:
		for line in rev:
			line = line.strip("\n")
			list_all_reviewed_uniprot_ids.append(line)
	
	line_num = 0

	with open(pdb_to_uniprot_anno, 'r') as mapping_fi:
		with open(reviewed_pdb_uniprot_out_fi, 'w') as outfi:
			
			for line in mapping_fi:
				if line_num < 2:
					outfi.write(line)
					line_num +=1
					continue
				uniprot_id = line.split("\t")[2]
				if uniprot_id in list_all_reviewed_uniprot_ids: #this protein id is reviewed, so we keep the entry
					outfi.write(line)
					outfi.flush()
				line_num +=1

def process_uniprot_to_pdb_blast(blast_results, blast_best_hits):
	"""
		1. qseqid: query or source (e.g., gene) sequence id
		2. sseqid: subject or target (e.g., reference genome) sequence id
		3. pident: percentage of identical matches
		4. length: alignment length (sequence overlap)
		5. mismatch: number of mismatches
		6. gapopen: number of gap openings
		7. qstart: start of alignment in query
		8. qend: end of alignment in query
		9. sstart: start of alignment in subject
		10. send: end of alignment in subject
		11. evalue: expect value
		12. bitscore: bit score

	"""
	dict_of_blast_results = {}
	with open(blast_results) as b:
		for line in b:
			info = line.split("\t")
			pdb_chain = info[0]
			uniprot_id = info[1]
			coverage = info[2]
			qstart = info[6]
			qend = info[7]
			sstart = info[8]
			send = info[9]
			e_val = info[10]
			if pdb_chain not in dict_of_blast_results:
				dict_of_blast_results[pdb_chain] = [uniprot_id, coverage, qstart, qend, sstart, send, e_val] #only take the first one
			
	with open(blast_best_hits, 'w') as best_hits_fi:
		best_hits_fi.write("pdb\tchain\tuniprotid\tpdb_start\tpdb_end\tuniprot_start\tuniprot_end\tcoverage\te_val\n")
		for pdbchain in dict_of_blast_results:
			hit = dict_of_blast_results[pdbchain]
			uniprotid = hit[0]
			coverage = hit[1]
			pdb_start = hit[2]
			pdb_end = hit[3]
			uniprot_start = hit[4]
			uniprot_end = hit[5]
			e_val = hit[6]
			best_hits_fi.write(pdbchain.split("_")[0] + "\t" + pdbchain.split("_")[1] 
							+ "\t" + uniprotid + "\t" + pdb_start + "\t" + pdb_end + "\t" 
							+ uniprot_start + "\t" + uniprot_end + "\t" + coverage + "\t" + e_val + "\n")
	
def main():
	script_dir = osp.dirname(__file__)
	pdb_to_uniprot_anno = osp.join(script_dir, '..', 'data', 'pdb_annotations', 'pdb_chain_uniprot.tsv')
	all_reviewed_uniprot_ids = osp.join(script_dir, '..', 'data', 'pdb_annotations', 'all_uniprot_reviewed_ids.txt')
	reviewed_pdb_uniprot_out_fi =  osp.join(script_dir, '..', 'data', 'pdb_annotations', 'pdb_chain_uniprot_reviewed_only.tsv')
	
	get_only_reviewed_uniprot_to_pdb_mapping(pdb_to_uniprot_anno, all_reviewed_uniprot_ids, reviewed_pdb_uniprot_out_fi)

	pdb_to_uniprot_blast_human = osp.join(script_dir, '..', 'data', 'blast', 'results', 'human_pdb_against_human_uniprot.blast_out.1e-10.txt')
	pdb_to_uniprot_blast_virus = osp.join(script_dir, '..', 'data', 'blast', 'results', 'virus_pdb_against_virus_uniprot.blast_out.1e-5.txt')

	best_hits_pdb_to_uniprot_blast_human = osp.join(script_dir, '..', 'data', 'pdb_annotations', 'human_pdb_chain_uniprot.blast_out.tsv')
	best_hits_pdb_to_uniprot_blast_virus = osp.join(script_dir, '..', 'data', 'pdb_annotations', 'virus_pdb_chain_uniprot.blast_out.tsv')

	process_uniprot_to_pdb_blast(pdb_to_uniprot_blast_human, best_hits_pdb_to_uniprot_blast_human)
	process_uniprot_to_pdb_blast(pdb_to_uniprot_blast_virus, best_hits_pdb_to_uniprot_blast_virus)
if __name__ == '__main__':
	main()

