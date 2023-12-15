'''

Short script to format fasta headers.
=================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

from pathlib import Path
from Bio import SeqIO
import os
from fasta_tools import format_uniprot_fasta_headers


def main(): 
	script_dir = os.path.dirname(__file__)
	fasta_fi_path = os.path.join(script_dir, '../data/pdb_annotations/fasta_files')
	format_uniprot_fasta_headers(os.path.join(fasta_fi_path, "human_reviewed_uniprot.fasta"), os.path.join(fasta_fi_path, "human_reviewed_uniprot.formatted_header.fasta"), '|', 1, 1)

	format_uniprot_fasta_headers(os.path.join(fasta_fi_path, "virus_affecting_humans_reviewed_uniprot.fasta"), os.path.join(fasta_fi_path, "virus_affecting_humans_reviewed_uniprot.formatted_header.fasta"), '|', 1, 1)

if __name__ == '__main__':
	main()
