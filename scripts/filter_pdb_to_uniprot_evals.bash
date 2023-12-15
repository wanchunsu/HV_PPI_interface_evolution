#!/bin/bash
#filter off e-values that are not within threshold for pdb to reviewed uniprot mapping

path_to_blast_results='../data/blast/results'

awk -F"\t" '$11<=10^-10' "$path_to_blast_results/human_pdb_against_human_uniprot.blast_out.txt" > "$path_to_blast_results/human_pdb_against_human_uniprot.blast_out.1e-10.txt"

awk -F"\t" '$11<=10^-5' "$path_to_blast_results/virus_pdb_against_virus_uniprot.blast_out.txt" > "$path_to_blast_results/virus_pdb_against_virus_uniprot.blast_out.1e-5.txt"
