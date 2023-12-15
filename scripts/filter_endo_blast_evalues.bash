#!/bin/bash
#filter off hits with e-values that are not within threshold

path_to_blast_results='../data/blast/results'

awk -F"\t" '$11<=10^-10' "$path_to_blast_results/human_endo_prots_against_pdb.blast_out.txt" > "$path_to_blast_results/human_endo_prots_against_pdb.blast_out.1e-10.txt"

