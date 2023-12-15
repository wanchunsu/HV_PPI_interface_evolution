#!/bin/bash
pdb_seqs="../data/PDB/human_only_pdb_seqres.txt"
path_to_uniprot_fasta="../data/pdb_annotations/fasta_files/human_reviewed_uniprot.formatted_header.fasta"

results_dir="../data/blast/results"



echo "Making uniprot seqs database ..."

#make database for human uniprot sequences
makeblastdb -in $path_to_uniprot_fasta -dbtype prot -out ../data/PDB/human_db

echo "Blasting PDB human seqs against human uniprot database ..."

#Blast human PDB chains against human uniprot database
blastp -db ../data/PDB/human_db -query $pdb_seqs -outfmt 6 -out "$results_dir/human_pdb_against_human_uniprot.blast_out.txt" -max_hsps 1 -max_target_seqs 3 #just examine best 2 targets