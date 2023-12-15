#!/bin/bash
pdb_seqs="../data/PDB/pdb_seqres.txt"
path_to_uniprot_fasta="../data/IntAct_human_virus/fasta_files"

human_uniprot_fa="$path_to_uniprot_fasta/human_uniprot_sequences.fasta"

virus_uniprot_fa="$path_to_uniprot_fasta/virus_uniprot_sequences.fasta"

results_dir="../data/blast/results"
mkdir -p $results_dir


echo "Making PDB database ..."

#make database for PDB structure sequences
makeblastdb -in "$pdb_seqs" -dbtype prot -out "../data/PDB/pdb_seqres_db"

echo "Blasting human seqs against PDB database ..."

# #blast human uniprot sequences against PDB database
blastp -db "../data/PDB/pdb_seqres_db" -query $human_uniprot_fa -outfmt 6 -out "$results_dir/human_against_pdb.blast_out.txt"

echo "Blasting virus seqs against PDB database ..."

# #blast virus uniprot sequences against PDB database
blastp -db "../data/PDB/pdb_seqres_db" -query $virus_uniprot_fa -outfmt 6 -out "$results_dir/virus_against_pdb.blast_out.txt"