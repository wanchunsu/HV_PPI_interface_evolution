#!/bin/bash
pdb_seqs="../data/PDB/virus_only_pdb_seqres.txt"
path_to_uniprot_fasta="../data/pdb_annotations/fasta_files/virus_affecting_humans_reviewed_uniprot.formatted_header.fasta"

results_dir="../data/blast/results"



echo "Making uniprot seqs database ..."

#make database for virus uniprot sequences
makeblastdb -in $path_to_uniprot_fasta -dbtype prot -out ../data/PDB/virus_affecting_humans_db

echo "Blasting PDB virus seqs against virus uniprot database ..."

#blast virus PDB chains against human uniprot database
blastp -db ../data/PDB/virus_affecting_humans_db -query $pdb_seqs -outfmt 6 -out "$results_dir/virus_pdb_against_virus_uniprot.blast_out.txt" -max_hsps 1 -max_target_seqs 3 #just examine best 2 targets

