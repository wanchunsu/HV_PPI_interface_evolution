#!/bin/bash
human_target_folder_path="../data/exo_and_endo/human_target_fasta_folder"
pdb_fasta_folder="../data/exo_and_endo/pdb_fasta_folder" 

results_dir="../data/blast/results/targ_prot_against_pdb_chains"

mkdir -p $results_dir	

database_dir="../data/exo_and_endo/blast_dbs"
mkdir -p $database_dir

for targ_fasta in $human_target_folder_path/*; do
	targ_fasta_with_ext=${targ_fasta##*/}
	targ_fasta_id="${targ_fasta_with_ext%.*}"
	echo $targ_fasta
	targ_prot_pdbs=$pdb_fasta_folder/$targ_fasta_id.pdbs.fasta
	echo $targ_prot_pdbs

	output_file=$results_dir/$targ_fasta_id".blast.out"
	echo $output_file	
	
	db_name=$database_dir/$targ_fasta_id"_pdbs_db"
	echo $db_name
	makeblastdb -in $targ_prot_pdbs -dbtype prot -out $db_name
	blastp -subject $targ_prot_pdbs -query $targ_fasta -outfmt "6 std qseq sseq" -out $output_file -max_hsps 1 -num_alignments 1000000 -seg no #get only the highest scoring pair for each hit and allow for as many alignments as possible because we want all the subjects to be present
	done