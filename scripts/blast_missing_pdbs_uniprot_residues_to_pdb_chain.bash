#!/bin/bash
human_target_folder_path="../data/exo_and_endo/human_target_fasta_folder"
pdb_fasta_folder="../data/exo_and_endo/pdb_fasta_folder" 

results_dir="../data/blast/results/targ_prot_against_pdb_chains"

mkdir -p $results_dir	

database_dir="../data/exo_and_endo/blast_dbs"
mkdir -p $database_dir

for missing_pdb in $pdb_fasta_folder/*.missing_pdbs.fasta; do
	targ_fasta_with_ext=${missing_pdb##*/}
	targ_fasta_id="${targ_fasta_with_ext%.*}"
	targ_fasta_id="${targ_fasta_id%.*}"
	echo $missing_pdb
	targ_prot_fasta=$human_target_folder_path/$targ_fasta_id.fasta
	echo $targ_prot_fasta

	output_file=$results_dir/$targ_fasta_id".missing_pdb.blast.out"
	echo $output_file	
	
	db_name=$database_dir/$targ_fasta_id"_missing_pdbs_db"
	echo $db_name
	makeblastdb -in $missing_pdb -dbtype prot -out $db_name
	blastp -subject $missing_pdb -query $targ_prot_fasta -outfmt "6 std qseq sseq" -out $output_file -max_hsps 1 -num_alignments 1000000 -seg no #get only the highest scoring pair for each hit and allow for as many alignments as possible because we want all the subjects to be present
	done