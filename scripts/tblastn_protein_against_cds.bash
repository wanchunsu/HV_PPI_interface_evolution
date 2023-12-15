#!/bin/bash

ortholog_folder_label=$1
orthologs_dir=../data/homologs_and_conservation/$ortholog_folder_label
human_targ_prots_fi=../data/homologs_and_conservation/orthologs/human_target_proteins.fasta
species_ortholog_proteins_dir=$orthologs_dir/species_ortholog_prots_fasta
all_cds_dir=$orthologs_dir/cds

blast_dbs=$species_ortholog_proteins_dir/blast_dbs
results_dir='../data/blast/results/protein_to_cds'

mkdir -p $blast_dbs
mkdir -p $results_dir

# Perform tblastn for human target proteins against human cds sequences (to get the corresponding nucleotide protein coding sequence)
if [ -z "$(ls -A $all_cds_dir)" ]; then
   echo No new organisms to perform tblastn on!
else
	for cds_fi in $all_cds_dir/*.cds.fa; do
		cds_fi_with_ext=${cds_fi##*/}
		organism="${cds_fi_with_ext%%.*}"
		
		echo Current organism: $organism . . .
		
		echo  . . . Making blastdb for $organism

		makeblastdb -in $cds_fi -dbtype nucl -out $blast_dbs/$organism"_db"

		echo  . . . Performing tblastn to map proteins to their cds sequences

		if [ $organism = "human" ]
		then
			tblastn -db $blast_dbs/$organism"_db" -query $human_targ_prots_fi -outfmt "6 std qseq sseq stitle" -out $results_dir/$organism".prot_against_cds.blast.out" -max_hsps 1 -max_target_seqs 5 -seg no
			

		else
			tblastn -db $blast_dbs/$organism"_db" -query $species_ortholog_proteins_dir/$organism"_orthologs.fasta" -outfmt "6 std qseq sseq stitle" -out $results_dir/$organism".prot_against_cds.blast.out" -max_hsps 1 -max_target_seqs 5 -seg no
			
		fi
		# stitle in outfmt 6 gives the entire header of the subject (including white space)
		# -seg no is to prevent blast from filtering out amino acids that are low complexity (we might lose some residues if we filter, so set seg to no)
	done
fi


