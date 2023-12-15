#!/bin/bash

# BLAST human proteome sequences against proteomes of other species
# and vice versa

orthologs_dir="../data/homologs_and_conservation/orthologs"
proteomes_dir="../data/homologs_and_conservation/orthologs/proteomes"
human_proteome=$orthologs_dir/human_proteome_reviewed.fasta
human_targets=$orthologs_dir/human_target_proteins.fasta

blast_dbs='../data/homologs_and_conservation/blast_dbs'
results_dir='../data/blast/results/reciprocal_best_hit'
mkdir -p $results_dir
mkdir -p $blast_dbs


#Make Blast database for human proteome:
makeblastdb -in $human_proteome -dbtype prot -out $blast_dbs/human_db


#Run through the proteomes for closely related species and do reciprocal blast 
# (1st blast: human targets --> species proteome db. 2nd blast: species proteome --> human_reviewed_db)
for proteome in $proteomes_dir/*_all.fasta
do 
        file_with_ext=${proteome##*/}
		
		organism_name=$(echo $file_with_ext | cut -d_ -f1)
		
		if [ ! -f $blast_dbs/$organism_name"_db.pdb" ]; then #if one of the blastdb files for this organism doesn't exist that means we haven't ran blast for this organism yet
			echo $organism_name

			# make blast database for the closely related species proteome
			makeblastdb -in $proteome -dbtype prot -out $blast_dbs/$organism_name"_db"

			# Blasting human targets against closely related species db (no need to use whole human proteome, because we're only interested in the human targets):
			echo Blasting human targets against $organism_name db
			blastp -db $blast_dbs/$organism_name"_db" -query $human_targets -outfmt "6 std qseq sseq" -out $results_dir/"human_against_"$organism_name".blast.out" -max_hsps 1 -max_target_seqs 5

			# Blasting species proteome against human db:
			echo Blasting $organism_name proteome against human db
			blastp -db $blast_dbs/human_db -query $proteome -outfmt "6 std qseq sseq" -out $results_dir/$organism_name"_against_human.blast.out" -max_hsps 1 -max_target_seqs 5
		fi

done
