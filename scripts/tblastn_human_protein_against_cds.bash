#!/bin/bash

targ_prot_fa_folder='../data/dnds/human_targ_prot_fasta'
targ_prot_cds_fa_folder='../data/dnds/human_targ_prot_cds_fasta_refseq'
results_dir='../data/blast/results/human_protein_to_cds'

mkdir -p $results_dir
for targ_fasta in $targ_prot_fa_folder/*; do
	targ_fasta_with_ext=${targ_fasta##*/}
	targ_fasta_id="${targ_fasta_with_ext%%.*}"
	cds_seq=$targ_prot_cds_fa_folder/$targ_fasta_id".cds.fasta"
	tblastn -query $targ_fasta -subject $cds_seq -outfmt "6 std qseq sseq" -out $results_dir/$targ_fasta_id".prot_against_cds.blast.out" -max_hsps 1 -seg no
done