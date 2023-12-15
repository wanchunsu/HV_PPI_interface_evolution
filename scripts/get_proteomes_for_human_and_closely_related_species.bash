#!/bin/bash
# Get proteomes for species closely related to humans 

orthologs_dir="../data/homologs_and_conservation/orthologs"
proteomes_dir="../data/homologs_and_conservation/orthologs/proteomes"
mkdir -p $orthologs_dir
mkdir -p $proteomes_dir #to store the proteomes for species other than human

#download proteomes (reviewed for human, both reviewed and unreviewed for other organisms)
echo Downloading proteomes . . .
wget -O $orthologs_dir/human_proteome_reviewed.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=(reviewed:true)%20AND%20(organism_id:9606)&format=fasta'
wget -O $proteomes_dir/mouse_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000000589&format=fasta'
wget -O $proteomes_dir/rat_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000002494&format=fasta'
wget -O $proteomes_dir/chimpanzee_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000002277&format=fasta'
wget -O $proteomes_dir/cow_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000009136&format=fasta'
wget -O $proteomes_dir/gorilla_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001519&format=fasta'
wget -O $proteomes_dir/orangutan_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001595&format=fasta'
wget -O $proteomes_dir/macaque_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000006718&format=fasta'
wget -O $proteomes_dir/dog_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000002254&format=fasta'
wget -O $proteomes_dir/rabbit_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001811&format=fasta'
wget -O $proteomes_dir/pig_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000008227&format=fasta'
wget -O $proteomes_dir/cat_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000011712&format=fasta'
wget -O $proteomes_dir/horse_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000002281&format=fasta'
wget -O $proteomes_dir/goat_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000291000&format=fasta'
wget -O $proteomes_dir/bonobo_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000240080&format=fasta'
wget -O $proteomes_dir/drill_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000233140&format=fasta'
wget -O $proteomes_dir/gibbon_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001073&format=fasta'
wget -O $proteomes_dir/squirrel-monkey_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000233220&format=fasta'
wget -O $proteomes_dir/baboon_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000028761&format=fasta'
wget -O $proteomes_dir/marmoset_proteome_all.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000008225&format=fasta'

