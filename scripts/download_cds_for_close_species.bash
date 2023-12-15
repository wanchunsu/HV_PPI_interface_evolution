#!/bin/bash

# Get coding sequences (CDS) for the closely related species (this is for blasting our orthologous proteins against so that we can find their corresponding nucleotide cds)
cds_dir="../data/homologs_and_conservation/orthologs/cds"

mkdir -p $cds_dir

echo Downloading cds fasta sequences . . . 

wget -O $cds_dir/human.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz'; gunzip $cds_dir/human.cds.fa.gz
wget -O $cds_dir/mouse.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz'; gunzip $cds_dir/mouse.cds.fa.gz
wget -O $cds_dir/squirrel-monkey.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/saimiri_boliviensis_boliviensis/cds/Saimiri_boliviensis_boliviensis.SaiBol1.0.cds.all.fa.gz'; gunzip $cds_dir/squirrel-monkey.cds.fa.gz
wget -O $cds_dir/marmoset.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/callithrix_jacchus/cds/Callithrix_jacchus.mCalJac1.pat.X.cds.all.fa.gz'; gunzip $cds_dir/marmoset.cds.fa.gz
wget -O $cds_dir/baboon.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/papio_anubis/cds/Papio_anubis.Panubis1.0.cds.all.fa.gz'; gunzip $cds_dir/baboon.cds.fa.gz
wget -O $cds_dir/gibbon.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/nomascus_leucogenys/cds/Nomascus_leucogenys.Nleu_3.0.cds.all.fa.gz'; gunzip $cds_dir/gibbon.cds.fa.gz
wget -O $cds_dir/drill.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/mandrillus_leucophaeus/cds/Mandrillus_leucophaeus.Mleu.le_1.0.cds.all.fa.gz'; gunzip $cds_dir/drill.cds.fa.gz
wget -O $cds_dir/bonobo.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/pan_paniscus/cds/Pan_paniscus.panpan1.1.cds.all.fa.gz'; gunzip $cds_dir/bonobo.cds.fa.gz
wget -O $cds_dir/goat.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/capra_hircus/cds/Capra_hircus.ARS1.cds.all.fa.gz'; gunzip $cds_dir/goat.cds.fa.gz
wget -O $cds_dir/horse.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/equus_caballus/cds/Equus_caballus.EquCab3.0.cds.all.fa.gz'; gunzip $cds_dir/horse.cds.fa.gz
wget -O $cds_dir/cat.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/felis_catus/cds/Felis_catus.Felis_catus_9.0.cds.all.fa.gz'; gunzip $cds_dir/cat.cds.fa.gz
wget -O $cds_dir/chimpanzee.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/pan_troglodytes/cds/Pan_troglodytes.Pan_tro_3.0.cds.all.fa.gz'; gunzip $cds_dir/chimpanzee.cds.fa.gz
wget -O $cds_dir/cow.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/bos_taurus/cds/Bos_taurus.ARS-UCD1.2.cds.all.fa.gz'; gunzip $cds_dir/cow.cds.fa.gz
wget -O $cds_dir/dog.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/canis_lupus_familiaris/cds/Canis_lupus_familiaris.ROS_Cfam_1.0.cds.all.fa.gz'; gunzip $cds_dir/dog.cds.fa.gz
wget -O $cds_dir/gorilla.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/gorilla_gorilla/cds/Gorilla_gorilla.gorGor4.cds.all.fa.gz'; gunzip $cds_dir/gorilla.cds.fa.gz
wget -O $cds_dir/macaque.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/macaca_mulatta/cds/Macaca_mulatta.Mmul_10.cds.all.fa.gz'; gunzip $cds_dir/macaque.cds.fa.gz
wget -O $cds_dir/orangutan.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/pongo_abelii/cds/Pongo_abelii.Susie_PABv2.cds.all.fa.gz'; gunzip $cds_dir/orangutan.cds.fa.gz
wget -O $cds_dir/pig.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/sus_scrofa/cds/Sus_scrofa.Sscrofa11.1.cds.all.fa.gz'; gunzip $cds_dir/pig.cds.fa.gz
wget -O $cds_dir/rabbit.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/oryctolagus_cuniculus/cds/Oryctolagus_cuniculus.OryCun2.0.cds.all.fa.gz'; gunzip $cds_dir/rabbit.cds.fa.gz
wget -O $cds_dir/rat.cds.fa.gz 'http://ftp.ensembl.org/pub/release-106/fasta/rattus_norvegicus/cds/Rattus_norvegicus.mRatBN7.2.cds.all.fa.gz'; gunzip $cds_dir/rat.cds.fa.gz
