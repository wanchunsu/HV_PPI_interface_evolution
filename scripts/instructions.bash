
## Instructions for building a human-virus SIN, mapping and categorizing interfaces, and performing evolutionary analyses

# Please read the instructions carefully. Some steps may require manual work.
# Data will be downloaded into the `../data/` directory. Results will also be stored in this directory.

#----------------- Building a human-virus SIN and determining interfaces -----------------#

### PART 1: Exogenous PPIs (human-virus PPIs) ###

### Downloading and processing IntAct data for human-virus binary PPIs

# 1. Go to http://www.ebi.ac.uk/legacy-intact/query/annot:%22dataset:virus%22?conversationContext=8. This is the virus PPI dataset provided by IntAct.

# 2. Make folder to save data: 
mkdir ../data/IntAct_human_virus

# 3. Download and save data in `../data/IntAct_human_virus` folder. 
# 'Select format to Download' --> MI-TAB 2.5 --> 'Download' --> save as `IntAct_virus_ppis.txt`

# 4. Get list of taxids corresponding to viruses that infect humans. 
#	i) Go to https://www.ncbi.nlm.nih.gov/taxonomy/?term=txid10239[Subtree]. This is a listing of subtree for 10239 taxID (virus superkingdom taxid).
#	ii) 'Send to' drop down --> 'File' --> Under 'Format' --> Choose 'Taxid list'
# 	iii) 'Create File' --> Save the file in the directory `../data/IntAct_human_virus` as `all_virus_taxonomy.txt`

# 5. Format and process IntAct file above. This produces lists of human and virus uniprot ids (`../data/IntAct_human_virus/unique_human_list.txt` and `../data/IntAct_human_virus/unique_virus_list.txt`)
python3 format_exogenous_intact_data.py

# 6. Feed `unique_human_list.txt` and `unique_virus_list.txt` separately into UniProt's ID mapping tool (https://www.uniprot.org/id-mapping)
# Download the information for these proteins in tsv format (with only the "Reviewed" status and "Sequence" columns selected). 
# Save the files as `../data/IntAct_human_virus/human_uniprot.tsv` and `../data/IntAct_human_virus/virus_uniprot.tsv`.

# 7. Get valid exogenous ppis and their corresponding sequences
python3 get_valid_exo_ppis_and_output_fa.py


### Sequence alignment for homology modeling (between human proteins and PDB sequences, and between virus proteins and PDB sequences)

# 1. Obtain all PDB sequence information
mkdir ../data/PDB
wget -O ../data/PDB/pdb_seqres.txt.gz 'https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz'
gunzip ../data/PDB/pdb_seqres.txt.gz

# 2. Run BLAST (if using BLAST for first time, see https://www.ncbi.nlm.nih.gov/books/NBK52640/ on how to install and set paths)
bash blast_human_and_virus_against_pdb.bash

# 3. Filter E-values (using 1e-10 for human protein and 1e-5 for virus proteins)
bash filter_exo_blast_evalues.bash


### Getting PDB structural information and finding interfacial residues

# 1. Download pdb resolution file: 
wget -O ../data/PDB/resolu.idx 'https://files.wwpdb.org/pub/pdb/derived_data/index/resolu.idx'

# 2. Format blast results and map pdb structures to human virus PPIs, and download corresponding PDB ids:  
python3 get_pdb_structures_for_hv_interacting_pairs.py

# 3. Get interfacial residues using SASA
#	i) Download DSSP along with dependencies (see https://github.com/PDB-REDO/dssp for more info): 
	bash install_dssp.bash

#	ii) Run DSSP to get the RSA values for the residues of a protein as a monomer and complex (in complex with partner chain) and calculate delta RSA (RSA for monomer – RSA for complex)
	python3 get_interfacial_residues_using_sasa.py \
   	-i ../data/exo_structural_annotation/human_virus_ppis_mapped_to_pdb_chains.tsv \
   	-o ../data/exo_structural_annotation/interface_residues/exo_hv_interface_residues.tsv \
   	-t exo


### Get PDB to uniprot mapping files (this will be used for finding additional H-V PPI data that may be missing in IntAct)

# 1. Get fasta file for reviewed human proteins.
mkdir -p ../data/pdb_annotations/fasta_files
wget -O ../data/pdb_annotations/fasta_files/human_reviewed_uniprot.fasta 'https://rest.uniprot.org/uniprotkb/stream?query=(reviewed:true)%20AND%20(organism_id:9606)&format=fasta'

# 2. Get list of reviewed human-affecting virus proteins.
#	i) Go to https://www.uniprot.org/uniprotkb?facets=reviewed%3Atrue&query=%28virus_host_id%3A9606%29  
#	ii) Click 'Download' (make sure 'Compressed' is set to 'No') and save into `../data/pdb_annotations/fasta_files` folder as `virus_affecting_humans_reviewed_uniprot.fasta`.

# 3. Format the fasta headers.
python3 format_fasta_headers_pdb_anno.py

# 4. Download and process taxonomy ID to PDB mappings.
wget -O ../data/pdb_annotations/pdb_chain_taxonomy.tsv.gz 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz'
gunzip ../data/pdb_annotations/pdb_chain_taxonomy.tsv.gz
bash process_pdb_chain_taxonomy_fi.bash

# 5. Download list of uniprot ID to PDB chain mappings.
wget -O ../data/pdb_annotations/pdb_chain_uniprot.tsv.gz 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz'
gunzip ../data/pdb_annotations/pdb_chain_uniprot.tsv.gz

# 6. Get list of taxids corresponding to human viruses.  
#	i) Go to: https://www.uniprot.org/taxonomy?query=(host:9606)
#	ii) Click 'Download' (make sure 'Compressed' is set to 'No') and save into `../data/pdb_annotations` as `taxonomy_ids_for_viruses_infecting_humans.tab`

# 7. Get list of all reviewed uniprot ids 
#	i) Go to https://www.uniprot.org/uniprotkb?query=reviewed:true  
#	ii) Click 'Download' (make sure 'Compressed' is set to 'No', and select 'List' Format) and save into `../data/pdb_annotations` as `all_uniprot_reviewed_ids.txt` 


### Get additional PDB structures that have human-virus annotations 

# 1. Retrieve human-derived and virus-derived chains
python3 get_pdb_structures_for_annnotations.py

# 2. Blast human PDB chains against all reviewed human proteins
bash blast_human_pdb_against_human_uniprot.bash

# 3. Blast virus pdb chains against all reviewed human-affecting virus proteins.  
bash blast_virus_pdb_against_virus_uniprot.bash

# 4. Filter BLAST results by E-value (10^-10 for human, 10^-5 for virus).
bash filter_pdb_to_uniprot_evals.bash

# 5. Process human and virus PDB-to-UniProt BLAST results
python3 make_uniprot_to_pdb_mapping_files.py

# 6. Obtain PDB chains with human and virus annotations using BLAST results and PDB-to-UniProt annotations.
python3 obtain_pdb_with_human_virus_annotations.py


### Finding interfacial residues for PDB structures with human-virus annotations and create combined data with homology data

# 1. Get interfacial residues using SASA:  
python3 get_interfacial_residues_using_sasa.py \
-i ../data/pdb_annotations/pdbs_with_human_virus_annotations.txt \
-o ../data/pdb_annotations/interface_residues/pdb_hv_annotations_interface_residues.tsv \
-t exo

# 2. Filter off ppis-pdb-pdbchain combinations that are also seen in the homology modeling results, then combine results from homology modeling and PDB annotations.
python3 filter_and_combine_homology_and_pdb_annotation_results.py \
-m ../data/exo_structural_annotation/interface_residues/exo_hv_interface_residues.tsv \
-p ../data/pdb_annotations/interface_residues/pdb_hv_annotations_interface_residues.tsv \
-f ../data/pdb_annotations/interface_residues/filtered_pdb_hv_annotations_interface_residues.tsv \
-c ../data/combined_homology_and_pdb_anno_human_virus_interface_residues.tsv \
-t exo


### PART 2: Endogenous PPIs (human-human ppis, using human target proteins found above) ###

### Downloading and processing IntAct data for endogenous (human-human) binary PPIs

# 1. Download all data from intact
mkdir ../data/IntAct_human_human
wget -O ../data/IntAct_human_human/intact.zip 'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip'
unzip ../data/IntAct_human_human/intact.zip -d ../data/IntAct_human_human
rm ../data/IntAct_human_human/intact.zip
rm ../data/IntAct_human_human/intact_negative.txt

# 2. Format and process the intact data above
python3 format_endogenous_intact_data.py

# 3. Feed `unique_human_list.txt` into Uniprot's id mapping tool (https://www.uniprot.org/id-mapping) 
# Download the information for these proteins in tsv format (with only the Reviewed" status and "Sequence" columns selected).  
# Save the file as `../data/IntAct_human_human/human_uniprot.tsv`.

# 4. Get valid endgenous ppis and their corresponding sequences
python3 get_valid_endo_ppis_and_output_fa.py


### Sequence alignment for homology modeling (between human proteins and PDB sequences)

# 1. Run BLAST
bash blast_human_endo_prots_against_pdb.bash

# 2. Filter E-values (using 1e-10 for human proteins)
bash filter_endo_blast_evalues.bash


### Getting PDB structural information and finding interfacial residues

# Note: Our bottleneck is exogenous PPIs, so we don't need to use PDB annotations to find additional endogenous PPIs here

# 1. Format blast results and map pdb structures to human-human PPIs, and download corresponding PDB ids:  
python3 get_pdb_structures_for_endo_interacting_pairs.py

# 2.  Get interfacial residues using SASA:
python3 get_interfacial_residues_using_sasa.py \
-i ../data/endo_structural_annotation/human_endogenous_ppis_mapped_to_pdb_chains.tsv \
-o ../data/endo_structural_annotation/interface_residues/endo_human_interface_residues.tsv \
-t endo


### PART 3: Annotate antibody/mhc complexes ###

### Annotate exogenous ppis

# The goal here is to annotate all exogenous PPIs as either immunological (human targeting virus; 'T') or pathogenic (virus targeting human; 'F')
# and to remove immunological exogenous PPIs since not related to viral targeting of host proteins

# Note: Most of the steps here are require manual work. 

# 1. Obtain list of human proteins with exogenous ppis
python3 get_lists_of_exo_uniprot_pdbs_for_manual_antibody_mhc_anno.py
# This produces a list of the human protein uniprot ids (one per line) 
# followed by a tab character (stored as `..data/antibody_mhc_complex_annotations/exogenous_human_targs_with_antibody_mhc_annotations.tsv`).

# 2.  Manual annotation of target proteins. 
# For each protein in the produced file (`../data/antibody_mhc_complex_annotations/exogenous_human_targs_with_antibody_mhc_annotations.tsv`), 
# search each protein up on UniProt and annotate with a "T" (antibody/mhc protein) or "F" (non-antibody/mhc protein) after the tab character.
# (an example file is provided)

# 3. Automatically annotate pdbs where all mapped human targets are non-antibody/mhc.
python3 annotate_human_targs_and_pdbs.py

# 4. Manual annotation of remaining pdbs (these pdbs either contain all antibody/mhc proteins, or some antibody/mhc and some non-antibody/mhc proteins). 
# For each pdb in the produced file (`../data/antibody_mhc_complex_annotations/exogenous_ppi_pdbs_with_antibody_mhc_annotations.tsv`), 
# annotate as follows (by looking at PDB entry and going through `../data/combined_homology_and_pdb_anno_human_virus_interface_residues.tsv` to see what proteins a pdb chain maps to)
# 	a) Annotate a PDB as "T": if all mapped human targets are T and are targeting virus protein(s) (immune response).
#	b) Annotate a PDB as "F": if virus protein is targeting all mapped human targets (even if some or all targets are antibody/mhc)
#	c) Annotate a PDB as "D": if PDB has both "T" and "F" proteins mapped to it. Virus protein targets F proteins and T proteins (antibody/mhc proteins) target virus protein. 

# 5. Annotate exogenous PPIs.
python3 annotate_exo_ppis_with_antibody_mhc_anno.py


### PART 4: Verification of potential exogenous PPIs (from pdb annotations) ###

# 1. Output lists of pathogenic and immune exogenous PPIs to verify: 
python3 verify_potential_exo_ppis_patho_and_immune.py
# Note: the lists are stored in `data/pdb_annotations/verify_ppis` directory.
# We are only interested in pathogenic PPIs, so will only need to verify PPIs in the pathogenic file

# 2. Go through PPIs in `data/pdb_annotations/verify_ppis/pathogenic_exogenous_ppis_to_verify.tsv` 
# and manually search their associated PDB structures to confirm that this PPI is a verified pathogenic interaction. 
# 	a) Anntoate with "T" if verified 
#	b) Annotate with "F" otherwise

# 3. Make a new file containing only the known and verified pathogenic exogenous PPIs. 
# Also makes a new endogenous file with only target proteins that exist in the new exogenous file
python3 extract_verified_and_pathogenic_ppis.py


#----------------- Mapping interfacial residues back onto human target proteins and categorizing them -----------------#

### PART 5: Categorize endogenous, exogenous specific, and mimicked interfacial residues ###

### Get corresponding interfacial residues on the uniprot protein sequence

# 1. Compile all endogenous and exogenous results and output relevant fasta files for mapping between uniprot sequence and corresponding PDB chains
python3 merge_int_res_on_same_chain.py

# 2. Get interfacial residues based on 'label' residues. 
# Note: Current interfacial residues are 'auth' residues and pdb seqres (and blast results) are in 'label' residues.
python3 get_interfacial_residues_based_on_pdb_label_res.py

# 3. Run BLAST to map uniprot protein seq residues to corresponding pdb chain seq residues:  
bash blast_match_uniprot_residues_to_pdb_chain.bash

# 4. Format the above BLAST results into a dict and get missing pdbs that were not aligned (because they were too short).
python3 format_int_res_pdb_to_uniprot_blast_results.py

# 5. Run BLAST for missing pdb_chains (try to recover if possible):  
bash blast_missing_pdbs_uniprot_residues_to_pdb_chain.bash

# 6. Get corresponding interfacial residues on uniprot protein.  
python3 get_corresponding_int_res_on_uniprot.py


### Categorize interfacial residues

# Categorize interfacial residues into exogenous-specific and endogenous-specific and mimicry groups:
python3 categorize_interfacial_residues.py \
-x "exogenous_ppis_with_int_res_on_uniprot_seq.tsv" \
-n "endogenous_ppis_with_int_res_on_uniprot_seq.tsv" \
-f "v_target_h" 


#----------------- Evolutionary conservation of different PPI interfaces -----------------#

### PART 6: Sequence mismatches ###

### Get orthologs for each of the target proteins of interest and find fraction 
# of mismatches between interfacial residues on the target protein and its orthologs.

# 1. Get fasta sequences for human target proteins.
python3 make_fasta_with_human_target_prot_seqs.py

# 2. Download proteomes for humans (reviewed) and for other species (both reviewed and unreviewed b/c other proteomes might not be as extensive as human):
bash get_proteomes_for_human_and_closely_related_species.bash

# 3. Blast human and other species proteomes against each other:
bash blast_human_with_closely_related_species_reciprocal.bash

# 4. Get reciprocal best hit orthologs in the other organisms for each target protein and 
# calculate % mismatch for interfacial residue categories.
python3 get_divergence_of_orthologs_in_closely_related_species.py -f "v_target_h"

# 5. Compute bootstrapped standard error for mismatch values (for plotting purposes):
python3 bootstrap_for_avg_divergence_standard_err.py -f "v_target_h" -o orthologs

# 6. Calculate % mismatches on a per-protein basis.
python3 get_per_protein_divergence_of_orthologs_in_closely_related_species.py -f "v_target_h"

# 7. Compute per-protein bootstrapped standard error for mismatch values (for plotting purposes):
python3 bootstrap_for_avg_per_protein_divergence_standard_err.py -f "v_target_h" -o orthologs


### PART 7: Site-specific dN/dS ###

### Get phylogenetic information and make MSAs using species specific orthologs from above 
# (these results will be used for both site-specific and binned dN/dS)

# 1. Make list of scientific names for the selected organisms
# and create a mapping file with common name and corresponding scientific name (tab separated).
# These files are already included in ../data/homologs_and_conservation/orthologs folder as
# `taxnames_19_orthologs.txt` and `taxnames_mappings.tsv`

# 2. Make species tree using Timetree (http://www.timetree.org/).
# Feed in list of scientific names for the selected organisms (`taxnames_19_orthologs.txt`)
# The species tree is already included in ../data/homologs_and_conservation/orthologs folder as
# `taxnames_19_orthologs.nwk` 

# 3. Construct MSAs, run Rate4Site and calculate overall Rate4Site scores.
bash construct_msas.bash


### Get site-specific dN/dS 

# 1. Download the refseq-to-uniprot mapping file: 
mkdir ../data/dnds
wget -O ../data/dnds/gene_refseq_uniprotkb_collab.gz 'https://ftp.ncbi.nlm.nih.gov/refseq/uniprotkb/gene_refseq_uniprotkb_collab.gz'
gunzip ../data/dnds/gene_refseq_uniprotkb_collab.gz

# 2. Map between UniProt and RefSeq protein IDs for human target proteins
python3 map_uniprotid_to_refseq.py

# 3. Run tBLASTn using human target proteins as query and their mapped RefSeq cds as subject to get correct positions on cds and to check whether 
# the cds matches the protein exactly, if not will need to take the translated amino acid when we parse tblastn results.  
bash tblastn_human_protein_against_cds.bash

# 4. Download CDS DNA sequences for human and orthologous species from Ensembl.
bash download_cds_for_close_species.bash

# 5. Make fasta files with orthologs (one fasta file per species): 
python3 make_fasta_for_species_orthologs.py

# 6. Blast human target proteins and ortholog proteins against corresponding species cds sequences (from Ensembl)
# to get corresponding nucleotide protein coding sequences:
bash tblastn_protein_against_cds.bash orthologs

# 7. Parse tblastn results and get cds sequences for each protein. Place results into separate files for each target protein: 
python3 parse_tblastn_results.py

# 8. Create protein and codon MSAs, compute dN/dS, and calculate overall dN/dS
bash calc_site_specific_dnds.bash


### PART 8: Binned dN/dS ###

# 1. Make folder to store binned dnds input and output files
mkdir ../data/binned_dnds

# 2. Run binned dN/dS scripts
bash get_binned_dnds.bash