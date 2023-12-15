#!/bin/bash

# Code for dN/dS pipeline (includes MSAs -- aa and codon, building trees, calculating and parsing dN/dS)
### we will be using just the species tree analysis (RAXML tree analysis is commented out)

### 1. Run MSAs for human and ortholog amino acid sequences ###
# ----- Uncomment to run -----	 
orthologs_fasta_folder="../data/dnds/aa_fasta"
msa_results_folder="../data/dnds/aa_msas"
mkdir $msa_results_folder

for targ_fasta in $orthologs_fasta_folder/*; do
	targ_fasta_with_ext=${targ_fasta##*/}
	targ_fasta_id="${targ_fasta_with_ext%%.*}"
	
	msa_out_fi=$msa_results_folder/$targ_fasta_id"_msa.fasta"

	# echo $targ_fasta
	# echo $msa_out_fi
	echo Running MSA for $targ_fasta_id
	mafft --auto $targ_fasta > $msa_out_fi
	
done
# ---------------------------- 


### 2. Download pal2nal (Suyama et al., 2006) for translating aa alignments into dna alignments ###
# ----- Uncomment to run -----
wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz; tar xzf pal2nal.v14.tar.gz; rm pal2nal.v14.tar.gz
# ---------------------------- 


### 3. Run pal2nal to translate the amino acid alignments into codon alignments ###
# ----- Uncomment to run -----
aa_msas_folder='../data/dnds/aa_msas'
dna_seqs_folder='../data/dnds/dna_fasta'
dna_msas_outfolder='../data/dnds/dna_msas'

mkdir -p $dna_msas_outfolder

for targ_fasta_msa in $aa_msas_folder/*; do
	targ_fasta_with_ext=${targ_fasta_msa##*/}
	targ_fasta_id="${targ_fasta_with_ext%.*}"
	targ_fasta_id="${targ_fasta_id%%_*}"

	dna_fasta=$dna_seqs_folder/$targ_fasta_id'.dna.fasta'
	outfi=$dna_msas_outfolder/$targ_fasta_id'_dna_msa.fasta'

	echo Translating aa alignments to codon for: $targ_fasta_msa
	
	# Translate amino acid alignments into codon alignments
	./pal2nal.v14/pal2nal.pl $targ_fasta_msa $dna_fasta -output fasta > $outfi 
	
done
# ---------------------------- 


### 4. clone repo for RAxML (Stamatakis, 2014) ###
# ----- Uncomment to run -----
git clone 'https://github.com/stamatak/standard-RAxML.git'
cd standard-RAxML; make -f Makefile.SSE3.gcc; rm *.o; cd ..
# ---------------------------- 


### 5. Two methods for making trees (a)Raxml (b)species tree

# ### 5a. Run RAxML on aa alignments ###
# # ----- Uncomment to run -----
# aa_msas_folder='../data/dnds/aa_msas'
# #out_tree_dir='../data/dnds/raxml_out'
# out_tree_dir="/mnt/c/Users/user/OneDrive - McGill University/Desktop/QLS/Xia_Lab/data/dnds/raxml_out" # output dir must be absolute path (RAxML doesn't work with relative paths)

# mkdir -p "$out_tree_dir"

# for targ_fasta_msa in $aa_msas_folder/*_msa.fasta; do
# 	targ_fasta_with_ext=${targ_fasta_msa##*/}
# 	targ_fasta_id="${targ_fasta_with_ext%.*}"
# 	targ_fasta_id="${targ_fasta_id%%_*}"

# 	tree_ext=$targ_fasta_id'_tree' #extension for tree output (note: output will be RAxML_bestTree.<tree_ext>)
	
# 	echo Running RAxML on aa alignment for $targ_fasta_id

# 	standard-RAxML/raxmlHPC-SSE3 -s $targ_fasta_msa -n $tree_ext -m PROTCATLG -p 12345 -w "$out_tree_dir"
	
# 	#check that the best tree was made, if not, something went wrong when making the tree, so redo
# 	if [ ! -f "$out_tree_dir"/'RAxML_bestTree.'$tree_ext ]; then 

# 		echo ----------- \n\n Could not find best tree for $targ_fasta_id, so re-running RAxML -----------\n\n

# 		rm "$out_tree_dir"/*.$tree_ext #remove prev files for this target id

# 		#Re-run RAxML
# 		standard-RAxML/raxmlHPC-SSE3 -s $targ_fasta_msa -n $tree_ext -m PROTCATLG -p 12345 -w "$out_tree_dir"
# 	fi

# done
# # ----------------------------

### 5b. Make species trees for aa alignments
aa_msas_folder='../data/dnds/aa_msas'
for targ_fasta in $aa_msas_folder/*; do
		targ_fasta_with_ext=${targ_fasta##*/}
		targ_fasta_id="${targ_fasta_with_ext%.*}"
		
		targ_fasta_id="${targ_fasta_id%%_*}"
		# -m: msa file (w/out pathname)
		# -s: species tree file (w/out pathname)
		# -o: output folder for relabeled species tree (w/out pathname)

		echo Making and relabeling species trees
		python3 make_and_label_species_tree_for_dnds.py -m $targ_fasta_id'_msa.fasta' -s taxnames_19_orthologs.nwk -o species_tree_out

done


### 6. Install HyPhy as detailed in https://github.com/veg/hyphy  ###
# ----- Uncomment to run -----
git clone https://github.com/veg/hyphy.git #I'm using hyphy 2.5.40
cd hyphy
cmake .
sudo make install
cd ..
# ----------------------------


### 7. Run hyphy FEL 

## Using Raxml tree
# ----- Uncomment to run -----
# dna_msas_folder='/mnt/c/Users/user/OneDrive - McGill University/Desktop/QLS/Xia_Lab/data/dnds/dna_msas'
# tree_folder='/mnt/c/Users/user/OneDrive - McGill University/Desktop/QLS/Xia_Lab/data/dnds/raxml_out'
# hyphy_output='/mnt/c/Users/user/OneDrive - McGill University/Desktop/QLS/Xia_Lab/data/dnds/hyphy_out'
# mkdir -p "$hyphy_output"

# for targ_dna_msa in "$dna_msas_folder"/*_dna_msa.fasta; do
# 	targ_fasta_with_ext=${targ_dna_msa##*/}
# 	targ_fasta_id="${targ_fasta_with_ext%.*}"
# 	targ_fasta_id="${targ_fasta_id%%_*}"


# 	tree_fi="$tree_folder"/'RAxML_bestTree.'$targ_fasta_id'_tree'
# 	output_fi="$hyphy_output"/$targ_fasta_id"_dna_msa.fasta.FEL.json"
	

# 	# Run Hyphy with command line inputs to infer rates 
	
# 	echo Running hyphy for $targ_fasta_id . . .
# 	hyphy fel --code Universal \
# 	 --alignment "$targ_dna_msa" \
# 	 --tree "$tree_fi" \
# 	 --branches All \
# 	 --srv No \
# 	 --pvalue 0.1 \
# 	 --output "$output_fi"
	
# 	# If hyphy output happens to be empty (b/c of any errors), re-run
# 	if [ ! -f "$output_fi" ]; then
# 		echo ----------- \n\n Re-running hyphy for $targ_fasta_id -----------\n\n
# 		rm "$output_fi"
# 		hyphy fel --code Universal \
# 		 --alignment "$targ_dna_msa" \
# 		 --tree "$tree_fi" \
# 		 --branches All \
# 		 --srv No \
# 		 --pvalue 0.1 \
# 		 --output "$output_fi"
# 	fi
# done
# ----------------------------

## Using species tree
# ----- Uncomment to run -----
dna_msas_folder='/mnt/c/Users/user/OneDrive - McGill University/Desktop/QLS/Xia_Lab/data/dnds/dna_msas'
tree_folder='/mnt/c/Users/user/OneDrive - McGill University/Desktop/QLS/Xia_Lab/data/dnds/species_tree_out'
hyphy_output='/mnt/c/Users/user/OneDrive - McGill University/Desktop/QLS/Xia_Lab/data/dnds/hyphy_out_species_tree'
mkdir -p "$hyphy_output"

for targ_dna_msa in "$dna_msas_folder"/*_dna_msa.fasta; do
	targ_fasta_with_ext=${targ_dna_msa##*/}
	targ_fasta_id="${targ_fasta_with_ext%.*}"
	targ_fasta_id="${targ_fasta_id%%_*}"


	tree_fi="$tree_folder"/$targ_fasta_id'.species_tree.nwk'
	output_fi="$hyphy_output"/$targ_fasta_id"_dna_msa.fasta.FEL.json"
	

	# Run Hyphy with command line inputs to infer rates 
	
	echo Running hyphy for $targ_fasta_id . . .
	hyphy fel --code Universal \
	 --alignment "$targ_dna_msa" \
	 --tree "$tree_fi" \
	 --branches All \
	 --srv No \
	 --pvalue 0.1 \
	 --output "$output_fi"
	
	# If hyphy output happens to be empty (b/c of any errors), re-run
	if [ ! -f "$output_fi" ]; then
		echo ----------- \n\n Re-running hyphy for $targ_fasta_id -----------\n\n
		rm "$output_fi"
		hyphy fel --code Universal \
		 --alignment "$targ_dna_msa" \
		 --tree "$tree_fi" \
		 --branches All \
		 --srv No \
		 --pvalue 0.1 \
		 --output "$output_fi"
	fi
done
# ----------------------------


### 8. Clone repo containing scripts for site-specific dn/ds calculations (Sydykova et al., 2018) ###
# ----- Uncomment to run -----
git clone 'https://github.com/clauswilke/proteinER.git'
# ----------------------------


## Note the following step makes use of custom script files from the proteinER repo (Sydykova et al., 2018) 

### 9. Parse Hyphy output  ###
# Using raxml hyphy output
# ----- Uncomment to run -----
# hyphy_output='../data/dnds/hyphy_out'
# aa_msas_folder='../data/dnds/aa_msas'
# for targ_fasta_hyphy in "$hyphy_output"/*_dna_msa.fasta.FEL.json; do
# 	targ_fasta_with_ext=${targ_fasta_hyphy##*/}
# 	targ_fasta_id="${targ_fasta_with_ext%%.*}"
# 	targ_fasta_id="${targ_fasta_id%%_*}"
	
	
# 	# a) Parse Hyphy output and extract dN and dS values
# 	parsed_data="$hyphy_output"/$targ_fasta_id'.parsed_dnds.csv'

# 	echo Parsing hyphy output for $targ_fasta_id
# 	python3 ./proteinER/src/parse_FEL.py -j "$targ_fasta_hyphy" -r "$parsed_data"

# 	#This is the older version (now, we process dN and dS and calculate dN/dS in the next step on our own -- b/c we want to treat dS=0 as special cases)
# 	# # b) Now take parsed hyphy data and calculate dN/dS values for each residue
# 	# processed_data="$hyphy_output"/$targ_fasta_id'.processed_dnds.csv'
# 	# aa_msa=$aa_msas_folder/$targ_fasta_id'_msa.fasta'

# 	# echo Processing and calculating site-wise dN/dS values for $targ_fasta_id
# 	# python3 ./proteinER/src/calc_dNdS.py -a $aa_msa -r "$parsed_data" -o "$processed_data"
# done
# ----------------------------

# Using species tree hyphy output
# ----- Uncomment to run -----
hyphy_output='../data/dnds/hyphy_out_species_tree'
aa_msas_folder='../data/dnds/aa_msas'
for targ_fasta_hyphy in "$hyphy_output"/*_dna_msa.fasta.FEL.json; do
	targ_fasta_with_ext=${targ_fasta_hyphy##*/}
	targ_fasta_id="${targ_fasta_with_ext%%.*}"
	targ_fasta_id="${targ_fasta_id%%_*}"
	
	
	# a) Parse Hyphy output and extract dN and dS values
	parsed_data="$hyphy_output"/$targ_fasta_id'.parsed_dnds.csv'

	echo Parsing hyphy output for $targ_fasta_id
	python3 ./proteinER/src/parse_FEL.py -j "$targ_fasta_hyphy" -r "$parsed_data"

	
done
# ----------------------------


### 10. Process and calculate dN/dS values and map back onto our actual target proteins of interest ###
# ----- Uncomment to run -----

# python3 process_dn_ds.py -f hyphy_out -o dnds_on_targ_prots #no srv, using raxml hyphy results

python3 process_dn_ds.py -f hyphy_out_species_tree -o dnds_on_targ_prots_species_tree #no srv, using species tree hyphy results
# ----------------------------


### 11. Categorize dN/dS values for our 5 different residues
# ----- Uncomment to run -----
echo Running categorize_dnds_and_special_cases_into_residue_groups.py

#python3 categorize_dnds_and_special_cases_into_residue_groups.py -f hyphy_out #no srv, using raxml hyphy results

python3 categorize_dnds_and_special_cases_into_residue_groups.py -f hyphy_out_species_tree #no srv, using species tree hyphy results

# ----------------------------