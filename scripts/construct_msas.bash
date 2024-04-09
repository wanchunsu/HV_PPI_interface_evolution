### Do Steps 1 and 2 manually    
### Step 1. Step 2. Make list of species (scientific name) and a mapping file (common name \t scientific name) 
### Step 2. Make species tree with the target protein and the orthologs of interest ### (Done manually before running this script)


### Step 3 (a - d)
output_folder_id=$1 # e.g. 19_orthologs 
species_tree_file=$2 # e.g. taxnames_19_orthologs.nwk 
ortholog_folder_label=$3 # e.g. orthologs 

### Step 3a. Run parse species orthologs into fasta file for MSA for "v_target_h" ###

# ----- Uncomment to run -----
python3 parse_human_selected_species_orthologs_into_fasta_file_for_MSA.py -f "v_target_h" -i $output_folder_id -o $ortholog_folder_label

# ----------------------------

### Step 3b. Run MAFFT to get MSA for the target sequence and its orthologs  ###
## Install MAFFT: See https://mafft.cbrc.jp/alignment/software/ for installation guidelines (I'm using https://mafft.cbrc.jp/alignment/software/ubuntu_on_windows.html)
# ----- Uncomment to run -----
type_analysis=v_target_h

echo \n Running analysis for $type_analysis	
orthologs_fasta_folder="../data/homologs_and_conservation/$ortholog_folder_label/$type_analysis/fasta_fis_for_MSA_$output_folder_id"
msa_results_folder="../data/homologs_and_conservation/$ortholog_folder_label/$type_analysis/MSA_$output_folder_id"
mkdir $msa_results_folder

for targ_fasta in $orthologs_fasta_folder/*; do
	targ_fasta_with_ext=${targ_fasta##*/}
	targ_fasta_id="${targ_fasta_with_ext%.*}"
	targ_fasta_id="${targ_fasta_id%%_*}"
	msa_out_fi=$msa_results_folder/$targ_fasta_id"_msa.fasta"

	# echo $targ_fasta
	# echo $msa_out_fi
	echo Running MSA for $targ_fasta_id
	mafft --auto $targ_fasta > $msa_out_fi
	
done
