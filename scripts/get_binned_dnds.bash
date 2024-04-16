#!/bin/bash

## This will generate two types of binned dn/ds results
# 1. all residues 
# 2. ds_neq_0: after removing special cases (residues with dS == 0, i.e., dn/ds is undefined) -- we will be taking these results


## 1. Download PAML ###
# ----- Uncomment to run -----	
# wget 'https://github.com/abacus-gene/paml/releases/download/v4.10.5/paml-4.10.5-linux-X86_64.tgz'; tar -xvzf paml-4.10.5-linux-X86_64.tgz
# ---------------------------- 


### 2. Bin codon alignments ###
# # ----- Uncomment to run -----
python3 bin_codon_alignments.py # all residues

python3 bin_codon_alignments_ds_neq_0.py # only residues where ds != 0 (based on ds values from site-specific method)
# # ---------------------------- 


### 3. Make tree file ###
# # ----- Uncomment to run -----
python3 make_and_label_species_tree_for_binned_dnds.py
# # ---------------------------- 


### 4. Unroot trees ###
# # ----- Uncomment to run -----
Rscript unroot_trees.R
# # ---------------------------- 


### 5. Make control files for PAML (edit example control file and output new control file) ###
# # ----- Uncomment to run -----

python3 make_ctl_files_for_codeml_unrooted_tr.py -p paml-4.10.5
python3 make_ctl_files_for_codeml_unrooted_tr_ds_neq_0.py -p paml-4.10.5 #this is done for only residues where (dS!=0, there are indeed subsitutions and more than one non-gap residue --i.e. we don't take special cases in site-specific dnds)
# # ---------------------------- 


### 6. Run codeml to get dN/dS for the 5 residue types ###
# ----- Uncomment to run -----	
codeml_ctl_fis_folder='../data/binned_dnds'
cd $codeml_ctl_fis_folder
for ctl_fi in ./*codeml_unrooted.ctl; do
	../../scripts/paml-4.10.5/bin/codeml $ctl_fi
done
cd ../../scripts


# Running for residues that aren't special cases (i.e. dS!=0)
codeml_ctl_fis_folder='../data/binned_dnds'
cd $codeml_ctl_fis_folder
for ctl_fi in ./*codeml_unrooted_ds_neq_0.ctl; do
	../../scripts/paml-4.10.5/bin/codeml $ctl_fi
done
cd ../../scripts

# ---------------------------- 


### 7. Extract dN/dS values from codeml output
# ----- Uncomment to run -----	
python3 extract_binned_dnds.py  -f codeml_unrooted
python3 extract_binned_dnds.py  -f codeml_unrooted_ds_neq_0
# ---------------------------- 