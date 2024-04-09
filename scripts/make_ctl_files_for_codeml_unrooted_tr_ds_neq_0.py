import shutil
import argparse
import os.path as osp
import os
from pathlib import Path


'''
Things to modify:

seqfile = stewart.aa * sequence data filename

treefile = stewart.trees      * tree structure file name

outfile = mlc           * main result file name

seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs

model = 2
           * models for codons:
               * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
           * models for AAs or codon-translated AAs:
               * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
               * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

               (model should be 0 for site models)

'''

def make_copies_of_ctl_fi(ex_ctl_fi, residue_type, binned_dnds_folder, tree_fi):
    """ Make control files for codeml (we will read in the original example control file and edit lines as necessary then output as a new control file)
    We will make one ctl file for each residue type 
    Args:

        ex_ctl_fi: original example control file
        residue_type: residue type that we want to make a control file for ('exo', 'endo', 'mimicry', 'surface', 'buried', 'all_exo', 'all_endo')
        codons_folder: folder containing codon sequences to be fed into codeml(seqfile)
        tree_fi: tree file to feed into codeml (treefile)
        modified_ctl_fi_folder: new modified control file to output

    """
    seqfi = osp.join(residue_type + '.binned_codons_ds_neq_0.fasta')
    modified_ctl_fi = osp.join(binned_dnds_folder, residue_type + '.codeml_unrooted_ds_neq_0.ctl')
    #modified_ctl_fi = osp.join(binned_dnds_folder, residue_type + '.codeml.ctl')

    codeml_outfi = residue_type + '.codeml_unrooted_ds_neq_0.out' #osp.join(modified_ctl_fi_folder, residue_type + '.codeml.out')
    #codeml_outfi = residue_type + '.codeml.out' #osp.join(modified_ctl_fi_folder, residue_type + '.codeml.out')
    
    with open(ex_ctl_fi) as c, open(modified_ctl_fi, 'w') as mc:
        for line in c:
            if 'seqfile = stewart.aa' in line: #replace input seqfile with our own codon fasta file
                modified_seqfi_line = line.replace('stewart.aa', seqfi) 
                mc.write(modified_seqfi_line)
            elif 'treefile = stewart.trees' in line:#replace treefile with our own treefile
                modified_treefi_line = line.replace('stewart.trees', tree_fi) 
                mc.write(modified_treefi_line)
            elif 'outfile = mlc' in line: #replace mlc with codeml output fi
                modified_outfi_line = line.replace('mlc', codeml_outfi) 
                mc.write(modified_outfi_line)
            elif 'seqtype = 2' in line: #modify to 1 for codons
                modified_seqtype = line.replace('= 2', '= 1')
                mc.write(modified_seqtype)
            elif 'model = 2' in line: #modify to model=0 for site models
                modified_model = line.replace('= 2', '= 0')
                mc.write(modified_model)
            elif 'getSE = 0' in line: # modify to getSE=1 so that we have SE value
                        modified_SE = line.replace('= 0', '= 1')
                        mc.write(modified_SE)
            else: #if not a line that we need to modify
                mc.write(line) #write all other lines unchanged



def main():
    script_dir = osp.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--paml_folder')

    args = parser.parse_args()
    paml_folder_name = args.paml_folder

    binned_dnds_folder = osp.join('..', 'data', 'binned_dnds')

    ex_ctl_fi = osp.join(paml_folder_name, 'examples', 'codeml.ctl')

    residue_types = ['exo', 'endo', 'mimicry', 'surface', 'buried', 'all_exo', 'all_endo']

    tree_fi = osp.join('unrooted_binned_dnds.species_tree.nwk')
    #tree_fi = osp.join('binned_dnds.species_tree.nwk')



    print('Making ctl files for codeml . . .')
    for residue_type in residue_types:
        make_copies_of_ctl_fi(ex_ctl_fi, residue_type, binned_dnds_folder, tree_fi)



if __name__ == '__main__':
    main()


