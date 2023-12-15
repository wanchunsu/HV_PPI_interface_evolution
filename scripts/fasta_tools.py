'''

Tools for formatting fasta files 
=================================================
Author: Wan-Chun Su (wan.su@mail.mcgill.ca)

'''

from Bio import SeqIO


def format_uniprot_fasta_headers(input_fasta, output_formatted_fasta, delimiter_to_split_by, first_ind_to_keep, last_ind_to_keep): 
    '''
    Formats the headers of a fasta file by first splitting it by the specified delimiter and then taking the elements starting from the first_ind_to_keep till the last_ind_to_keep (note these indices start from 0),
    
    Args:
        input_fasta (Path): path to fasta file to be formatted
        output_formatted_fast (Path): path to output of formatting
        delimiter_to_split_by (char): delimiter to split fasta file header by
        first_ind_to_keep (int): start index of split header to keep
        last_ind_to_keep (int): end index of split header to keep


    '''

    s = list(SeqIO.parse(str(input_fasta), 'fasta'))
    with open(output_formatted_fasta, "w") as fout:
        for _, row in enumerate(s):
                    idsplit = list(map(str.strip, row.id.split(delimiter_to_split_by)))
                    shortened_id = '|'.join(idsplit[ first_ind_to_keep  : last_ind_to_keep +1 ])
                    fout.write('>' + shortened_id + '\n')
                    fout.write(str(row.seq) + '\n')

def output_fasta_fi_from_list_of_fasta(fasta_list, outfi):
    ''' Outputs a fasta file from a list of SeqRecord objects

    Args:
        fasta_list (list): list of SeqRecord objects to output to file
        outfi: fasta file output path
    '''
    with open(outfi, "w") as of:
        SeqIO.write(fasta_list, of, "fasta")

def parse_fasta_file_into_dict(input_fasta, change_header=False, split_header_by='|', ind_to_keep = 1):
    
    dict_of_fasta = {}
    fastas = SeqIO.parse(open(input_fasta),'fasta') 
    for fasta in fastas:
        name, sequence = fasta.description, str(fasta.seq)
        if change_header == True:
            name = name.split(split_header_by)[ind_to_keep] #e.g. We want "Q86XN6" from: >sp|Q86XN6|ZN761_HUMAN Zinc finger protein ..." 

        dict_of_fasta[name] = sequence

    return dict_of_fasta

def retrieve_first_fasta_seq_from_fi(fasta_file):
    """ Get sequence for first or only fasta entry in a fasta file

    Args:
        fasta_file: input fasta file to retrieve sequence from
    
    Return:
        first sequence in given fasta file
    """
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        break
    return sequence    




