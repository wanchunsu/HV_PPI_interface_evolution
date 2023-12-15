#!/bin/bash
#process the file, so that we only have the relevant data (pdb, chain, taxid) and unique entries

cut -f 1-3 ../data/pdb_annotations/pdb_chain_taxonomy.tsv | uniq  > ../data/pdb_annotations/pdb_chain_taxonomy_rel.tsv