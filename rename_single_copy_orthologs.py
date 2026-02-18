#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import pandas as pd
from Bio import SeqIO

def _get_args():
    parser = argparse.ArgumentParser(description='Rename single-copy orthologue sequences identified by OrthoFinder2')
    parser.add_argument("--input","-i", "--in",   metavar="DIR", help="Single_Copy_Orthologue_Sequences directory", required=True)
    parser.add_argument("--output","-o", "--out",metavar="DIR",help="output directory", required=True)
    parser.add_argument("-t", "--orthogroups", type=str ,help="Orthogroups.tsv", required=True)
    parser.add_argument("-s", "--single_copy_orthologues", type=str ,help="Orthogroups_SingleCopyOrthologues.txt", required=True)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def rename_fa(in_dir, out_dir, og_id, og_table):
    in_fa  = '{}/{}.fa'.format(in_dir,og_id)
    out_fa = '{}/{}.fa'.format(out_dir,og_id)
    with open(in_fa, 'r') as infh, open(out_fa, 'w') as outfh :
        records = SeqIO.parse(infh,'fasta')
        for r in records:
            r_name = r.id
            sp_name = og_table.columns[(og_table.query(f'Orthogroup == "{og_id}"').apply(lambda x: x.str.contains(rf'\b{r_name}\b', na=False)).any())][0]
            r.id = sp_name
            r.description = ''
            SeqIO.write(r, outfh, 'fasta')

def main():
    args = _get_args()
    in_dir = args.input
    out_dir = args.output
    orthogroups = args.orthogroups
    single_copy_orthologues = args.single_copy_orthologues
    df_og = pd.read_table(orthogroups)
    with open(single_copy_orthologues) as sog:
        OGs = sog.read().splitlines()
        for OG in OGs:
            rename_fa(in_dir, out_dir, OG, df_og)

if __name__ == '__main__':
    main()
