#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
from contextlib import redirect_stdout

def get_args():
    parser = argparse.ArgumentParser(
        description='convert BLASTN/TBLASTN/TBLASTX output into BED format')
    parser.add_argument(
        '-i',
        '--input',
        help='tab-separated blast output (required) with "-outfmt "[6|7]"',
        type=str,
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        help='output BED format file (default: stdout)',
        type=str)
    parser.add_argument(
        '-s',
        '--score',
        help='score (default: 0)',
        type=int, default = 0)
    parser.add_argument(
        '-e',
        '--evalue',
        help='E-value threshold (default: 1e-30)',
        type=float, default = 1e-30)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def blast_to_bed(blast_out, score, evalue):
    with open(blast_out) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            else:
                line = line.split('\t')
                if float(line[10]) > evalue:
                    continue
                else:
                    if (float(line[8]) > float(line[9])):
                        start = str(int(line[9]) -1)
                        end = str(int(line[8]))
                        strand = "-"
                    else:
                        start = str(int(line[8])-1)
                        end = str(int(line[9]))
                        strand = "+"
                    print("{}\t{}\t{}\t{}\t{}\t{}".format(line[1], start, end, strand, score, line[0]))

def main():
    args = get_args()
    in_blast = args.input
    out_file = args.output
    score = args.score
    evalue = args.evalue
    if out_file:
        with open(out_file, "w") as output_table:
            with redirect_stdout(output_table):
                blast_to_bed(in_blast, score, evalue)
    else:
        blast_to_bed(in_blast, score, evalue)
        
if __name__ == "__main__":
    main()
