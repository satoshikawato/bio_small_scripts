#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def _get_args():
    parser = argparse.ArgumentParser(
        description='Split scaffolds into contigs')
    parser.add_argument(
        "--input",
        "-i",
        "--in",
        type=str,
        help="Input FASTA file",
        required=True)
    parser.add_argument(
        "--output",
        "-o",
        "--out",
        type=str,
        help="output FASTA file (default: stdout)")
    parser.add_argument(
        "-d",
        "--digit",
        type=int,
        default=3,
        help="number of digits for zero-padding (default:3)")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def split_scaffolds(in_records, digit):
    out_records = {}
    for record in in_records:
        record_id = record.id
        contig_list = re.sub('[nN]+','\n',str(record.seq)).split('\n')
        for i in range(1,len(contig_list)):
            contig_id = "{}_{}".format(record_id, str(i).zfill(digit))
            out_records[contig_id] = contig_list[i]
    return out_records

def output_contigs(out_records, out_fa):
    if out_fa:
        with open(out_fa, "w") as f:
            for key in out_records.keys():
                print(">{}\n".format(key), file=f, end='')
                print("{}\n".format(out_records[key]), file=f, end='')
    else:
        for key in out_records.keys():
            print(">{}\n".format(key), end='')
            print("{}\n".format(out_records[key]),end='')

def main():
    args = _get_args()
    in_fa = args.input
    in_records = SeqIO.parse(in_fa, "fasta")
    out_fa = args.output
    digit = args.digit
    out_records = split_scaffolds(in_records, digit)
    output_contigs(out_records, out_fa)

if __name__ == "__main__":
    main()
