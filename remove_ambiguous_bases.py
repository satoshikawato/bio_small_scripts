#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

FORMAT_FASTA = "fasta"

def get_arguments():
    parser = argparse.ArgumentParser(
        description='Discard DNA sequences with ambiguous bases')
    parser.add_argument(
        '-i',
        '--input',
        help='FASTA file (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        help='output fasta file (default: out.fa)',
        type=str,
        default="out.fa")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def filter_sequences(in_records):
    valid_bases = set('ATGC')
    out_records = []
    for record in in_records:
        if set(record.seq.upper()) <= valid_bases:
            out_records.append(SeqRecord(record.seq, id=record.id, description=record.description))
    return out_records

def write_output(out_fa, out_records):
    with open(out_fa, "w") as output_handle:
        SeqIO.write(out_records, output_handle, FORMAT_FASTA)

def main():
    args = get_arguments()
    in_records = SeqIO.parse(args.input, FORMAT_FASTA)
    out_records = filter_sequences(in_records)
    write_output(args.output, out_records)

if __name__ == "__main__":
    main()
