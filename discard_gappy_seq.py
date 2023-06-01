#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def _get_args():
    parser = argparse.ArgumentParser(
        description='Discard FASTA entries with gaps')
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
    args = parser.parse_args()
    return args

def discard_gaps(in_records):
    out_records_tmp = list()
    out_records_dict = dict()
    out_records = list()
    for record in in_records:
        if record.seq.count("N") >0:
            continue
        else:
            out_record = SeqRecord(record.seq, id=record.id, description=record.description)
            out_records_tmp.append(out_record)
    for record in out_records_tmp:
        out_records_dict[record.id] = record
    out_records_tmp = sorted(out_records_dict.items())
    for record in out_records_tmp:
        out_records.append(record[1])
    return out_records

def output_fasta(out_fa, out_records):
    with open(out_fa, "w") as output_handle:
        for out_record in out_records:
            SeqIO.write(out_record, output_handle, "fasta")

def main():
    args = _get_args()
    in_fa = args.input
    out_fa = args.output
    in_records = SeqIO.parse(in_fa, "fasta")
    out_records = discard_gaps(in_records)
    output_fasta(out_fa, out_records)


if __name__ == "__main__":
    main()
