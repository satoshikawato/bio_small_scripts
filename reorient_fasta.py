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
        description='Reorient FASTA entries based on BLAST results')
    parser.add_argument(
        '-i',
        '--input',
        help='FASTA file (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-b',
        '--blast',
        help='Tab-separated BLAST output file (required)',
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

def blast_to_dict(blast_out):
    '''
    Read a tab-delimited blast output (-outfmt 6 or 7).
    Comment lines are ignored.
    Returns a dictionary (key: query; value: whole line, tab-delimited).
    '''
    blast_out_dict = defaultdict(list)
    with open(blast_out) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            else:
                line = line.split('\t')
                contig_id = line[0]
                blast_out_dict[contig_id].append(line)
    return blast_out_dict

def reorient_fasta(in_records, blast_out_dict):
    out_records_tmp = list()
    out_records_dict = dict()
    out_records = list()
    for record in in_records:
        line = blast_out_dict[record.id][0]
        if int(line[6]) > int(line[7]):
            out_record = SeqRecord(record.seq.reverse_complement(), id=line[1], description=record.description)
            out_records_tmp.append(out_record)
        else:
            out_record = SeqRecord(record.seq, id=line[1], description=record.description)
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
    in_blast = args.blast
    out_fa = args.output
    in_records = SeqIO.parse(in_fa, "fasta")
    blast_out_dict = blast_to_dict(in_blast)
    out_records = reorient_fasta(in_records, blast_out_dict)
    output_fasta(out_fa, out_records)

if __name__ == "__main__":
    main()
