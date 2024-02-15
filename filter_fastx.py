#!/usr/bin/env python
# coding: utf-8

import argparse
import gzip
from Bio import SeqIO
from Bio.SeqUtils import GC
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract and filter FASTA/FASTQ entries by various criteria")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA/FASTQ file (optionally gzipped)")
    parser.add_argument("-f", "--id_file", help="File containing list of IDs to include")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA/FASTQ file (gzipped if ending with .gz)")
    parser.add_argument("-m", "--min_length", type=int, help="Minimum sequence length to include")
    parser.add_argument("-M", "--max_length", type=int, help="Maximum sequence length to include")
    parser.add_argument("--gc_min", type=float, help="Minimum GC content to include (as a percentage)")
    parser.add_argument("--gc_max", type=float, help="Maximum GC content to include (as a percentage)")
    parser.add_argument("--q_min", type=float, help="Minimum average quality to include (FASTQ only)")
    parser.add_argument("--q_max", type=float, help="Maximum average quality to include (FASTQ only)")
    return parser.parse_args()

def seq_passes_filters(record, args):
    if args.desired_ids and not record.id in args.desired_ids:
        return False
    if args.min_length and len(record.seq) < args.min_length:
        return False
    if args.max_length and len(record.seq) > args.max_length:
        return False
    if args.gc_min and GC(record.seq) < args.gc_min:
        return False
    if args.gc_max and GC(record.seq) > args.gc_max:
        return False
    if args.q_min and sum(record.letter_annotations["phred_quality"]) / len(record) < args.q_min:
        return False
    if args.q_max and sum(record.letter_annotations["phred_quality"]) / len(record) > args.q_max:
        return False
    return True

def read_id_list(id_file):
    with open(id_file, 'r') as f:
        return {line.strip() for line in f}

def filter_sequences(input_file, output_file, args):
    format = "fasta" if input_file.endswith(".fasta") or input_file.endswith(".fa") else "fastq"
    opener = gzip.open if input_file.endswith(".gz") else open
    with opener(input_file, 'rt') as in_handle, opener(output_file, 'wt') as out_handle:
        records = (record for record in SeqIO.parse(in_handle, format) if seq_passes_filters(record, args))
        count = SeqIO.write(records, out_handle, format)
        print(f"Extracted {count} records")

def main():
    args = parse_arguments()

    # If an ID file is specified, update the desired_ids set
    if args.id_file:
        args.desired_ids = read_id_list(args.id_file)
    else:
        args.desired_ids = None

    filter_sequences(args.input, args.output, args)

if __name__ == "__main__":
    main()

