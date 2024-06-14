#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys
import gzip
from Bio import SeqIO

def parse_arguments(raw_args=None):
    parser = argparse.ArgumentParser(description="Extract and filter FASTA/FASTQ entries by various criteria")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA/FASTQ file (optionally gzipped)")
    parser.add_argument("-o", "--output", required=True, help="Output directory for split FASTA/FASTQ files")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(raw_args)

def open_input_file(input_file):
    if input_file.endswith('.gz'):
        return gzip.open(input_file, 'rt')
    else:
        return open(input_file, 'r')

def split_fasta(input_fasta, output_dir):
    fasta_records = SeqIO.parse(input_fasta, 'fasta')
    os.makedirs(output_dir, exist_ok=True)

    for record in fasta_records:
        output_file = os.path.join(output_dir, f'{record.id}.fa')
        with open(output_file, 'w') as out_fasta:
            SeqIO.write(record, out_fasta, 'fasta')

def split_fastq(input_fastq, output_dir):
    fastq_records = SeqIO.parse(input_fastq, 'fastq')
    os.makedirs(output_dir, exist_ok=True)

    for record in fastq_records:
        output_file = os.path.join(output_dir, f'{record.id}.fq')
        with open(output_file, 'w') as out_fastq:
            SeqIO.write(record, out_fastq, 'fastq')

def main(raw_args=None):
    args = parse_arguments(raw_args)

    with open_input_file(args.input) as input_file:
        if args.input.endswith('.fa') or args.input.endswith('.fasta') or args.input.endswith('.fa.gz') or args.input.endswith('.fasta.gz'):
            split_fasta(input_file, args.output)
        elif args.input.endswith('.fq') or args.input.endswith('.fastq') or args.input.endswith('.fq.gz') or args.input.endswith('.fastq.gz'):
            split_fastq(input_file, args.output)
        else:
            sys.stderr.write("Input file must be in FASTA or FASTQ format (optionally gzipped).\n")
            sys.exit(1)

if __name__ == "__main__":
    main()
