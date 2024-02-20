#!/usr/bin/env python
# coding: utf-8

import argparse
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import GC
import gzip
import logging
import time

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract and filter FASTA/FASTQ entries by various criteria")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA/FASTQ file (optionally gzipped)")
    parser.add_argument("-f", "--id_file", help="File containing list of IDs to include")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA/FASTQ file (gzipped if ending with .gz)")
    parser.add_argument("-m", "--min_length", type=int, help="Minimum sequence length to include")
    parser.add_argument("-M", "--max_length", type=int, help="Maximum sequence length to include")
    parser.add_argument("--gc_min", type=float, help="Minimum GC content to include (as a percentage)")
    parser.add_argument("--gc_max", type=float, help="Maximum GC content to include (as a percentage)")
    parser.add_argument("-q", "--qual_min", type=int, help="Minimum average quality score to include")
    parser.add_argument("-Q", "--qual_max", type=int, help="Maximum average quality score to include")
    parser.add_argument("-c", "--compress_level", type=int, default=3, help="Gzip compression level (1-9)")
    return parser.parse_args()

def open_file(file_path, mode, compresslevel):
    if file_path.endswith('.gz'):
        if 'w' in mode:
            return gzip.open(file_path, mode + 't', compresslevel=compresslevel)  # Add 't' to ensure text mode
        else:
            return gzip.open(file_path, mode + 't')  # For reading, compresslevel is not needed
    else:
        return open(file_path, mode)

def average_quality(qualities):
    return sum([ord(q) - 33 for q in qualities]) / len(qualities)

def filter_sequences(input_file, output_file, id_file=None, min_length=None, max_length=None, gc_min=None, gc_max=None, qual_min=None, compresslevel=3):
    processed_sequences = 0
    start_time = time.time()
    output_fasta = output_file.endswith('.fa') or output_file.endswith('.fasta')
    with open_file(input_file, 'rt', compresslevel) as infile, open_file(output_file, 'wt', compresslevel=compresslevel) as outfile:
        if input_file.endswith('.fastq') or input_file.endswith('.fq') or input_file.endswith('.fq.gz') or input_file.endswith('.fastq.gz'):
            for title, seq, qual in FastqGeneralIterator(infile):
                processed_sequences += 1
                if processed_sequences % 10000 == 0:  # Log progress every 10000 sequences
                    logging.info(f"Processed {processed_sequences} sequences") 
                if id_file and title.split()[0] not in id_file:
                    continue
                if min_length and len(seq) < min_length:
                    continue
                if max_length and len(seq) > max_length:
                    continue
                if gc_min and GC(seq) < gc_min:
                    continue
                if gc_max and GC(seq) > gc_max:
                    continue
                if qual_min and average_quality(qual) < qual_min:
                    continue
                if output_fasta:
                    # Write in FASTA format if output is .fa or .fasta
                    outfile.write(f">{title}\n{seq}\n")
                else:
                    # Write in FASTQ format if output is not .fa or .fasta
                    outfile.write(f"@{title}\n{seq}\n+\n{qual}\n")
               
        else:  # Assume FASTA format
            for title, seq in SimpleFastaParser(infile):
                processed_sequences += 1
                if processed_sequences % 10000 == 0:  # Log progress every 10000 sequences
                    logging.info(f"Processed {processed_sequences} sequences") 
                if id_file and title.split()[0] not in id_file:
                    continue
                if min_length and len(seq) < min_length:
                    continue
                if max_length and len(seq) > max_length:
                    continue
                if gc_min and GC(seq) < gc_min:
                    continue
                if gc_max and GC(seq) > gc_max:
                    continue
                outfile.write(f">{title}\n{seq}\n")
    end_time = time.time()
def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    args = parse_arguments()
    id_set = set()
    if args.id_file:
        with open(args.id_file) as f:
            for line in f:
                id_set.add(line.strip())
    filter_sequences(args.input, args.output, id_set if args.id_file else None, args.min_length, args.max_length, args.gc_min, args.gc_max, args.qual_min, compresslevel=args.compress_level)

if __name__ == "__main__":
    main()
