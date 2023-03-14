#!/usr/bin/env python
# coding: utf-8

import sys,os
import pysam
import argparse
from contextlib import redirect_stdout

def _get_args():
    parser = argparse.ArgumentParser(description='Extract flanking reads from a BAM file')
    parser.add_argument("--input","-i", "--in",metavar="FILE",help="Input BAM file",required=True)
    parser.add_argument("--output","-o","--out", "--output", metavar="FILE", help="output txt file")
    parser.add_argument("-r","--ref",type=str,help="reference entry name",required=True)
    parser.add_argument("-w","--window",type=int,help="window (defalt:100)",default=100)
    parser.add_argument("-m","--min",type=int,help="minimum outut read lenth threshold (defalt:200)",default=200)
    parser.add_argument("-f","--flank",type=int,help="flanking bases (defalt:100)",default=100)
    parser.add_argument("-p","--prime",type=int,help="5'/3'-end",default=5)
    parser.add_argument("-q","--fastq",help="output fastq",action='store_true')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

class ReadObject:
    def __init__(self, read_name, read_sequence, read_quality):
        self.name = read_name
        self.sequence = read_sequence
        self.quality = read_quality

def calculate_read_coordinates(read, window, flank, prime):
    if prime == 5:
        read_start = read.query_alignment_start - flank
        read_end = read.query_alignment_start + window 
        read_length = len(read.query_sequence[read_start:read_end])
    else:
        read_start = read.query_alignment_end - window
        read_end = read.query_alignment_end + flank 
        read_length = len(read.query_sequence[read_start:read_end])
    return read_start, read_end, read_length
    
def create_read_object(read, read_start, read_end):
    new_read_name = "{}:{}-{}".format(read.qname, read_start, read_end)
    new_read_sequence = read.query_sequence[read_start:read_end]
    new_read_quality = str(''.join(map(lambda x: chr( x+33 ), read.query_qualities[read_start:read_end])))
    read_object = ReadObject(new_read_name, new_read_sequence, new_read_quality)
    return read_object

def print_fastq_entry(read):
    print('@{}'.format(read.name))
    print(read.sequence)
    print("+")
    print(read.quality)
    
def print_fasta_entry(read):
    print('>{}'.format(read.name))
    print(read.sequence)

def main():
    args = _get_args()
    in_bam = args.input
    out_file = args.output
    flank= args.flank
    window = args.window
    read_length_threshold = args.min
    prime=args.prime
    ref=args.ref
    fastq = args.fastq
    bamfile = pysam.AlignmentFile(in_bam, "rb")
    reference_length = bamfile.get_reference_length(ref)
    if prime == 5:
        start = 1
        end = window
    else: 
        start = reference_length - window
        end = reference_length
    reads = list()
    for read in bamfile.fetch(ref, start, end):
        read_start, read_end, read_length = calculate_read_coordinates(read, window, flank, prime)
        if read_length < read_length_threshold:
            continue
        read_object = create_read_object(read, read_start, read_end)
        reads.append(read_object)
    if out_file is not None:
        with open(out_file, 'w') as f:
            with redirect_stdout(f):
                for read in reads:
                    if fastq == True: 
                        print_fastq_entry(read)
                    else: 
                        print_fasta_entry(read)
    else:
        for read in reads:
            if fastq == True: 
                print_fastq_entry(read)
            else: 
                print_fasta_entry(read)

if __name__ == "__main__":
    main()
