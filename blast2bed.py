#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from contextlib import redirect_stdout

def get_args():
    parser = argparse.ArgumentParser(
        description='convert BLASTN/BLASTX/TBLASTX output into BED format')
    parser.add_argument(
        '-i', '--input',
        help='tab-separated blast output (required) with "-outfmt 6"',
        type=str, required=True)
    parser.add_argument(
        '-o', '--output',
        help='output BED format file (default: stdout)',
        type=str)
    parser.add_argument(
        '-e', '--evalue',
        help='E-value threshold (default: 1e-30)',
        type=float, default=1e-30)
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def blast_to_bed(blast_out, evalue_threshold):
    with open(blast_out) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue
            
            cols = line.split('\t')
            if len(cols) < 12:
                continue

            # BLAST fmt 6 columns:
            # 0:qseqid, 1:sseqid, 2:pident, 3:length, 4:mismatch, 5:gapopen, 
            # 6:qstart, 7:qend, 8:sstart, 9:send, 10:evalue, 11:bitscore
            
            q_start = int(cols[6])
            q_end = int(cols[7])
            evalue = float(cols[10])
            bit_score = cols[11]
            subject_id = cols[1]
            chrom = cols[0]

            if evalue > evalue_threshold:
                continue

            # BED coordinates: 0-based start, 1-based end
            if q_start < q_end:
                start = q_start - 1
                end = q_end
                strand = "+"
            else:
                start = q_end - 1
                end = q_start
                strand = "-"

            # 1.chrom, 2.start, 3.end, 4.name, 5.score, 6.strand
            print(f"{chrom}\t{start}\t{end}\t{subject_id}\t{bit_score}\t{strand}")

def main():
    args = get_args()
    if args.output:
        with open(args.output, "w") as output_table:
            with redirect_stdout(output_table):
                blast_to_bed(args.input, args.evalue)
    else:
        blast_to_bed(args.input, args.evalue)
        
if __name__ == "__main__":
    main()