#!/usr/bin/env python
# coding: utf-8

import sys
import os
import re
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from contextlib import redirect_stdout


def _get_args():
    parser = argparse.ArgumentParser(
        description='Extract the longest isoforms of protein-coding genes from a NCBI RefSeq euaryotic genome assembly')
    parser.add_argument(
        '-i', '--input', help='GenBank flat file format of the genomic sequence(s) (required)', type=str, required=True)
    parser.add_argument(
        '-o', '--output', help='output bed file', type=str, default="out.bed")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def _gb_record_to_bed(gb_record):
    for feature in gb_record.features:
        if feature.type == 'gene':
            if feature.location.strand == 1:
                strand  = "+"
            else:
                strand  = "-"
            start = re.sub(r'[^0-9]', '', str(feature.location.start))
            end = re.sub(r'[^0-9]', '', str(feature.location.end))
            if int(start) <= 1: 
                start ==  0
            else:
                 start == feature.location.start - 1
            print('{}\t{}\t{}\t{}\t{}\t{}'.format(gb_record.id, start, end, strand, "1000", feature.qualifiers['gene'][0]))

def main():
    args = _get_args()
    in_gb = args.input
    out_file = args.output
    gb_records = SeqIO.parse(in_gb, 'genbank')
    with open(out_file, "w") as output_table:
        with redirect_stdout(output_table):
            for record in gb_records:
                _gb_record_to_bed(record)


if __name__ == "__main__":
    main()
