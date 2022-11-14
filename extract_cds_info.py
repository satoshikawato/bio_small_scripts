#!/usr/bin/env python
# coding: utf-8

import argparse
import collections
import logging
import pathlib
import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from contextlib import redirect_stdout

def _get_args():
    parser = argparse.ArgumentParser(
        description='Extract CDS information from genbank file')
    parser.add_argument(
        "--input",
        "-i",
        "--in",
        metavar="FILE",
        help="Input genbank file",
        required=True)
    parser.add_argument(
        "--output",
        "-o",
        "--out",
        "--output",
        metavar="FILE",
        help="output txt file")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    #args = parser.parse_args(args=['--input','CN01.gb',  '-o','out.txt'])
    return args

def print_record(records):
    for record in records:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format("start", "end", "strand", "Accession no.", "length", "product", "note"))
            for feature in record.features:
                if feature.type == 'CDS':
                    strand = ""
                    record_id = ""
                    protein_length =""
                    product = ""
                    note = ""
                    location = feature.location
                    if location.strand == 1:
                        strand = "+"
                    else: 
                        strand = "-"
                    if 'protein_id' in feature.qualifiers.keys():
                        record_id = feature.qualifiers['protein_id'][0]
                    
                    if 'translation' in feature.qualifiers.keys():
                        protein_length = len(feature.qualifiers['translation'][0])
                    if 'product' in feature.qualifiers.keys():            
                        product = feature.qualifiers['product'][0]
                    if 'note' in feature.qualifiers.keys():
                        note = feature.qualifiers['note'][0]
                    print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(location.start, location.end, strand, record_id, protein_length, product, note))


def main():
    args = _get_args()
    in_gb = args.input
    out_file = args.output
    records = SeqIO.parse(in_gb, "genbank")

    if out_file is not None:
        with open(out_file, 'w') as f:
            with redirect_stdout(f):
                 print_record(records)
    else:
        print_record(records)

if __name__ == "__main__":
    main()
