#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

def _get_args():
    parser = argparse.ArgumentParser(
        description='Extract protein sequences from a genbank file')
    parser.add_argument(
        '-i',
        '--input',
        help='GenBank flat file format of the genomic sequence(s) (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        help='output fasta file (default: out.faa)',
        type=str,
        default="out.faa")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def main():
    args = _get_args()
    in_gb = args.input
    out_faa = args.output
    records = SeqIO.parse(in_gb, 'genbank')
    out_records = []
    for record in records:
        for feature in record.features:
            if feature.type == 'CDS':
                if 'protein_id' in feature.qualifiers.keys():
                    record_id = feature.qualifiers['protein_id'][0]
                else: 
                    record_id = "{}_{}-{}".format(record.id,feature.location.start,feature.location.end)
                out_record = SeqRecord(
                Seq(
                        feature.qualifiers['translation'][0]),
                    id=record_id,
                    description="{} [{}]".format(
                        feature.qualifiers['product'][0],
                        record.annotations['organism']))
                out_records.append(out_record)
    with open(out_faa, "w") as output_handle:
        for out_record in out_records:
            SeqIO.write(out_record, output_handle, "fasta")

if __name__ == "__main__":
    main()
