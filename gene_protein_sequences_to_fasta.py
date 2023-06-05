#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(
        description='Extract protein sequences from genbank files and write them into separate files by gene name')
    parser.add_argument(
        '-i',
        '--input',
        help='GenBank flat file format of the genomic sequence(s) (required)',
        type=str,
        required=True,nargs='*')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def load_gbks(gbk_list):
    """Load GenBank records as SeqRecords.

    Args:
        gbk_list (list): List of GenBank flat file paths.

    Returns:
        list: List of SeqRecord objects.
    """
    record_list = []
    for gbk_file in gbk_list:
        records = SeqIO.parse(gbk_file, 'genbank')
        for record in records:
            record_list.append(record)
    return record_list

def main():
    args = get_args()
    gb_files = args.input
    out_records = defaultdict(list)
    records = load_gbks(gb_files)
    for record in records:
        record_id = record.annotations['organism'].replace(" ", "_")
        for feature in record.features:
            if feature.type == 'CDS' and 'gene' in feature.qualifiers and 'translation' in feature.qualifiers:
                gene = feature.qualifiers['gene'][0]
                out_record = SeqRecord(Seq(feature.qualifiers['translation'][0]), id=record_id,description="")
                out_records[gene].append(out_record)
    for key, records in out_records.items():
        out_faa = "{}.faa".format(key)
        with open(out_faa, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

if __name__ == "__main__":
    main()
