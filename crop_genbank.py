#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def _get_args():
    parser = argparse.ArgumentParser(description='Crop genbank file. ')
    parser.add_argument(
        '-i',
        '--input',
        help='Genbank/DDBJ flatfile (required)',
        type=str,
        required=True)
    parser.add_argument(
        "--output",
        "-o",
        "--out",
        "--output",
        metavar="FILE",
        help="output FASTA file",
        required=True)
    parser.add_argument(
        "-s",
        "--start",
        type=int,
        help="start position",
        required=True)
    parser.add_argument(
        "-e",
        "--end",
        type=int,
        help="end position",
        required=True)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args


def gbk_to_seqrecord(in_gbk):
    """Load the first instance of GenBank as SerRecord
    Args:
    """
    records = SeqIO.parse(gbk_file, 'genbank')
    record = next(records)
    return record


def check_start_end_coords(record, start, end):
    record_len = len(record.seq)
    start = start - 1
    if start < 0:
        start = 0
    if end > record_len:
        end = record_len
    return start, end


def main():
    args = _get_args()
    in_gbk = args.input
    out_gbk = args.output
    start = args.start
    end = args.end
    record = gbk_to_seqrecord(in_gbk)
    start, end = check_start_end_coords(record, start, end)
    new_record = record[start:end]
    SeqIO.write(new_record, out_gbk, "genbank")


if __name__ == "__main__":
    main()
