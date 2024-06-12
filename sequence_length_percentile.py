#!/usr/bin/env python
# coding: utf-8

import sys
import dnaio
import os
import argparse
from numpy import array, percentile


def _get_args(raw_args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Get the selected percentiles of sequence lengths in a FASTA or FASTQ file")
    parser.add_argument("-i", "--input", required=True, type=str, help="Input FASTA or FASTQ file")
    parser.add_argument("-p", "--percentile", help="Comma-separated list of percentiles to calculate (e.g. '5.0,95.0')", type=str, required=True)
    parser.add_argument("-o", "--output", dest="output", default=None, type=str, metavar="path", help="Output file name (default: stdout)")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(raw_args)

def calculate_lengths(records):
    length_array = array([len(record.sequence) for record in records])
    return length_array

def process_file(file_path):
    if file_path == "-":
        input_file = dnaio.open(sys.stdin.buffer)
    else:
        input_file = dnaio.open(file_path)
    records = list(input_file)
    return records

def main(raw_args=None):
    args = _get_args(raw_args)
    in_fa = args.input
    percentiles = [float(p) for p in args.percentile.split(",")]
    records = process_file(in_fa)
    length_array = calculate_lengths(records)
    for pctl in percentiles:
        print(f"{pctl}%: {int(percentile(length_array, pctl))}")

if __name__ == "__main__":
    main()