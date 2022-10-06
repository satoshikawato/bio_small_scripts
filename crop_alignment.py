#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def _get_args():
    parser = argparse.ArgumentParser(description='Crop multiple sequence alignment gap-exclusively.')
    parser.add_argument("--input","-i", "--in",   metavar="FILE", help="Input FASTA file", required=True)
    parser.add_argument("--output","-o", "--out", "--output",metavar="FILE",help="output FASTA file", required=True)
    parser.add_argument("-r", "--ref", type=str ,help="reference entry name", required=True)
    parser.add_argument("-s", "--start",  type=int,help="start position", required=True)
    parser.add_argument("-e", "--end",  type=int, help="end position", required=True)
    parser.add_argument("-g", "--gap", type=str, help="gap character (default: \"-\")", default="-")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def get_ref_record(records, ref_name):
    record_dict = dict()
    for record in records:
        record_dict[record.id] = record
    ref_record = record_dict[ref_name]
    return ref_record

def check_start_end_coords(ref_record, start, end, gap_character):
    ref_record_len_gap_exclusive = len(ref_record.seq.ungap(gap_character))
    if start < 1 : 
        start = 1
    if end > ref_record_len_gap_exclusive:
        end = ref_record_len_gap_exclusive 
    return start, end

def residue_count(ref_record, gap_character):
    gap_character = gap_character
    gap_exclusive_residue_count = 1
    gap_exclusive_as_key_inclusive_value = dict()
    seq_length = len(ref_record.seq)
    for i in range(seq_length):
        if ref_record.seq[i] == gap_character:
            pass
        else:
            gap_exclusive_as_key_inclusive_value[gap_exclusive_residue_count]  =  i
            gap_exclusive_residue_count += 1
    return gap_exclusive_as_key_inclusive_value

def main():
   args = _get_args()
   in_fa = args.input
   out_fa = args.output
   ref_name = args.ref
   start = args.start
   end = args.end
   gap_character = args.gap
   
   records = AlignIO.read(in_fa, "fasta")
   
   ref_record = get_ref_record(records, ref_name)
   start, end = check_start_end_coords(ref_record, start, end, gap_character)
   residue_dict = residue_count(ref_record, gap_character)
   aln_out = records[:, residue_dict[start]:residue_dict[end] + 1]
   AlignIO.write(aln_out, out_fa, "fasta")

main()
