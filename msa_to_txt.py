#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
from collections import defaultdict
from collections import Counter
from contextlib import redirect_stdout
from Bio import AlignIO

def _get_args():
    parser = argparse.ArgumentParser(
        description='Convert FASTA-format multiple sequence alignment into a txt file. Assumes Courier New')
    parser.add_argument(
        "--input",
        "-i",
        "--in",
        metavar="FILE",
        help="Input FASTA file",
        required=True)
    parser.add_argument(
        "--output",
        "-o",
        "--out",
        "--output",
        metavar="FILE",
        help="output txt file")
    parser.add_argument(
        "-r",
        "--ref",
        type=str,
        help="reference entry name")
    parser.add_argument(
        "-s",
        "--start",
        type=int,
        help="start position")
    parser.add_argument(
        "-e",
        "--end",
        type=int,
        help="end position")
    parser.add_argument(
        "-g",
        "--gap",
        type=str,
        help='gap character (default: "-")',
        default="-")
    parser.add_argument(
        "-w",
        "--wrap",
        type=int,
        help='line width (default: 100)',
        default=100)
    parser.add_argument(
        '--gap_inclusive',
        help='Gap inclusive (default: False). ',
        action='store_true')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args


def get_ref_record(records, ref_name):
    if ref_name is not None:
        record_dict = dict()
        for record in records:
            record_dict[record.id] = record
        ref_record = record_dict[ref_name]
    else:
        ref_record = records[0]
    return ref_record


def check_start_end_coords(
        ref_record,
        start,
        end,
        gap_character):
    ref_record_len_gap_exclusive = len(ref_record.seq.ungap(gap_character))
    if start < 1:
        start = 1
    if end > ref_record_len_gap_exclusive:
        end = ref_record_len_gap_exclusive
    return start, end

def residue_count(ref_record, gap_character, gap_inclusive):
    gap_character = gap_character
    residue_count = 1
    residue_count_dict = dict()
    seq_length = len(ref_record.seq)
    if gap_inclusive is True:
        for i in range(seq_length):
            residue_count_dict[residue_count] = i
            residue_count += 1
    else:
        for i in range(seq_length):
            if ref_record.seq[i] == gap_character:
                pass
            else:
                residue_count_dict[residue_count] = i
                residue_count += 1
    return residue_count_dict

def get_residues(records):
    alignment_length = int(records.get_alignment_length())
    residue_dir = {}
    num_records = len(records)
    for column in range(0, alignment_length):
        residue_list = []
        column_residues = records[:, column:column + 1]
        for column_residue in column_residues:
            residue_list.append(str(column_residue.seq))
        c = Counter(residue_list)
        counter_out = [(i, c[i] / num_records * 100.0)
                       for i, count in c.most_common()]
        if counter_out[0][0] == "-":
            residue_dir[column] = " "
        else:
            if 80 <= float(counter_out[0][1]) < 100:
                residue_dir[column] = "."
            elif float(counter_out[0][1]) == 100:
                residue_dir[column] = "*"
            else:
                residue_dir[column] = " "
    residues = []
    for key in residue_dir.keys():
        residue = residue_dir[key]
        residues.append(residue)
    return residues

def add_conservation(residues):
    residues = "".join(residues)
    record_id = "".ljust(15, " ")
    out_text = "{} {} {} {}".format(
        record_id, " ".rjust(4), residues, " ".rjust(4))
    return out_text

def make_text(record, start, end):
    record_id = str(record.id)[:15].ljust(15, " ")
    out_text = "{} {} {} {}".format(
        record_id,
        str(start).rjust(4),
        record.seq,
        str(end).rjust(4))
    return out_text

def make_start_end_dict(records, residue_dict, start, end):
    record_ids = []
    for record in records:
        record_ids.append(len(record.id))
    max_name_length = max(record_ids)
    start_end_dict = {}
    for record in records:
        name_space = max_name_length - len(record.id)
        start_loc = len(record.seq[:residue_dict[start]].replace(gap_character, "")) + 1
        end_loc = len(record.seq[:residue_dict[end] + 1].replace(gap_character, ""))
        start_end_dict[record.id] = (name_space, start_loc, end_loc)
    return start_end_dict

def define_preceding_residues(
        records,
        ref_name,
        ref_record,
        start,
        end,
        gap_character,
        gap_inclusive):
    preceding_residues = defaultdict(dict)
    if ref_name is not None:
        start, end = check_start_end_coords(
            ref_record, start, end, gap_character, gap_inclusive)
        residue_dict = residue_count(ref_record, gap_character, gap_inclusive)
        aln_out = records[:, residue_dict[start]:residue_dict[end] + 1]
        aln_preceding = records[:, :residue_dict[start] - 1]
    else:
        start, end = 1, len(ref_record.seq)
        aln_preceding = records[:, 0:0]
        aln_out = records
    for record in aln_preceding:
        if ref_name is not None:
            preceding_residues[record.id] = (len(record.seq.replace(gap_character, "")) + 2)
        else:
            preceding_residues[record.id] = (len(record.seq.replace(gap_character, "")) + 1)
    return aln_out, preceding_residues

def print_records(records, start, end, wrap, preceding_residues, gap_character):
    count = 0
    residues = get_residues(records)

    residues_chunks = [residues[i:i + wrap]
                       for i in range(0, len(residues), wrap)]
    record_chunks = [records[:, i:i + wrap]
                     for i in range(0, len(residues), wrap)]
    num_chunks = len(range(0, len(residues), wrap))

    chunk_seq_len = defaultdict(int)
    for key in preceding_residues.keys():
        chunk_seq_len[key] = preceding_residues[key]
    for chunk in range(0, num_chunks):
        residues = residues_chunks[chunk]
        residue_group = add_conservation(residues)
        print(residue_group)
        record_chunk = record_chunks[chunk]
        for record in record_chunk:
            count += 1
            chunk_seq_len[record.id] += len(record.seq.replace(gap_character, ""))
            start = (chunk_seq_len[record.id] - len(record.seq.replace(gap_character, "")))
            if len(record.seq.replace(gap_character, "")) == 0:
                end = start
            else:
                end = (chunk_seq_len[record.id] - 1)
            record_group = make_text(record, start, end)
            print(record_group)
        print()

def main():
    args = _get_args()
    in_fa = args.input
    out_file = args.output
    ref_name = args.ref
    start = args.start
    end = args.end
    wrap = args.wrap
    gap_character = args.gap
    gap_inclusive = args.gap_inclusive
    records = AlignIO.read(in_fa, "fasta")
    ref_record = get_ref_record(records, ref_name)

    if start is None:
        start = 1
    if end is None or end > len(ref_record.seq):
        end = len(ref_record.seq)
    aln_out, preceding_residues = define_preceding_residues(
        records, ref_name, ref_record, start, end, gap_character, gap_inclusive)

    if out_file is not None:
        with open(out_file, 'w') as f:
            with redirect_stdout(f):
                print_records(aln_out, start, end, wrap, preceding_residues, gap_character)
    else:
        print_records(aln_out, start, end, wrap, preceding_residues, gap_character)

if __name__ == "__main__":
    main()
