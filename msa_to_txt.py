#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
from contextlib import redirect_stdout


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
        help="output txt file",
        required=True)
    parser.add_argument(
        "-r",
        "--ref",
        type=str,
        help="reference entry name",
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
    parser.add_argument(
        "-g",
        "--gap",
        type=str,
        help='gap character (default: "-")',
        default="-")
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
    if start < 1:
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
            gap_exclusive_as_key_inclusive_value[gap_exclusive_residue_count] = i
            gap_exclusive_residue_count += 1
    return gap_exclusive_as_key_inclusive_value


def add_record_to_track(record, start_end_dict):

    start_x = 0
    start_y = 11
    record_id = str(record.id)[:15]
    record_id = record_id.ljust(15, " ")
    out_text = "{} {} {}{}".format(record_id,
                                   str(start_end_dict[record.id][1]).rjust(4),
                                   record.seq,
                                   str(start_end_dict[record.id][2]).rjust(4))
    return out_text


def create_residue_dir(records):
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
    return residue_dir


def add_conservation_to_track(residue_dir):
    residues = []
    for key in residue_dir.keys():
        residue = residue_dir[key]
        residues.append(residue)
    residues = "".join(residues)
    record_id = " ".ljust(15, " ")
    out_text = "{} {} {}{}".format(
        record_id, " ".rjust(4), residues, " ".rjust(4))
    return out_text


def add_records_on_canvas(records, start_end_dict):
    count = 0
    residue_dir = create_residue_dir(records)
    residue_group = add_conservation_to_track(residue_dir)
    print(residue_group)
    for record in records:
        count += 1
        record_group = add_record_to_track(record, start_end_dict)
        print(record_group)


def make_start_end_dict(records, residue_dict, start, end):
    record_ids = []
    for record in records:
        record_ids.append(len(record.id))
    max_name_length = max(record_ids)
    start_end_dict = {}
    for record in records:
        name_space = max_name_length - len(record.id)
        start_loc = len(record.seq[:residue_dict[start]].ungap()) + 1
        end_loc = len(record.seq[:residue_dict[end] + 1].ungap())
        start_end_dict[record.id] = (name_space, start_loc, end_loc)
    return start_end_dict


def main():
    args = _get_args()
    in_fa = args.input
    out_file = args.output
    ref_name = args.ref
    start = args.start
    end = args.end
    gap_character = args.gap
    records = AlignIO.read(in_fa, "fasta")
    out_file_prefix = "test"
    num_of_entries = len(records)
    ref_record = get_ref_record(records, ref_name)
    start, end = check_start_end_coords(ref_record, start, end, gap_character)
    residue_dict = residue_count(ref_record, gap_character)
    aln_out = records[:, residue_dict[start]:residue_dict[end] + 1]
    start_end_dict = make_start_end_dict(records, residue_dict, start, end)
    with open(out_file, 'w') as f:
        with redirect_stdout(f):
            add_records_on_canvas(aln_out, start_end_dict)


if __name__ == "__main__":
    main()
