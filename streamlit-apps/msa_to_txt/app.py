#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
import os
import traceback
from collections import defaultdict
from collections import Counter
from contextlib import redirect_stdout
from Bio import AlignIO
import streamlit as st
from io import StringIO
def parse_arguments():
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
    ref_record_len_gap_exclusive = len(ref_record.seq.replace(gap_character, ""))
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

def residue_to_column_map(ref_record, gap_character, gap_inclusive):
    """
    ref_record の座標を alignment のカラム番号にマッピング
    """
    residue_to_col = {}
    col = 0
    pos = 1  # refの1-based座標
    for base in ref_record.seq:
        if gap_inclusive or base != gap_character:
            residue_to_col[pos] = col
            pos += 1
        col += 1
    return residue_to_col

def define_preceding_residues(
        records,
        ref_name,
        ref_record,
        start,
        end,
        gap_character,
        gap_inclusive):
    preceding_residues = defaultdict(int)
    if ref_name is not None:
        start, end = check_start_end_coords(
            ref_record, start, end, gap_character)
        residue_to_col = residue_to_column_map(ref_record, gap_character, gap_inclusive)
        aln_start = residue_to_col[start]
        aln_end = residue_to_col[end]
        aln_out = records[:, aln_start:aln_end + 1]
        aln_preceding = records[:, :aln_start]
    else:
        aln_out = records
        aln_preceding = records[:, 0:0]

    for record in aln_preceding:
        preceding_residues[record.id] = len(record.seq.replace(gap_character, "")) + 1
    return aln_out, preceding_residues

def get_records_text(records, start, end, wrap, preceding_residues, gap_character):
    buffer = StringIO()
    count = 0
    residues = get_residues(records)

    residues_chunks = [residues[i:i + wrap]
                       for i in range(0, len(residues), wrap)]
    record_chunks = [records[:, i:i + wrap]
                     for i in range(0, len(residues), wrap)]
    num_chunks = len(residues_chunks)

    chunk_seq_len = defaultdict(int)
    for key in preceding_residues.keys():
        chunk_seq_len[key] = preceding_residues[key]

    for chunk in range(num_chunks):

        record_chunk = record_chunks[chunk]
        for record in record_chunk:
            count += 1
            chunk_seq_len[record.id] += len(record.seq.replace(gap_character, ""))
            if len(record.seq.replace(gap_character, "")) == 0:
                start_pos = int(chunk_seq_len[record.id] - 1)
                end_pos = start_pos
            else:
                start_pos = (chunk_seq_len[record.id] - len(record.seq.replace(gap_character, "")))
                end_pos = (chunk_seq_len[record.id] - 1)
            record_line = make_text(record, start_pos, end_pos)
            buffer.write(record_line + "\n")
        residues_line = add_conservation(residues_chunks[chunk])
        buffer.write(residues_line + "\n")
        buffer.write("\n")
    return buffer.getvalue()

def main(raw_args=None):
    args = parse_arguments(raw_args)
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


st.title("🧬 MSA to TXT Converter")

# === Upload ===
uploaded_file = st.file_uploader("Upload MSA (FASTA format)", type=["fa", "fasta"])
if uploaded_file:
    st.write("✅ File uploaded:")
    st.code(uploaded_file.name)

# === Options ===
ref_name = st.text_input("Reference sequence name (optional)")
output_filename = st.text_input("Output file name (optional)")
start_str = st.text_input("Start position (leave blank for 0)", value="")
end_str = st.text_input("End position (leave blank for end)", value="")
gap_character = st.text_input("Gap character", value="-")
wrap = st.number_input("Line width", min_value=10, max_value=200, value=100)
gap_inclusive = st.checkbox("Gap inclusive counting", value=False)
    

# === Convert ===
if uploaded_file and st.button("Convert"):
    try:
        if output_filename:
            output_filename = output_filename
        else:
            output_basename = os.path.splitext(uploaded_file.name)[0]
            output_filename = f'{output_basename}.txt'
        # Read 
        input_text = uploaded_file.getvalue().decode("utf-8")

        records = AlignIO.read(StringIO(input_text), "fasta")
        if ref_name.strip():
            if ref_name.strip() not in [r.id for r in records]:
                st.error(f"❌ Reference sequence '{ref_name.strip()}' not found in MSA.")
                st.stop()
            ref_record = get_ref_record(records, ref_name.strip())
        else:
            ref_record = records[0]
        try:
            start = int(start_str) if start_str.strip() else 0
        except ValueError:
            st.error("❌ Start position must be an integer.")
            st.stop()

        try:
            end = int(end_str) if end_str.strip() else len(ref_record.seq)
        except ValueError:
            st.error("❌ End position must be an integer.")
            st.stop()

        aln_out, preceding_residues = define_preceding_residues(
            records, ref_name, ref_record, start, end, gap_character, gap_inclusive)

        # MSA proccessing
        txt_output = get_records_text(aln_out, start, end, wrap, preceding_residues, gap_character)

        # Output
        st.success("✅ Conversion complete!")
        st.code(txt_output, language=None)  
        st.download_button(
            label=f"Download {output_filename}",
            data=txt_output,
            file_name=output_filename,
            mime="text/plain"
        )

    except Exception as e:
        st.error(f"❌ {type(e).__name__}: {e if e else '(no message)'}")
        st.text(traceback.format_exc())

st.markdown("---")
st.markdown(
    "Author: [Satoshi Kawato](https://github.com/satoshikawato)  |  "
    "Source: [bio_small_scripts](https://github.com/satoshikawato/bio_small_scripts)",
    unsafe_allow_html=True
)