#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import streamlit as st
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from io import StringIO, BytesIO

def get_args():
    parser = argparse.ArgumentParser(description='Predict ORFs')
    parser.add_argument(
        '-i',
        '--input',
        help='sequence file in FASTA format (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-o',
        '--out_gff',
        help="output annotation in gff3 format (default: stdout)",
        type=str)
    parser.add_argument(
        '-a',
        '--out_faa',
        help="output protein sequences in FASTA format (optional)",
        type=str)
    parser.add_argument(
        '-f',
        '--out_fna',
        help="output CDS sequences in FASTA format  (optional)",
        type=str)
    parser.add_argument(
        '-g',
        '--trans_table',
        help='translation table (default: 1)',
        type=int,
        default='1')
    parser.add_argument(
        '-m',
        '--min_aa_len',
        help='minimum protein length (default: 50)',
        type=int,
        default='50')
    parser.add_argument(
        '--keep_nested',
        action='store_true',
        help='keep nested ORFs (default: False)')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return args

def fasta_to_records(in_fa):
    seq_records = [record for record in SeqIO.parse(in_fa, 'fasta')]
    return seq_records

def remove_nested_orfs(orf_list):
    if not orf_list:
        return []
    orf_list.sort(key=lambda x: (x[0], -x[1]))
    
    non_nested_orfs = []
    if not orf_list:
        return non_nested_orfs
        
    last_added = orf_list[0]
    non_nested_orfs.append(last_added)
    
    for i in range(1, len(orf_list)):
        current_orf = orf_list[i]
        if current_orf[1] <= last_added[1]:
            continue 
        else:
            non_nested_orfs.append(current_orf)
            last_added = current_orf
            
    return non_nested_orfs

def get_orf(record, trans_table, min_protein_length, keep_nested):
    orf_list = []
    prefix = record.id
    seq = record.seq
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            nuc_frame = nuc[frame:]
            if len(nuc_frame) % 3 != 0:
                nuc_frame = nuc_frame[:-(len(nuc_frame) % 3)]
            trans = nuc_frame.translate(trans_table)
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            aa_start = trans.find("M", aa_start)
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                        cds_seq = seq[start:end]
                    else:
                        start =  max(0, (seq_len - (frame + aa_end * 3 + 3)))
                        end = seq_len - (frame + aa_start * 3)
                        cds_seq = seq[start:end].reverse_complement()
                    orf_list.append(
                        [start, end, strand, trans[aa_start:aa_end], cds_seq])
                aa_start = trans.find("M", aa_end + 1)
                if aa_start < 0:
                    break
    if not keep_nested:
        orf_list = remove_nested_orfs(orf_list)
    
    orf_list.sort(key=lambda x: x[0])
    
    zfill_len = int(len(str(len(orf_list)))+1)
    count = 0
    for orf in orf_list:
        count += 1
        orf_id = "{}_{}".format(prefix, str(count).zfill(zfill_len))
        orf.insert(0, orf_id)
    return orf_list


def orf_to_faa(orf_dict):
    faa_records = []
    for key in orf_dict.keys():
        for orf in orf_dict[key]:
            orf_seq = orf[4]
            faa_record = SeqRecord(orf_seq, id=orf[0], name="", description="")
            faa_records.append(faa_record)
    return faa_records


def orf_to_fna(orf_dict):
    faa_records = []
    for key in orf_dict.keys():
        for orf in orf_dict[key]:
            orf_seq = orf[5]
            faa_record = SeqRecord(orf_seq, id=orf[0], name="", description="")
            faa_records.append(faa_record)
    return faa_records

def get_gff3_features(orf_dict):
    annot_lines_list = []
    for key in orf_dict.keys():
        for orf in orf_dict[key]:
            seqid = key
            source = "orfind.py"
            score = "."
            feature_type = "CDS"
            start = int(orf[1]) + 1
            end = int(orf[2])
            feature_strand = orf[3]
            if feature_strand == 1:
                strand = "+"
            else:
                strand = "-"
            phase = "0"
            feature_id = orf[0]
            qualifers = "ID={};".format(feature_id)
            annot_lines_list.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                seqid, source, feature_type, start, end, score, strand, phase, qualifers))
    
    return "".join(annot_lines_list) 


def main():
    args = get_args()
    records = fasta_to_records(args.input)
    out_gff3 = args.out_gff
    out_faa = args.out_faa
    out_fna = args.out_fna
    trans_table = args.trans_table
    min_protein_length = args.min_aa_len
    keep_nested = args.keep_nested
    header = "##gff-version  3"
    comments = ''
    comments += '{}\n'.format(header)
    orf_dict = {}

    for record in records:
        record_id = record.id
        orf_dict[record_id] = get_orf(record, trans_table, min_protein_length, keep_nested)
    gff3_features = get_gff3_features(orf_dict)
    if out_gff3:
        with open(out_gff3, "w") as f:
            print(comments, file=f, end='')
            print(gff3_features, file=f, end='')
    else:
        print(comments, end='')
        print(gff3_features, end='')
    if out_faa:
        faa_records = orf_to_faa(orf_dict)
        with open(out_faa, "w") as f:
            SeqIO.write(faa_records, f, "fasta")
    if out_fna:
        fna_records = orf_to_fna(orf_dict)
        with open(out_fna, "w") as f:
            SeqIO.write(fna_records, f, "fasta")

# === Predict and cache ORFs ===
def predict_orfs(fasta_text, trans_table, min_protein_length, keep_nested):
    stringio = StringIO(fasta_text)
    seq_records = list(SeqIO.parse(stringio, 'fasta'))

    orf_dict = {}
    for record in seq_records:
        orf_dict[record.id] = get_orf(record, trans_table, min_protein_length, keep_nested)

    gff3_str = get_gff3_features(orf_dict)
    gff3_bytes = BytesIO(gff3_str.encode())

    faa_records = orf_to_faa(orf_dict)
    faa_text = StringIO()
    SeqIO.write(faa_records, faa_text, "fasta")
    faa_bytes = BytesIO(faa_text.getvalue().encode("utf-8"))

    fna_records = orf_to_fna(orf_dict)
    fna_text = StringIO()
    SeqIO.write(fna_records, fna_text, "fasta")
    fna_bytes = BytesIO(fna_text.getvalue().encode("utf-8"))

    orf_count = sum(len(v) for v in orf_dict.values())

    return gff3_bytes, faa_bytes, fna_bytes, gff3_str, orf_count

# === Generate Codon Table Choices ===
codon_table_choices = [
    f"{table_id} ({table.names[0]})"
    for table_id, table in CodonTable.generic_by_id.items()
]
codon_table_ids = list(CodonTable.generic_by_id.keys())

# === Streamlit UI ===
st.title("🧬 ORFIND v0.1.0")

def reset_results():
    keys_to_remove = ['gff3_bytes', 'faa_bytes', 'fna_bytes', 'gff3_str', 'orf_count']
    for key in keys_to_remove:
        if key in st.session_state:
            del st.session_state[key]

uploaded_file = st.file_uploader("Upload FASTA file", on_change=reset_results)

if uploaded_file:
    input_basename = os.path.splitext(uploaded_file.name)[0]
else:
    input_basename = "out"

gff3_name = st.text_input("GFF3 output file name", value=f"{input_basename}.gff3")
faa_name = st.text_input("FAA output file name", value=f"{input_basename}.faa")
fna_name = st.text_input("FNA output file name", value=f"{input_basename}.fna")
selected_codon_table = st.selectbox("Codon Table", codon_table_choices, index=0)
trans_table = codon_table_ids[codon_table_choices.index(selected_codon_table)]
min_len = st.number_input("Minimum ORF length (aa)", min_value=1, value=50)
keep_nested = st.checkbox("Retain nested ORFs", value=False)

if st.button("Run") and uploaded_file:
    with st.spinner("Predicting ORFs..."):
        fasta_text = uploaded_file.getvalue().decode("utf-8")
        if not fasta_text.lstrip().startswith(">"):
            st.error("❌ The uploaded file does not appear to be a valid FASTA file.")
        else:
            gff3_bytes, faa_bytes, fna_bytes, gff3_str, orf_count = predict_orfs(
                fasta_text=fasta_text,
                trans_table=trans_table,
                min_protein_length=min_len,
                keep_nested=keep_nested
            )
    
            st.session_state.gff3_bytes = gff3_bytes
            st.session_state.faa_bytes = faa_bytes
            st.session_state.fna_bytes = fna_bytes
            st.session_state.gff3_str = gff3_str
            st.session_state.orf_count = orf_count
    
            st.success(f"Done! Found a total of {orf_count} ORFs.")

st.subheader("📄 Output files")
if 'gff3_bytes' in st.session_state:
    st.text_area("GFF3 Output", st.session_state.gff3_str, height=300)
    st.download_button("Download GFF3", st.session_state.gff3_bytes,
                       file_name=gff3_name, mime="text/plain")
if 'faa_bytes' in st.session_state:
    st.text_area("FAA Output", st.session_state.faa_bytes.getvalue().decode("utf-8"), height=300)
    st.download_button("Download proteins (FAA)", st.session_state.faa_bytes,
                file_name=faa_name, mime="text/plain")
if 'fna_bytes' in st.session_state:
    st.text_area("FNA Output", st.session_state.fna_bytes.getvalue().decode("utf-8"), height=300)
    st.download_button("Download CDSs (FNA)", st.session_state.fna_bytes,
                file_name=fna_name, mime="text/plain")

st.markdown("---")
st.markdown(
    "Author: [Satoshi Kawato](https://github.com/satoshikawato)  |  "
    "Source: [bio_small_scripts](https://github.com/satoshikawato/bio_small_scripts)",
    unsafe_allow_html=True
)
