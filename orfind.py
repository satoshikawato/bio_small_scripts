#!/usr/bin/env python
# coding: utf-8

# Basically rewrite of this:
# Biopython Tutorial and Cookbook 20.1.13 Identifying open reading frames
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec384

import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


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
    non_nested_orfs = []
    for i in range(len(orf_list)):
        is_nested = False
        for j in range(len(orf_list)):
            if i != j:
                if orf_list[i][0] >= orf_list[j][0] and orf_list[i][1] <= orf_list[j][1]:
                    is_nested = True
                    break
        if not is_nested:
            non_nested_orfs.append(orf_list[i])
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
    orf_list.sort()
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
    annot_lines = ''
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
            annot_lines += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                seqid, source, feature_type, start, end, score, strand, phase, qualifers)
    return annot_lines


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


if __name__ == "__main__":
    main()
