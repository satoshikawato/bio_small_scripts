#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from collections import defaultdict
import csv
from typing import List, Dict, Tuple


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Add functional annotation to the 9th column of a GFF3 file')
    parser.add_argument(
        '-g', '--gff', help='Annotation in GFF3 format (required)', type=str, required=True)
    parser.add_argument(
        '-b', '--blast', help='Tab-separated BLAST output (required) with "-outfmt "[6|7] qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs""', type=str, required=True)
    parser.add_argument(
        '-f', '--func', help='Tab-separated function table (required)', type=str, required=True)
    parser.add_argument(
        '-o', '--out', help='Annotation in GFF3 format (default: stdout)', type=str)
    parser.add_argument(
        '-p', '--prefix', help='Locus ID prefix (default: gene)', type=str, default='gene')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def read_gff(gff_file: str, min_cds_len: int = 75, min_confidence: int = 70) -> Tuple[str, List[List[str]]]:
    '''
    Read a GFF3 file.
    Features shorter than min_cds_len or confidence lower than min_confidence are discarded.
    Returns comment lines and feature lines.
    '''
    comments = ''
    features = []
    with open(gff_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                comments += f'{line}\n'
            else:
                feature = line.split('\t')
                attributes = feature[8].split(';')
                attributes = list(filter(None, attributes))
                tags = dict([tmp.split('=') for tmp in attributes])
                feat_len = int(feature[4]) - int(feature[3]) + 1
                feature[8] = tags
                features.append(feature)
    return comments, features

def tbl_to_dict(function_table: str) -> Dict[str, Tuple[str, str]]:
    '''
    Read a two-column or three-column, tab-delimited table.
    Returns a dictionary (key: left; value: tuple of middle and right columns (if available)).
    '''
    func_dict = {}
    with open(function_table) as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) == 2:
                func_dict[row[0]] = (row[1], "")
            elif len(row) == 3:
                func_dict[row[0]] = (row[1], row[2])
    return func_dict



def blast_to_dict(blast_output: str, evalue_cutoff: float = 1e-2, min_identity: float = 80) -> Dict[str, List[List[str]]]:
    '''
    Read a tab-delimited BLAST output (-outfmt 6 or 7).
    Comment lines are removed automatically.
    Returns a dictionary (key: query; value: list of BLAST hits).
    '''
    blast_dict = defaultdict(list)
    with open(blast_output) as f:
        for row in csv.reader(f, delimiter='\t'):
            if row[0].startswith('#'):
                continue
            if float(row[2]) < min_identity or float(row[10]) > evalue_cutoff:
                continue
            blast_dict[row[0]].append(row)
    return blast_dict


def select_best_hits(hits: List[List[str]]) -> Tuple[str, List[str]]:
    '''
    Select the best BLAST hits based on query coverage and subject coverage.
    Returns the CDS name and a list of notes.
    '''
    hits_over_80 = {}
    hits_under_80 = []
    for hit in hits:
        aln_len, q_start, q_end, q_len, s_len = map(int, [hit[3], hit[6], hit[7], hit[12], hit[13]])
        q_cov = 100 * (aln_len / q_len)
        s_cov = 100 * (aln_len / s_len)
        s_cov_ratio = 100 * (q_len / s_len)
        if q_cov > 80:
            hits_over_80[hit[1]] = {'aln_len': aln_len, 'q_start': q_start, 'q_end': q_end,
                                    'q_len': q_len, 's_len': s_len, 'evalue': float(hit[10]),
                                    'q_cov': q_cov, 's_cov': s_cov, 's_cov_ratio': s_cov_ratio}
        else:
            hits_under_80.append(hit[1])

    notes = []
    if hits_over_80:
        cds_kept = max(hits_over_80, key=lambda x: hits_over_80[x]['s_cov'])
        cds_name = cds_kept
        if hits_over_80[cds_kept]['s_cov'] < 80:
            notes.append(f"{cds_name}%3Btruncated")
        else:
            notes.append(cds_name)
    else:
        if len(hits_under_80) > 1:
            cds_name = "-".join(hits_under_80)
            notes.append(f"{cds_name}%3Bchimeric")
        elif len(hits_under_80) == 1:
            cds_name = "-".join(hits_under_80)
            notes.append(cds_name)
        else:
            cds_name = ""
            notes.append("hypothetical protein")
    return cds_name, notes


def add_func_to_cds(features: List[List[str]], blast_dict: Dict[str, List[List[str]]],
                    func_tbl: Dict[str, Tuple[str, str]], prefix: str) -> str:
    '''
    Add annotation to the 9th column of the GFF3 file.
    Discards CDS if no significant hits to target sequences.
    Returns the updated GFF3 features as a string.
    '''
    out_features = ''
    new_count = 0
    for feature in features:
        if 'CDS' in feature[2]:
            qualifiers = feature[8]
            if qualifiers["ID"] in blast_dict:
                cds_name, notes = select_best_hits(blast_dict[qualifiers["ID"]])
                if cds_name:
                    new_count += 1
                    new_id = f"{prefix}_{str(new_count).zfill(3)}"
                    product, note = func_tbl.get(cds_name, ("hypothetical protein", ""))
                    if note:
                        notes.append(note.replace(";", "%3B"))
                    note_strs = "%3B".join(f"{note}" for note in notes)
                    out_qualifiers = f"ID={new_id};note={note_strs};product={product}"
                    out_feature_line = "\t".join(feature[:8] + [out_qualifiers])
                    out_features += f'{out_feature_line};\n'
    return out_features


def main() -> None:
    '''
    Main function to run the script.
    '''
    args = get_args()
    comments, features = read_gff(args.gff)
    func_tbl = tbl_to_dict(args.func)
    blast_dict = blast_to_dict(args.blast)
    prefix = args.prefix
    out_file = args.out
    out_features = add_func_to_cds(features, blast_dict, func_tbl, prefix)
    if out_file:
        with open(out_file, "w") as f:
            f.write(comments)
            f.write(out_features)
    else:
        print(comments, end='')
        print(out_features, end='')


if __name__ == "__main__":
    main()
