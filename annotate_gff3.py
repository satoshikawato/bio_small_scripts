#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
from collections import defaultdict


def get_args():
    parser = argparse.ArgumentParser(
        description='add functional annotation to the 9th column of a gff3 file')
    parser.add_argument(
        '-g',
        '--gff',
        help='Annotation in gff3 format (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-b',
        '--blast',
        help='tab-separated blast output (required) with "-outfmt "[6|7] qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs""',
        type=str,
        required=True)
    parser.add_argument(
        '-f',
        '--func',
        help='tab-separated function table (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-o',
        '--out',
        help='Annotation in gff3 format (default: stdout)',
        type=str)
    parser.add_argument(
        '-p',
        '--prefix',
        help='locus ID prefix (default: gene)',
        type=str,
        default='gene')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args


def read_gff(in_gff):
    '''
    Read a gff3 file.
    Features shorter than 75 nt or confidence lower than 70 are discarded.
    Returns: comment lines and feature lines.

    Also see:
    GFF3 File Format - Definition and supported options (http://gmod.org/wiki/GFF3)
    '''
    comments = ''
    features = []
    min_cds_len = 75
    with open(in_gff) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                comments += '{}\n'.format(line)
            else:
                feature = line.split('\t')
                attributes = feature[8].split(';')
                attributes = list(filter(None, attributes))
                tags = dict([tmp.split('=') for tmp in attributes])
                feat_len = (int(feature[4]) - int(feature[3]) + 1)
                feature[8] = tags
                features.append(feature)
    return comments, features


def tbl_to_dict(in_table):
    '''
    Read a two-column, tab-delimited table.
    Returns a dictionary (key: left; value: right).
    '''
    out_dict = {}
    with open(in_table) as f:
        for line in f:
            line = line.rstrip().split('\t')
            out_dict[line[0]] = line[1]
    return out_dict


def blast_to_dict(blast_out):
    '''
    Read a tab-delimited blast output (-outfmt 6 or 7).
    Comment lines removed automatically.
    Returns a dictionary (key: query; value: whole line, tab-delimited).
    '''
    blast_out_dict = defaultdict(list)
    with open(blast_out) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            else:
                line = line.split('\t')
                cds_id = line[0]
                if (float(line[2]) < 80) or (
                        float(line[10]) > 1e-2) or (float(line[11]) < 80):
                    continue
                else:
                    blast_out_dict[cds_id].append(line)
    return blast_out_dict


def judge_overlaps(retained_hits):
    '''

    1)  The subject name is assigned that covers 80% or more of the query.
        Multiple hits are ranked based on s_cov (subject coverage), which favors the CDS of the closest length over longer/shorter ones.
        If the query length is less than 80% of the selected subject length, the query is regarded as truncated.

    2)  If none of the subjects cover 80% of the query, multiple hits are concatenated with hyphen (-) to denote a fusion gene.

    Returns a CDS name.

    '''
    hits2 = {}
    hits2_over80 = {}
    hits2_under80 = []
    cds_note_list = []
    for hit in retained_hits:
        hit_dict = {}
        aln_len = int(hit[3])
        q_start = int(hit[6])
        q_end = int(hit[7])
        q_len = int(hit[12])
        s_len = int(hit[13])
        hit_dict["aln_len"] = aln_len
        hit_dict["q_start"] = q_start
        hit_dict["q_end"] = q_end
        hit_dict["q_len"] = q_len
        hit_dict["s_len"] = s_len
        hit_dict["evalue"] = float(hit[10])
        hit_dict["q_cov"] = float(100 * (aln_len / q_len))
        hit_dict["s_cov"] = float(100 * (aln_len / s_len))
        hit_dict["s_cov_ratio"] = float(100 * (q_len / s_len))
        hits2[hit[1]] = hit_dict

    hits2_under80 = []

    for key in hits2.keys():
        if hits2[key]['q_cov'] > 80:
            hits2_over80[key] = hits2[key]
        else:
            hits2_under80.append(key)

    if len(hits2_over80) > 0:
        cds_kept = max(hits2_over80, key=lambda x: hits2_over80[x]['s_cov'])
        cds_name = "{}".format(cds_kept)
        if hits2_over80[cds_kept]['s_cov'] < 80:
            cds_note_list.append("{}, truncated".format(cds_name))
        else:
            cds_note_list.append("{}".format(cds_name))
    else:
        if len(hits2_under80) > 1:
            cds_name = "-".join(hits2_under80)
            cds_note_list.append("{}, chimeric".format(cds_name))
        elif len(hits2_under80) == 1:
            cds_name = "-".join(hits2_under80)
            cds_note_list.append(cds_name)
        else:
            cds_note_list.append("hypothetical protein")
    return cds_name, cds_note_list


def add_func_to_cds(features, blast_dict, func_tbl, prefix):
    '''
    Add annotation to the 9th column of the gff3 file.
    Discards CDS if no significant hits to target sequences.
    '''
    out_features = ''
    new_count = 0
    for feature in features:
        if 'CDS' in feature[2]:
            qualifiers = feature[8]
            new_qualifiers = defaultdict(list)
            if qualifiers["ID"] in blast_dict.keys():
                cds_name, cds_note_list = judge_overlaps(
                    blast_dict[qualifiers["ID"]])
                new_count += 1
                new_id = "{}_{}".format(prefix, str(new_count).zfill(3))
                new_qualifiers["ID"] = new_id
                for note in cds_note_list:
                    new_qualifiers["note"].append(note)
                if cds_name in func_tbl.keys():
                    new_qualifiers["product"] = func_tbl[cds_name]
                else:
                    new_qualifiers["product"] = "hypothetical protein"
            else:
                continue
            note_strs = ""
            for note in new_qualifiers["note"]:
                note_strs += "note={};".format(note)
            out_qualifiers = "ID={};{}product={};".format(
                new_qualifiers["ID"], note_strs, new_qualifiers["product"])
            out_feature_line = feature[0:8]
            out_feature_line.append(out_qualifiers)
            out_feature_line = '\t'.join(out_feature_line)
            out_features += '{}\n'.format(out_feature_line)
    return out_features


def main():
    args = get_args()
    comments, features = read_gff(args.gff)
    func_tbl = tbl_to_dict(args.func)
    blast_dict = blast_to_dict(args.blast)
    prefix = args.prefix
    out_file = args.out
    out_features = add_func_to_cds(features, blast_dict, func_tbl, prefix)
    if out_file:
        with open(out_file, "w") as f:
            print(comments, file=f, end='')
            print(out_features, file=f, end='')
    else:
        print(comments, end='')
        print(out_features, end='')


if __name__ == "__main__":
    main()
