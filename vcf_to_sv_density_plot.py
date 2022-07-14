#!/usr/bin/env python
# coding: utf-8

import sys
import os
import re
import gzip
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _get_args():
    parser = argparse.ArgumentParser(
        description='plot the occurences of structural variations (SV) over sliding windows')
    parser.add_argument(
        '-i',
        '--input',
        help='VCF (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-w',
        '--window',
        help='window size (default: 1000)',
        type=int,
        default=1000)
    parser.add_argument(
        '-s',
        '--step',
        help='step size (default: 100)',
        type=int,
        default=100)
    parser.add_argument(
        '-m',
        '--max_sv_len',
        help='maximum SV length (default: 10)',
        type=int,
        default=10)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args


def get_header(in_file):
    with open(in_file, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                header = [x for x in line.replace('\n', '').split('\t')]
                break
    return header


def get_contig_len(in_file):
    with open(in_file, "rt") as ifile:
        contig_dict = {}
        for line in ifile:
            if line.startswith("##contig="):
                s = line
                m = re.search(r'<.+>', s)
                contig = m.group(0)
                contig_info = line.split(',')
                contig_id = contig_info[0].split('=')[2]
                contig_length = contig_info[1].split('=')[1].replace('>\n', '')
                contig_dict[contig_id] = contig_length
    return contig_dict


def get_strain_names(header):
    strains = header[9:]
    strains_out = [i.split('.')[0] for i in strains]
    return strains_out


def filter_sv(df, max_sv_len):
    filtered_sv = []
    for index, row in df.iterrows():
        len_list = []
        ref_seq = row["REF"]
        len_list.append(len(ref_seq))
        alt_seq = [x for x in row["ALT"].split(',')]
        len_list.append(len(max(alt_seq)))
        if max(len_list) <= max_sv_len:
            filtered_sv.append(row)
    df_sv = pd.DataFrame(filtered_sv)
    return df_sv


def get_snp_density(df_snp, length_of_sequence, window_size, interval):
    df_snp2 = df_snp["POS"].to_frame()
    df_snp2.insert(1, "val", 1, True)
    df_tmp = pd.DataFrame({"POS": range(1, length_of_sequence + 1)})
    df_merge = pd.merge(df_tmp, df_snp2, how="outer")
    df_merge.fillna(0, inplace=True)
    out_dict = {}
    for start in range(1, length_of_sequence + 1, interval):
        end = start + window_size
        if end > length_of_sequence:
            break
        out_dict[start] = df_merge["val"].iloc[start:end].sum(axis=0)
    snp_density_df = pd.DataFrame(
        list(
            out_dict.items()),
        columns=[
            'POS',
            'SV_COUNT'])
    return snp_density_df


def generate_plot(snp_density_df):
    fig = plt.figure(figsize=(25, 1))
    ax = fig.add_subplot(1, 1, 1)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    x = snp_density_df["POS"]
    y = snp_density_df["SV_COUNT"]
    ax.set_xlim(0, max(x))
    ax.fill_between(x, y, y2=0, color="#9cd2f7")
    plt.rcParams['svg.fonttype'] = 'none'
    return fig


def main():
    args = _get_args()
    in_vcf = args.input
    window = args.window
    step = args.step
    max_sv_len = args.max_sv_len
    header_line = get_header(in_vcf)
    contig_dict = get_contig_len(in_vcf)
    df = pd.read_csv(
        in_vcf,
        comment='#',
        sep='\t',
        header=None,
        names=header_line)
    df_snp = filter_sv(df, max_sv_len)
    for contig_id in contig_dict.keys():
        contig_length = int(contig_dict[contig_id])
        df_snp_density = get_snp_density(df_snp, contig_length, window, step)
        figure = generate_plot(df_snp_density)
        fig_name = '{}.svg'.format(contig_id)
        figure.savefig(fig_name, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
