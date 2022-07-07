#!/usr/bin/env python
# coding: utf-8


import sys
import os
import argparse
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def _get_args():
    parser = argparse.ArgumentParser(
        description='Generate dinucleotide skew plot(s) of FASTA format DNA sequences in SVG format. Plots are saved separately for each entry in a multifasta file')
    parser.add_argument(
        '-i',
        '--input',
        help='Fasta (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-n',
        '--nt',
        help='dinucleotide (default: GC). ',
        type=str,
        default="GC")
    parser.add_argument(
        '-w',
        '--window',
        help='window size (default: 1000) ',
        type=int,
        default="1000")
    parser.add_argument(
        '-s',
        '--step',
        help='step size (default: 100) ',
        type=int,
        default="100")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    #args = parser.parse_args(args=['--input','input.fa'])
    return args


def fasta_to_records(in_fa):
    seq_records = [record for record in SeqIO.parse(in_fa, 'fasta')]
    return seq_records


def calculate_dinucleotide_skew(seq, base1, base2):
    "Return dinucleotide skew in a given sequence"
    base1_count = seq.count(base1)
    base2_count = seq.count(base2)
    skew = ((base1_count - base2_count) / (base1_count + base2_count))
    return skew


def sliding_window(seq, window, step):
    for start in range(0, len(seq), step):
        end = start + window
        if end > len(seq):
            break
        out_seq = seq[start:end]
        yield start, out_seq


def skew_df(record, window, step, nt):
    nt_list = list(nt)
    nt_1 = nt_list[0]
    nt_2 = nt_list[1]
    skew_sum = 0
    skew_dict = {}
    skew_cumulative_dict = {}
    for start, seq_part in sliding_window(record.seq, window, step):
        skew = calculate_dinucleotide_skew(seq_part, nt_1, nt_2)
        skew_dict[start] = skew
        skew_sum = (skew_sum + skew)
        skew_cumulative_dict[start] = (skew_sum)
    max_skew_abs = abs(max(skew_dict.values(), key=abs))
    max_skew_cumulative_abs = abs(max(skew_cumulative_dict.values(), key=abs))
    factor = (max_skew_abs / max_skew_cumulative_abs)
    skew_cumulative_dict.update((x, (y * factor))
                                for x, y in skew_cumulative_dict.items())
    skew_legend = "{} skew".format(nt)
    cumulative_skew_legend = "Cumulative {} skew, normalized".format(nt)
    df = pd.DataFrame({skew_legend: pd.Series(skew_dict),
                       cumulative_skew_legend: pd.Series(skew_cumulative_dict)})
    return df


def draw_plot(record, df):
    fig = plt.figure(figsize=(30, 5))
    ax = fig.add_subplot(111)
    record_id = record.id
    ax.set_title(record_id, fontsize=20)
    ax.set_xlim([0, len(record.seq)])
    ax.ticklabel_format(style='plain', axis='both')
    ax.tick_params(axis="both", labelsize=20)
    coordinate_of_max_cumulative_skew = df.iloc[:, 1].idxmax()
    max_cumulative_skew = df.iloc[:, 1].max()
    fig.tight_layout()
    ax.plot(df)
    ax.axvline(x=coordinate_of_max_cumulative_skew, color="red")
    ax.text(
        x=coordinate_of_max_cumulative_skew,
        y=max_cumulative_skew * 0.8,
        s=coordinate_of_max_cumulative_skew,
        fontsize=20)
    ax.legend(df.columns, fontsize=20)
    filename = "{}.svg".format(record_id)
    fig.savefig(filename)


def main():
    args = _get_args()
    in_fa = args.input
    window = args.window
    step = args.step
    nt = args.nt
    if len(nt) != 2:
        nt = "GC"
    nt_list = list(nt)
    window = 1000
    step = 100
    records = fasta_to_records(in_fa)
    for record in records:
        df = skew_df(record, window, step, nt)
        draw_plot(record, df)


if __name__ == "__main__":
    main()
