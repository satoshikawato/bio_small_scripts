#!/usr/bin/env python
# coding: utf-8

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pysam
import argparse

def _get_args():
    parser = argparse.ArgumentParser(description='Visualize strand-specific read coverage in SVG using matplotlib')
    parser.add_argument(
        '-f',
        '--forward',
        help='BAM file (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-r',
        '--reverse',
        help='BAM file (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-w',
        '--window',
        help='window size (default: 10)',
        type=int,
        default="10")
    parser.add_argument(
        '--ref',
        help='Reference name in BAM file (required)',
        type=str,
        required=True)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def process_bam(forward_bamfile, reverse_bamfile, window, ref_name):
    bam_name = os.path.splitext(os.path.basename(str(forward_bamfile.filename.decode('UTF-8'))))[0]
    ref_length = forward_bamfile.get_reference_length(ref_name)
    df_forward = pd.DataFrame(0, index=range(ref_length), columns=['depth'])
    df_reverse = pd.DataFrame(0, index=range(ref_length), columns=['depth'])
    for pileupcolumn in forward_bamfile.pileup(ref_name, min_mapping_quality=60):
        df_forward.at[pileupcolumn.pos, 'depth'] = - pileupcolumn.n
    df_forward_average = df_forward.rolling(window).mean().fillna(0)
    for pileupcolumn in reverse_bamfile.pileup(ref_name, min_mapping_quality=60):
        df_reverse.at[pileupcolumn.pos, 'depth'] = pileupcolumn.n
    df_forward_average = df_forward.rolling(window).mean().fillna(0)
    df_reverse_average = df_reverse.rolling(window).mean().fillna(0)
    fig, ax = plt.subplots(figsize=(20, 2))
    fig.tight_layout()
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.spines.bottom.set_visible(False)
    ax.set_xlim(0, ref_length)
    df_forward_average_label = 'depth (average depth per {} bp window)'.format(window)
    ax.fill_between(df_reverse_average.index, 0, df_reverse_average['depth'], color='#1f77b4', alpha=1)
    ax.fill_between(df_forward_average.index, 0, df_forward_average['depth'], color='#ff7f0e', alpha=1)

    plot_file = "{}_{}.svg".format(ref_name, bam_name)
    fig.savefig(plot_file, format="svg", bbox_inches="tight")

def main():
    args = _get_args()
    forward_bam = args.forward
    reverse_bam = args.reverse
    window = args.window
    ref_name = args.ref
    forward_bamfile = pysam.AlignmentFile(forward_bam, "rb")
    reverse_bamfile = pysam.AlignmentFile(reverse_bam, "rb")
    process_bam(forward_bamfile, reverse_bamfile, window, ref_name)

if __name__ == "__main__":
    main()
