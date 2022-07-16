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
    parser = argparse.ArgumentParser(description='BAM depths and read alignment breakpoints in PNG')
    parser.add_argument(
        '-i',
        '--input',
        help='BAM file (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-w',
        '--window',
        help='window size (default: 100) ',
        type=int,
        default="100")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def process_bam(bamfile, window):
    bam_name = os.path.splitext(os.path.basename(
        str(bamfile.filename.decode('UTF-8'))))[0]
    for ref_name in bamfile.references:
        ref_length = bamfile.get_reference_length(ref_name)
        df1 = pd.DataFrame(
            0,
            index=range(ref_length),
            columns=['count']).astype('int')
        for read in bamfile.fetch(ref_name):
            if read.mapping_quality >= 60:
                cord_start = read.reference_start
                cord_end = read.reference_end - 1
                df1.at[cord_start, 'count'] += 1
                df1.at[cord_end, 'count'] += 1
            else:
                continue
        df2 = df1.rolling(window).sum().fillna(0)

        df5 = pd.DataFrame(0, index=range(ref_length), columns=['depth'])
        for pileupcolumn in bamfile.pileup(ref_name, min_mapping_quality=60):
            df5.at[pileupcolumn.pos, 'depth'] = pileupcolumn.n
        df6 = df5.rolling(window).mean().fillna(0)

        fig, ax = plt.subplots(figsize=(10, 4))
        fig.tight_layout()
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_xlim(0, ref_length)
        df2_label = 'breakpoint (sum per {} bp window)'.format(window)
        df6_label = 'depth (average depth per {} bp window)'.format(window)
        df2.plot(y='count', ax=ax, label=df2_label)
        df6.plot(y='depth', ax=ax, label=df6_label)
        plot_file = "{}_{}.svg".format(ref_name, bam_name)
        fig.savefig(plot_file, format="svg", bbox_inches="tight")


def main():
    args = _get_args()
    in_bam = args.input
    window = args.window
    bamfile = pysam.AlignmentFile(in_bam, "rb")
    process_bam(bamfile, window)


if __name__ == "__main__":
    main()
