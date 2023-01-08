#!/usr/bin/env python
# coding: utf-8

import sys
import pysam
import argparse

def _get_args():
    parser = argparse.ArgumentParser(
        description='Extract partial reads mapping to particular region from BAM')
    parser.add_argument(
        "--input",
        "-i",
        "--in",
        metavar="FILE",
        help="Input BAM file",
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

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def get_fasta(bamfile, start, end, ref):
    ref_length = bamfile.get_reference_length(ref)
    
    if end > ref_length  :
        end = ref_length 
    for read in bamfile.fetch(ref, start , end ):
        read_tuple = read.get_aligned_pairs()
        read_tbl = []
        possible_ends = []
        read_start = 0
        read_end = 0
        for tpl in read_tuple:
            if tpl[1] is not None and tpl[0] is not None:
                possible_ends.append(tpl[0])
            if start == tpl[1]:
                if tpl[0] != None:
                    read_start = tpl[0]
                else:
                    read_start = max(possible_ends)
            if end == tpl[1]:
                read_end = tpl[0]
        if read_end == None:
            read_end = max(possible_ends)
        if read_start > read_end:
            read_start_tmp = read_end
            read_end = read_start
            read_start = read_start_tmp
        for i in range(read_start, read_end):
            if i <= len(read.seq)-1:
                read_tbl.append(read.seq[i])
            else:
                break
        read_name = "{}_{}-{}".format(read.qname,read_start,read_end)
        print(">{}".format(read_name))
        print("".join(read_tbl))

def main():
    args = _get_args()
    bamfile = pysam.AlignmentFile(args.input, "rb")
    start = args.start
    end = args.end
    ref = args.ref
    get_fasta(bamfile, start, end, ref)

if __name__ == "__main__":
    main()

