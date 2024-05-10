#!/usr/bin/env python
# coding: utf-8
import os
import sys
import argparse
import random
import pysam


def _get_args(raw_args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Subsample a fixed number of reads per contig in a BAM file. Currently supports single-end reads only")
    parser.add_argument("-i", "--input", required=True,type=str,
                        help="Input BAM file")
    parser.add_argument("-o", "--output", dest="output",
                        default=None,type=str,
                        metavar="path", help="Output file name (default: basname + '.subsampled.bam')")
    parser.add_argument("-c", "--count", dest="count", default=1000,
                        type=int, help="Number of retained reads per contig")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(raw_args)


def process_bam(input_bam: str, output_bam: str, reads_per_contig: int) -> None:

    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        # Copy the header from the input BAM file
        header = in_bam.header.copy()

        # Iterate over the contigs and select a fixed number of reads per contig. If the contig has fewer reads than the specified number, all reads are selected.
        all_selected_reads = []
        for contig in in_bam.references:
            contig_reads = [read for read in in_bam.fetch(contig)]
            selected_reads = random.sample(contig_reads, min(reads_per_contig, len(contig_reads)))
            all_selected_reads.extend(selected_reads)

        # Sort all the selected reads
        all_selected_reads.sort(key=lambda x: (x.reference_name, x.reference_start))

        # Create the output BAM file
        with pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
            # Write the sorted selected reads to the output BAM file
            for read in all_selected_reads:
                out_bam.write(read)

        # Create the .bai file for the output BAM file
        pysam.index(output_bam)

def main(raw_args=None):
    args = _get_args(raw_args)
    input = args.input
    output = args.output
    count = args.count
    input_basename = '.'.join(os.path.basename(args.input).split(".")[:-1])
    if not output:
        output = f"{input_basename}.subsampled.bam"
    process_bam(input, output, count)



if __name__ == "__main__":
    main()