#!/usr/bin/env python
# coding: utf-8

import argparse

def detect_overlap(sequence):
    half_length = len(sequence) // 2
    for overlap_length in range(half_length, 0, -1):
        if sequence[:overlap_length] == sequence[-overlap_length:]:
            return overlap_length
    return 0

def detect_self_loops(input_file, output_file, minlen):
    sequences = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith("S"):
                parts = line.split("\t")
                sequences[parts[1]] = parts[2]

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("L"):
                parts = line.split("\t")
                if parts[1] == parts[3]:  # Check if source == target
                    # Check for overlaps at the ends and trim the sequence if necessary
                    if parts[1] in sequences:
                        seq = sequences[parts[1]]
                        overlap = detect_overlap(seq)
                        trimmed_seq = seq[:-overlap] if overlap else seq
                        if len(trimmed_seq) >= minlen:
                            outfile.write(">" + parts[1] + "\n" + trimmed_seq + "\n")

def main():
    parser = argparse.ArgumentParser(description="Detect self-loops in a GFA file, trim overlaps, and output them in a FASTA format")
    parser.add_argument('-i', '--input', required=True, help='Input GFA file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    parser.add_argument('--minlen', type=int, default=0, help='Minimum sequence length to be included in the output (after trimming overlaps)')
    args = parser.parse_args()

    detect_self_loops(args.input, args.output, args.minlen)

if __name__ == "__main__":
    main()
