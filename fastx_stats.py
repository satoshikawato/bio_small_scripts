#!/usr/bin/env python
# coding: utf-8

import argparse
import logging
import sys
import dnaio
import os
import glob


logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
logging.basicConfig(level=logging.INFO, format='%(message)s')


def parse_arguments(raw_args=None):
    parser = argparse.ArgumentParser(description="Calculate statistics for FASTA/FASTQ files")
    parser.add_argument("-i", "--input", nargs="+", help="Input FASTA/FASTQ file(s) (optionally gzipped)")
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")
    args = parser.parse_args(raw_args)

    if not args.input:
        parser.print_help()
        sys.exit(1)

    return args


def GC(sequence: str) -> float:
    gc_count = sequence.count('G') + sequence.count('C')
    gc_percent = 100 * (gc_count / len(sequence))
    return round(gc_percent, 2)


def calculate_stats(records):
    read_lengths = [len(record.sequence) for record in records]
    gc_contents = [GC(record.sequence) for record in records]
    
    stats = {
        'num_seqs': len(records),
        'sum_len': sum(read_lengths),
        'min_len': min(read_lengths),
        'avg_len': round(sum(read_lengths) / len(read_lengths), 2),
        'max_len': max(read_lengths),
        'average_gc_content': round(sum(gc_contents) / len(gc_contents), 2),
    }
    
    return stats



def process_file(file_path):
    try:
        if file_path == "-":
            input_file = dnaio.open(sys.stdin.buffer)
        else:
            input_file = dnaio.open(file_path)

        records = list(input_file)
        stats = calculate_stats(records)
        
        file_format = "FASTQ" if records[0].qualities else "FASTA"
        file_type = "DNA"  # Assuming DNA sequences
        
        return [file_path, file_format, file_type, stats['num_seqs'], stats['sum_len'], stats['min_len'], stats['avg_len'], stats['max_len']]
    except Exception as e:
        logger.warning(f"Skipping file '{file_path}' due to an error: {str(e)}")
        return None


def main(raw_args=None):
    args = parse_arguments(raw_args)
    file_paths = []
    for file_path in args.input:
        file_path_expanded = os.path.expanduser(file_path)
        globbed_file_paths = glob.glob(file_path_expanded)
        file_paths = file_paths + globbed_file_paths

    data = [process_file(file_path) for file_path in file_paths]
    
    headers = ["file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len"]
    
    # Print the header
    header_line = "\t".join(headers)
    output_lines = [header_line]
    
    # Print the data rows
    for row in data:
        if row is not None:
            row_line = "\t".join(str(item) for item in row)
            output_lines.append(row_line)

    output = "\n".join(output_lines)

    if args.output:
        if args.verbose:
            for i in range(1, len(output_lines)):
                printed_line = output_lines[i].split("\t")
                print("Input file: ", printed_line[0])
                print("Format: ", printed_line[1])
                print("Seuqnce type: ", printed_line[2])
                print("Number of records: ", printed_line[3])
                print("Total length (bp): ", printed_line[4])
                print("Minimum length (bp): ", printed_line[5])
                print("Average length (bp): ", printed_line[6])
                print("Maximum length (bp): ", printed_line[7])
        with open(args.output, "w") as output_file:
            output_file.write(output)
    else:
        if args.verbose:
            for i in range(1, len(output_lines)):
                printed_line = output_lines[i].split("\t")
                print("Input file: ", printed_line[0])
                print("Format: ", printed_line[1])
                print("Seuqnce type: ", printed_line[2])
                print("Number of records: ", printed_line[3])
                print("Total length (bp): ", printed_line[4])
                print("Minimum length (bp): ", printed_line[5])
                print("Average length (bp): ", printed_line[6])
                print("Maximum length (bp): ", printed_line[7])
        else:
            print(output)


if __name__ == "__main__":
    main()