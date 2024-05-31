#!/usr/bin/env python
# coding: utf-8

import argparse
import logging
import os
import sys
import dnaio
from dnaio.singleend import _detect_format_from_name, _detect_format_from_content
from xopen import xopen

logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def parse_arguments(raw_args=None):
   parser = argparse.ArgumentParser(description="Rename FASTA sequence IDs with a given prefix and numbering")
   parser.add_argument("-i", "--input", required=True, help="Input FASTA file (optionally gzipped)")
   parser.add_argument("-o", "--output", required=True, help="Output FASTA file (gzipped if ending with .gz)")
   parser.add_argument("-p", "--prefix", required=True, help="Prefix to use for renaming")
   parser.add_argument("-n", "--start_number", type=int, default=1, help="Starting number for renaming (default: 1)")
   parser.add_argument("-c", "--compression_level", type=int, default=3, help="Gzip compression level (1-9)")
   parser.add_argument("-t", "--num_threads", type=int, default=1, help="Number of threads to use for parallel processing")
   if len(sys.argv) == 1:
       parser.print_help(sys.stderr)
       sys.exit(1)
   return parser.parse_args(raw_args)


def rename_sequence(record, prefix, start_number):
   record.name = f"{prefix}_{start_number}"
   return record


def main(raw_args=None):
   args = parse_arguments(raw_args)
   logging.info(f"Renaming FASTA IDs in {args.input} and writing to {args.output}")
   num_threads = args.num_threads
   if num_threads == 0:
       num_threads = os.cpu_count()
       logging.info(f"Using {num_threads} threads")
   if num_threads == 1:
       dnaio_open_threads = 0
   else:
       dnaio_open_threads = num_threads
   output_format = _detect_format_from_name(args.output)
   fasta_extensions = [".fasta", ".fa", ".fas", ".fsa", ".faa", ".fna", ".ffn", ".frn", ".mpfa"]
   with xopen(args.input, mode="rb", threads=num_threads) as input_file:
       number_of_renamed_sequences = 0
       if output_format is None:
           if any(args.output.endswith(ext) for ext in fasta_extensions):
               output_format = "fasta"
           else:
               output_format = _detect_format_from_content(input_file)
       with dnaio.open(args.output, mode="w", compression_level=args.compression_level, open_threads=dnaio_open_threads, fileformat=output_format) as output_handle:
           seq_number = args.start_number
           for record in dnaio.open(input_file):
               renamed_record = rename_sequence(record, args.prefix, seq_number)
               output_handle.write(renamed_record)
               seq_number += 1
               number_of_renamed_sequences += 1
       logging.info(f"A total of {number_of_renamed_sequences} sequences were renamed and saved to {args.output}")


if __name__ == "__main__":
   main()