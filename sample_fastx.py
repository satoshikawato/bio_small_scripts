#!/usr/bin/env python
# coding: utf-8

import argparse
import logging
import os
import sys
import dnaio
from dnaio.singleend import _detect_format_from_name, _detect_format_from_content
from xopen import xopen
import random
import io
import concurrent.futures

logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract and filter FASTA/FASTQ entries by various criteria")
    parser.add_argument("-i", "--input", required=True, type=str, help="Input FASTA/FASTQ file (optionally gzipped)")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output FASTA/FASTQ file (gzipped if ending with .gz)")
    parser.add_argument("-t", "--num_threads", type=int, default=0, help="Number of threads to use (default: use all available cores)")
    parser.add_argument("-c", "--compression_level", type=int, default=6, help="Compression level for gzip (default: 6)")
    
    # the following arguments are mutually exclusive. To do so, we can use the add_mutually_exclusive_group() method
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "--fraction", type=float, help="Fraction of reads to keep")
    group.add_argument("-n", "--num_reads", type=int, help="Number of reads to keep")
    return parser.parse_args()

def get_ids(input_file):
    ids = set()
    for record in dnaio.open(input_file):
        ids.add(record.name.split()[0])
    return ids


#
def randomly_select_ids(ids, num_reads):
    # return random.sample(ids, num_reads) # DeprecationWarning: Sampling from a set deprecated since Python 3.9 and will be removed in a subsequent version.
    # if the number of reads to be selected is greater than the number of reads in the file, then return all the reads
    if num_reads > len(ids):
        return ids
    else:
        return random.sample(list(ids), num_reads)


def determine_read_number_from_fraction(ids, fraction):
    return int(len(ids) * fraction)


def is_selected(record, selected_ids):
    if record.name.split()[0] in selected_ids:
        return record
    else:
        return None

def process_chunk(input_chunk, selected_ids):
    output_chunk = []
    with io.BytesIO(input_chunk) as temp_file:
        records = list(dnaio.open(temp_file))
        output_chunk = [record for record in records if is_selected(record, selected_ids)]

    return output_chunk

def main():
    args = parse_arguments()
    logging.info(f"Filtering {args.input} and writing to {args.output}")
    id_set = set()
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
    num_workers = num_threads + 5
    if output_format is None:
        if args.output.endswith(".gz"):
            with xopen(args.output, "rb") as f:
                output_format = _detect_format_from_content(f.read(1000))
        else:
            if any(args.output.endswith(ext) for ext in fasta_extensions):
                output_format = "fasta"
            else:
                output_format = "fastq"
    if output_format is None:
        raise ValueError(f"Could not determine output format from file name {args.output}")
    id_set = get_ids(args.input)

    logging.info(f"Read {len(id_set)} reads from {args.input}")


    if args.fraction:
        number_of_saved_reads = determine_read_number_from_fraction(id_set, args.fraction)
        selected_ids = randomly_select_ids(id_set, number_of_saved_reads)
    elif args.num_reads:
        number_of_saved_reads = args.num_reads
        selected_ids = randomly_select_ids(id_set, number_of_saved_reads)
    else:
        raise ValueError("Either --fraction or --num_reads must be given")
    number_of_saved_reads = 0
    with xopen(args.input, mode="rb", threads=num_threads) as input_file:
        with dnaio.open(args.output, mode="w", compression_level=args.compression_level, open_threads=dnaio_open_threads, fileformat=output_format) as output_handle:
            chunk_count = 0
            for chunk in dnaio.read_chunks(input_file):  # dnaio.read_chunks() is a generator that yields chunks of sequences from the input file. What happens when we use concurrent.futures.ThreadPoolExecutor within the iterator? The answer is that the iterator will yield chunks of sequences from the input file. The ThreadPoolExecutor will process the chunks of sequences in parallel.
                chunk_count += 1
                chunk_in_bytes = chunk.tobytes()
                futures = []
                with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
                    future = executor.submit(process_chunk, chunk_in_bytes, selected_ids)
                    futures.append(future)
                for future in concurrent.futures.as_completed(futures):
                    number_of_saved_reads += len(future.result())
                    for record in future.result():
                        output_handle.write(record)

                    
    logging.info(f"Done. Saved {number_of_saved_reads} reads to {args.output}")
if __name__ == "__main__":
    main()
