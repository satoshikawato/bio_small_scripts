#!/usr/bin/env python
# coding: utf-8

import argparse
import logging
import concurrent.futures
import os
import sys
import dnaio
from dnaio.singleend import _detect_format_from_name, _detect_format_from_content
from xopen import xopen
from typing import Tuple
import io
import numpy as np
from statistics import median
# from collections import Counter

logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def parse_arguments(raw_args=None):
    parser = argparse.ArgumentParser(description="Extract and filter FASTA/FASTQ entries by various criteria")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA/FASTQ file (optionally gzipped)")
    parser.add_argument("-f", "--id_file", help="File containing list of IDs to include")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA/FASTQ file (gzipped if ending with .gz)")
    parser.add_argument("-m", "--min_length", type=int, help="Minimum sequence length to include")
    parser.add_argument("-M", "--max_length", type=int, help="Maximum sequence length to include")
    parser.add_argument("--gc_min", type=float, help="Minimum GC content to include (as a percentage)")
    parser.add_argument("--gc_max", type=float, help="Maximum GC content to include (as a percentage)")
    parser.add_argument("-q", "--qual_min", type=int, help="Minimum average quality score to include")
    parser.add_argument("-Q", "--qual_max", type=int, help="Maximum average quality score to include")
    parser.add_argument("-c", "--compression_level", type=int, default=3, help="Gzip compression level (1-9)")
    parser.add_argument("-t", "--num_threads", type=int, default=1, help="Number of threads to use for parallel processing")
    parser.add_argument("--reverse_complement", action="store_true", help="Reverse complement the sequence")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(raw_args)


complement = {ord('A'): 'T', ord('T'): 'A', ord('C'): 'G', ord('G'): 'C', ord('N'): 'N', ord('R'): 'Y', ord('Y'): 'R', ord('S'): 'S', ord('W'): 'W', ord('K'): 'M', ord('M'): 'K', ord('B'): 'V', ord('V'): 'B', ord('D'): 'H', ord('H'): 'D'}


def get_reverse_complement(sequence: str, qualities: str) -> Tuple[str, str]:
    out_sequence = sequence.translate(complement)[::-1]  # Reverse the sequence
    out_quality = qualities[::-1]  # Reverse the quality scores
    return out_sequence, out_quality


def GC(sequence: str) -> float:

    gc_count = sequence.count('G') + sequence.count('C')
    gc_percent = 100 * (gc_count / len(sequence))
    return round(gc_percent, 2)


def average_quality(qualities) -> float:
    quality_scores = np.frombuffer(qualities.encode(), dtype=np.uint8) - 33
    return np.mean(quality_scores)

def calculate_stats(records):
    read_lengths = [len(record.sequence) for record in records]
    gc_contents = [GC(record.sequence) for record in records]
    
    stats = {
        'number_of_reads': len(records),
        'average_read_length': round(sum(read_lengths) / len(read_lengths), 2),
        'median_read_length': median(read_lengths),
        'read_lengths': read_lengths,
        'average_gc_content': round(sum(gc_contents) / len(gc_contents), 2),
        'gc_contents': gc_contents
    }
    
    return stats


def calculate_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target_length = total_length * 0.5
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= target_length:
            return length
    return 0

def filter_sequence(record, ids, min_length, max_length, gc_min, gc_max, qual_min, qual_max, reverse_complement):
    seq = record.sequence
    qual = record.qualities

    if ids and record.name.split()[0] not in ids:
        return None
    if min_length is not None or max_length is not None:
        seq_length = len(seq)
        if min_length is not None and seq_length < min_length:
            return None
        if max_length is not None and seq_length > max_length:
            return None

    if gc_min is not None or gc_max is not None:
        gc = GC(seq)
        if gc_min is not None and gc < gc_min:
            return None
        if gc_max is not None and gc > gc_max:
            return None

    if qual_min is not None or qual_max is not None:
        avg_qual = average_quality(qual)
        if qual_min is not None and avg_qual < qual_min:
            return None
        if qual_max is not None and avg_qual > qual_max:
            return None

    if reverse_complement:
        record.sequence, record.qualities = get_reverse_complement(seq, qual)

    return record


def process_chunk(input_chunk, ids, min_length, max_length, gc_min, gc_max, qual_min, qual_max, reverse_complement):
    with io.BytesIO(input_chunk) as temp_file:
        records = list(dnaio.open(temp_file))
        output_chunk = [filtered_record for record in records if (filtered_record := filter_sequence(record, ids, min_length, max_length, gc_min, gc_max, qual_min, qual_max, reverse_complement)) is not None]
    return output_chunk

def main(raw_args=None):
    args = parse_arguments(raw_args)
    logging.debug(f"Filtering {args.input} and writing to {args.output}")
    id_set = set()
    num_threads = args.num_threads
    if num_threads == 0:
        num_threads = os.cpu_count()
        logging.debug(f"Using {num_threads} threads")
    if args.id_file:
        if not os.path.isfile(args.id_file) or os.path.getsize(args.id_file) == 0:
            logging.error(f"ID file {args.id_file} is empty or does not exist")
            sys.exit(1)
        else:
            with open(args.id_file) as f:
             for line in f:
                    id_set.add(line.strip())
            logging.debug(f"Loaded {len(id_set)} IDs from {args.id_file}")

    if num_threads == 1:
        dnaio_open_threads = 0
    else:
        dnaio_open_threads = num_threads
    output_format = _detect_format_from_name(args.output)
    fasta_extensions = [".fasta", ".fa", ".fas", ".fsa", ".faa", ".fna", ".ffn", ".frn", ".mpfa"]
    num_workers = num_threads + 5
    
    # Calculate input file statistics
    with dnaio.open(args.input) as input_file:
        input_records = list(input_file)
    input_stats = calculate_stats(input_records)
    input_stats['read_n50'] = calculate_n50(input_stats['read_lengths'])
    
    logging.info("Input file statistics:")
    logging.info(f"Number of reads: {input_stats['number_of_reads']}")
    logging.info(f"Average read length: {input_stats['average_read_length']}")
    logging.info(f"Median read length: {input_stats['median_read_length']}")
    logging.info(f"Read N50: {input_stats['read_n50']}")
    logging.info(f"Average GC content: {input_stats['average_gc_content']}%")
    
    with xopen(args.input, mode="rb", threads=num_threads) as input_file:
        number_of_saved_reads = 0
        filtered_records = []
        if output_format is None:
            if any(args.output.endswith(ext) for ext in fasta_extensions):
                output_format = "fasta"
            else:
                output_format = _detect_format_from_content(input_file)
        with dnaio.open(args.output, mode="w", compression_level=args.compression_level, open_threads=dnaio_open_threads, fileformat=output_format) as output_handle:
            chunk_count = 0
            for chunk in dnaio.read_chunks(input_file):
                chunk_count += 1
                chunk_in_bytes = chunk.tobytes()
                futures = []
                with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
                    future = executor.submit(process_chunk, chunk_in_bytes, id_set, args.min_length, args.max_length, args.gc_min, args.gc_max, args.qual_min, args.qual_max, args.reverse_complement)
                    futures.append(future)
                for future in concurrent.futures.as_completed(futures):
                    chunk_filtered_records = future.result()
                    number_of_saved_reads += len(chunk_filtered_records)
                    filtered_records.extend(chunk_filtered_records)
                    for record in chunk_filtered_records:
                        output_handle.write(record)
        
        output_stats = calculate_stats(filtered_records)
        output_stats['read_n50'] = calculate_n50(output_stats['read_lengths'])
        
        logging.info(f"A total of {number_of_saved_reads} reads were saved to {args.output}")
        logging.info("Filtered reads statistics:")
        logging.info(f"Number of reads: {output_stats['number_of_reads']}")
        logging.info(f"Average read length: {output_stats['average_read_length']}")
        logging.info(f"Median read length: {output_stats['median_read_length']}")
        logging.info(f"Read N50: {output_stats['read_n50']}")
        logging.info(f"Average GC content: {output_stats['average_gc_content']}%")

if __name__ == "__main__":
    main()
