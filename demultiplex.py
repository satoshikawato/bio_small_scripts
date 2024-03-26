#!/usr/bin/env python
# coding: utf-8
import cProfile
import edlib
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
from collections import defaultdict
import numpy as np


logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

complement = {ord('A'): 'T', ord('T'): 'A', ord('C'): 'G', ord('G'): 'C', ord('N'): 'N', ord('R'): 'Y', ord('Y'): 'R', ord('S'): 'S', ord('W'): 'W', ord('K'): 'M', ord('M'): 'K', ord('B'): 'V', ord('V'): 'B', ord('D'): 'H', ord('H'): 'D'}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract and filter FASTA/FASTQ entries by various criteria")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA/FASTQ file (optionally gzipped)")
    parser.add_argument("-f", "--adapter_file", help="FASTA File containing barcodes")
    parser.add_argument("-t", "--num_threads", type=int, default=0, help="Number of threads to use (default: 0, meaning all available cores)")
    parser.add_argument("-c", "--compression_level", type=int, default=3, help="Compression level for gzip output (default: 3)")
    return parser.parse_args()

def get_reverse_complement(sequence: str, qualities: str) -> Tuple[str, str]:
    out_sequence = sequence.translate(complement)[::-1]  # Reverse the sequence
    out_quality = qualities[::-1]  # Reverse the quality scores
    return out_sequence, out_quality


def get_N(sequence: str) -> float:
    N_count = sequence.count('N')
    N_percent = 100 * (N_count / len(sequence))
    return round(N_percent, 2)

def GC(sequence: str) -> float:
    gc_count = sequence.count('G') + sequence.count('C')
    gc_percent = 100 * (gc_count / len(sequence))
    return round(gc_percent, 2)


def average_quality(qualities) -> float:
    quality_scores = np.frombuffer(qualities.encode(), dtype=np.uint8) - 33
    return np.mean(quality_scores)

def filter_sequence(record, primers):
    primer_sites = {}
    for primer in primers:
        query_f = primer.sequence
        query_r = primer.sequence.translate(complement)[::-1]
        
        forward_alignment = edlib.align(query_f, record.sequence, mode="HW", task="path", additionalEqualities=[("R", "A"), ("R", "G")])
        if forward_alignment["editDistance"] >2:
            forward_alignment = None
        reverse_alignment = edlib.align(query_r, record.sequence, mode="HW", task="path", additionalEqualities=[("R", "A"), ("R", "G")])
        if reverse_alignment["editDistance"] >2:
            reverse_alignment = None
        if forward_alignment != None or reverse_alignment != None:
            primer_sites[primer.id] = [forward_alignment, reverse_alignment]
    if len(primer_sites)>1: #  flag sequences with more than one primer as failed reads. If only a single primer is found and both forward and reverse alignments are found, then the sequence is a pass. Since we have to demultiplex reads, records should be classified according to the primer that it aligns to.
        return record, "chimeric"
    elif len(primer_sites) == 1:
        key = list(primer_sites.keys())[0]
        # if both forward and reverse alignments are found, then the sequence is a pass
        if primer_sites[key][0] != None and primer_sites[key][1] != None:
            return record, key
        else:
            return record, "failed" 
    else:
        return record, "failed"


def process_chunk(input_chunk, primers):
    output_chunk = defaultdict(list)
    with io.BytesIO(input_chunk) as temp_file:  # this function is used to read the input_chunk without creating a temporary file
        for record in dnaio.open(temp_file):
            record, destination = filter_sequence(record, primers)
            output_chunk[destination].append(record)
    return output_chunk

def main():
    args = parse_arguments()
    logging.info(f"Filtering {args.input}")
    id_set = set()
    num_threads = args.num_threads
    if num_threads == 0:
        num_threads = os.cpu_count()
        logging.info(f"Using {num_threads} threads")
    if num_threads == 1:
        dnaio_open_threads = 0
    else:
        dnaio_open_threads = num_threads
    input_basename = os.path.basename(args.input).split(".")[0]
    num_workers = num_threads
    primers = []
    with dnaio.open(args.adapter_file, mode="r", open_threads=0) as reader:  # type: ignore
        for number_of_records, record in enumerate(reader, start=1):
            primers.append(record)
    with xopen(args.input, mode="rb", threads=num_threads) as input_file:
        number_of_saved_reads = 0
        output_format = _detect_format_from_content(input_file)
        open_flag = defaultdict(bool)  # Flag to track if an output file has been opened

        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = []
            for chunk in dnaio.read_chunks(input_file, buffer_size=4194304):
                future = executor.submit(process_chunk, chunk.tobytes(), primers)
                futures.append(future)

            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                for key in result.keys():
                    if not open_flag[key]:  # If the output file hasn't been opened yet
                        open_mode = "w"  # Open in write mode
                        open_flag[key] = True  # Set the flag to True
                    else:
                        open_mode = "a"  # Open in append mode

                    with dnaio.open(f"{input_basename}_{key}.{output_format}.gz", mode=open_mode, compression_level=args.compression_level, open_threads=dnaio_open_threads, fileformat=output_format) as output_handle:
                        for record in result[key]:
                            output_handle.write(record)
                            number_of_saved_reads += 1

        logging.info(f"A total of {number_of_saved_reads} reads were saved")

if __name__ == "__main__":
    cProfile.run('main()', filename='main.prof')
