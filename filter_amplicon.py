#!/usr/bin/env python
# coding: utf-8
import edlib
import argparse
import logging
import concurrent.futures
import os
import sys
import dnaio
from dnaio.singleend import _detect_format_from_content
from xopen import xopen
from typing import Tuple
import io
from collections import defaultdict


logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

complement = {ord('A'): 'T', ord('T'): 'A', ord('C'): 'G', ord('G'): 'C', ord('N'): 'N', ord('R'): 'Y', ord('Y'): 'R', ord('S'): 'S', ord('W'): 'W', ord('K'): 'M', ord('M'): 'K', ord('B'): 'V', ord('V'): 'B', ord('D'): 'H', ord('H'): 'D'}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract and filter FASTQ entries")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file (optionally gzipped)")
    parser.add_argument('-f', '--fwd_seq', help='forward primer sequence', type=str, default="")
    parser.add_argument('-r', '--rev_seq', help='reverse primer sequence', type=str, default="")
    parser.add_argument("-t", "--num_threads", type=int, default=0, help="Number of threads to use (default: 0, meaning all available cores)")
    parser.add_argument("-c", "--compression_level", type=int, default=3, help="Compression level for gzip output (default: 3)")
    return parser.parse_args()

def get_reverse_complement(sequence: str, qualities: str) -> Tuple[str, str]:
    out_sequence = sequence.translate(complement)[::-1]  # Reverse the sequence
    out_quality = qualities[::-1]  # Reverse the quality scores
    return out_sequence, out_quality


def filter_sequence(record, forward_primer, reverse_primer):
    primer_sites = {}

    primer_f = forward_primer
    primer_f_rc = forward_primer.translate(complement)[::-1]
    forward_alignment = None
    reverse_alignment = None

    forward_alignment = edlib.align(primer_f, record.sequence, mode="HW", task="path", additionalEqualities=[("R", "A"), ("R", "G")])
    if forward_alignment["editDistance"] >2:
        forward_alignment = None
    reverse_alignment = edlib.align(primer_f_rc, record.sequence, mode="HW", task="path", additionalEqualities=[("R", "A"), ("R", "G")])
    if reverse_alignment["editDistance"] >2:
        reverse_alignment = None
    if forward_alignment != None and reverse_alignment != None:
        return record, "failed"
    elif forward_alignment == None and reverse_alignment == None:
        return record, "failed"
    elif forward_alignment != None:
        primer_sites["forward"] = [forward_alignment, "forward"]
    elif reverse_alignment != None:
        primer_sites["forward"] = [reverse_alignment, "reverse"]

    forward_alignment = None
    reverse_alignment = None
    primer_r = reverse_primer
    primer_r_rc = reverse_primer.translate(complement)[::-1]
    forward_alignment = edlib.align(primer_r, record.sequence, mode="HW", task="path", additionalEqualities=[("R", "A"), ("R", "G"), ("Y", "C"), ("Y", "T"), ("S", "G"), ("S", "C"), ("W", "A"), ("W", "T"), ("K", "G"), ("K", "T"), ("M", "A"), ("M", "C"), ("B", "C"), ("B", "G"), ("B", "T"), ("V", "A"), ("V", "C"), ("V", "G"), ("D", "A"), ("D", "G"), ("D", "T"), ("H", "A"), ("H", "C"), ("H", "T"), ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")])
    if forward_alignment["editDistance"] >2:
        forward_alignment = None
    reverse_alignment = edlib.align(primer_r_rc, record.sequence, mode="HW", task="path", additionalEqualities=[("R", "A"), ("R", "G"), ("Y", "C"), ("Y", "T"), ("S", "G"), ("S", "C"), ("W", "A"), ("W", "T"), ("K", "G"), ("K", "T"), ("M", "A"), ("M", "C"), ("B", "C"), ("B", "G"), ("B", "T"), ("V", "A"), ("V", "C"), ("V", "G"), ("D", "A"), ("D", "G"), ("D", "T"), ("H", "A"), ("H", "C"), ("H", "T"), ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")])
    if reverse_alignment["editDistance"] >2:
        reverse_alignment = None
    if forward_alignment != None and reverse_alignment != None:
        return record, "failed"
    elif forward_alignment == None and reverse_alignment == None:
        return record, "failed"
    elif forward_alignment != None:
        primer_sites["reverse"] = [forward_alignment, "forward"]
    elif reverse_alignment != None:
        primer_sites["reverse"] = [reverse_alignment, "reverse"]
    
    if primer_sites["forward"][1] != primer_sites["reverse"][1]:
        if primer_sites["forward"][0]["locations"][0][0] < primer_sites["reverse"][0]["locations"][0][1]:
            pass
        elif primer_sites["forward"][0]["locations"][0][0] > primer_sites["reverse"][0]["locations"][0][1]:
            record.sequence, record.qualities = get_reverse_complement(record.sequence, record.qualities)
        return record, "success"
    else:
        return record, "failed"

def process_chunk(input_chunk, forward_primer, reverse_primer):
    output_chunk = defaultdict(list)
    with io.BytesIO(input_chunk) as temp_file: 
        for record in dnaio.open(temp_file):
            record, destination = filter_sequence(record, forward_primer, reverse_primer)
            output_chunk[destination].append(record)
    return output_chunk

def main():
    args = parse_arguments()
    forward_primer = args.fwd_seq
    reverse_primer = args.rev_seq
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
    # To ensure that the output files are named correctly, we need to remove the file extension from the input file. To eliminate "fq.gz" or ".fastq" etc, we split the basename by "." and take the first element
    input_basename = os.path.basename(args.input).split(".")
    if len(input_basename) >=3:
        input_basename = ".".join(input_basename[:-2])
    else:
        input_basename = ".".join(input_basename[:-1])
    num_workers = num_threads
    with xopen(args.input, mode="rb", threads=num_threads) as input_file:
        number_of_saved_reads = 0
        output_format = _detect_format_from_content(input_file)
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = []
            open_flag = defaultdict(bool) 
            for chunk in dnaio.read_chunks(input_file, buffer_size=4194304):
                future = executor.submit(process_chunk, chunk.tobytes(), forward_primer, reverse_primer)
                futures.append(future)
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                for key in result.keys():
                    if open_flag[key] == False:
                        open_mode = "w"
                        open_flag[key] = True
                    else:
                        open_mode = "a"
                    with dnaio.open(f"{input_basename}_{key}.{output_format}.gz", mode=open_mode, compression_level=args.compression_level, open_threads=dnaio_open_threads, fileformat=output_format) as output_handle:
                        for record in result[key]:
                            output_handle.write(record)
                            number_of_saved_reads += 1
        logging.info(f"A total of {number_of_saved_reads} reads were saved")

if __name__ == "__main__":
    main()