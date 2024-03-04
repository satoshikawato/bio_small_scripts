#!/usr/bin/env python
# coding: utf-8

import sys
import logging
import argparse

# Setup the logging system. Configures a stream handler to output log messages to stdout.
# Default logging level is set to INFO.
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


class PafLine:
    def __init__(self, paf_line):
        self.fields = paf_line.strip().split('\t')
        self.query_name, self.query_length, self.query_start, self.query_end, self.relative_strand, self.target_name, self.target_length, self.target_start, self.target_end, self.num_residue_matches, self.block_length, self.mapping_quality, *_ = self.fields


def arg_parser():
    parser = argparse.ArgumentParser(description='Process PAF file and target reference list.')
    parser.add_argument('-i', '--input', help='Path to the PAF file')
    parser.add_argument('-t', '--target', help='Path to the target reference list file')
    parser.add_argument('-o', '--output', help='Path to the output file')
    parser.add_argument('-s', '--relative_strand', help='Relative strand ("plus", "minus", "both")', choices=['plus', 'minus', 'both'], default='both')
    parser.add_argument('-m', '--min_query_length', help='Minimum query length', type=int, default=0)
    parser.add_argument('-M', '--max_query_length', help='Maximum query length', type=int, default=None)
    parser.add_argument('-b', '--min_block_length', help='Minimum block length', type=int, default=0)
    parser.add_argument('-B', '--max_block_length', help='Maximum block length', type=int, default=None)
    parser.add_argument('-g', '--min_target_length', help='Minimum target length', type=int, default=0)
    parser.add_argument('-G', '--max_target_length', help='Maximum target length', type=int, default=None)
    parser.add_argument('-q', '--minimum_mapping_quality', help='Minimum mapping quality', type=int, default=0)
    parser.add_argument('-Q', '--maximum_mapping_quality', help='Maximum mapping quality', type=int, default=None)
    parser.add_argument('--min_identity', help='Minimum nucleotide identity (0-1)', type=float, default=0.0)
    parser.add_argument('--max_identity', help='Maximum nucleotide identity (0-1)', type=float, default=None)
    parser.add_argument('--min_query_target_ratio', help='Minimum query to target length ratio (0-1)', type=float, default=0.0)
    parser.add_argument('--max_query_target_ratio', help='Maximum query to target length ratio (0-1)', type=float, default=None)
    parser.add_argument('--min_block_query_ratio', help='Minimum block to query length ratio (0-1)', type=float, default=0.0)
    parser.add_argument('--max_block_query_ratio', help='Maximum block to query length ratio (0-1)', type=float, default=None)
    parser.add_argument('--min_block_target_ratio', help='Minimum block to target length ratio (0-1)', type=float, default=0.0)
    parser.add_argument('--max_block_target_ratio', help='Maximum block to target length ratio (0-1)', type=float, default=None)
    parser.add_argument('--minimuim_residue_matches', help='Minimum residue matches', type=int, default=0)
    parser.add_argument('--maximum_residue_matches', help='Maximum residue matches', type=int, default=None)
    parser.add_argument('--max_left_overhang', help='Maximum left overhang (0-1)', type=float, default=None)
    parser.add_argument('--max_right_overhang', help='Maximum right overhang (0-1)', type=float, default=None)
    return parser.parse_args()


def read_paf(in_file):
    with open(in_file, "r") as f:
        return [PafLine(line) for line in f]

def meet_identity_thresholds(paf_line, min_identity, max_identity):
    nucleotide_identity = int(paf_line.num_residue_matches) / int(paf_line.block_length)
    if max_identity is None:
        return min_identity <= nucleotide_identity
    return min_identity <= nucleotide_identity <= max_identity

def meet_min_identity_thresholds(paf_line, min_identity):
    if min_identity is None:
        return True
    return int(paf_line.num_residue_matches) / int(paf_line.block_length) >= min_identity

def meet_max_identity_thresholds(paf_line, max_identity):
    if max_identity is None:
        return True
    return int(paf_line.num_residue_matches) / int(paf_line.block_length) <= max_identity

def meet_min_block_length_thresholds(paf_line, min_block_length):
    if min_block_length is None:
        return True
    return int(paf_line.block_length) >= min_block_length

def meet_max_block_length_thresholds(paf_line, max_block_length):
    if max_block_length is None:
        return True
    return int(paf_line.block_length) <= max_block_length

def meet_min_query_target_ratio_thresholds(paf_line, min_query_target_ratio):
    query_target_ratio = int(paf_line.query_length) / int(paf_line.target_length)
    if min_query_target_ratio is None:
        return True
    return min_query_target_ratio <= query_target_ratio

def meet_max_query_target_ratio_thresholds(paf_line, max_query_target_ratio):
    query_target_ratio = int(paf_line.query_length) / int(paf_line.target_length)
    if max_query_target_ratio is None:
        return True
    return query_target_ratio <= max_query_target_ratio

def meet_min_query_length_thresholds(paf_line, min_query_length):
    if min_query_length is None:
        return True
    return int(paf_line.query_length) >= min_query_length

def meet_max_query_length_thresholds(paf_line, max_query_length):
    if max_query_length is None:
        return True
    return int(paf_line.query_length) <= max_query_length

def meet_min_target_length_thresholds(paf_line, min_target_length):
    if min_target_length is None:
        return True
    return int(paf_line.target_length) >= min_target_length

def meet_max_target_length_thresholds(paf_line, max_target_length):
    if max_target_length is None:
        return True
    return int(paf_line.target_length) <= max_target_length

def meet_min_mapping_quality_thresholds(paf_line, min_mapping_quality):
    if min_mapping_quality is None:
        return True
    return int(paf_line.mapping_quality) >= min_mapping_quality

def meet_max_mapping_quality_thresholds(paf_line, max_mapping_quality):
    if max_mapping_quality is None:
        return True
    return int(paf_line.mapping_quality) <= max_mapping_quality

def meet_min_block_query_ratio_thresholds(paf_line, min_block_query_ratio):
    block_query_ratio = int(paf_line.block_length) / int(paf_line.query_length)
    if min_block_query_ratio is None:
        return True
    return min_block_query_ratio <= block_query_ratio

def meet_max_block_query_ratio_thresholds(paf_line, max_block_query_ratio):
    block_query_ratio = int(paf_line.block_length) / int(paf_line.query_length)
    if max_block_query_ratio is None:
        return True
    return block_query_ratio <= max_block_query_ratio

def meet_min_block_target_ratio_thresholds(paf_line, min_block_target_ratio):
    block_target_ratio = int(paf_line.block_length) / int(paf_line.target_length)
    if min_block_target_ratio is None:
        return True
    return min_block_target_ratio <= block_target_ratio

def meet_max_block_target_ratio_thresholds(paf_line, max_block_target_ratio):
    block_target_ratio = int(paf_line.block_length) / int(paf_line.target_length)
    if max_block_target_ratio is None:
        return True
    return block_target_ratio <= max_block_target_ratio

def meet_min_residue_matches_thresholds(paf_line, min_residue_matches):
    if min_residue_matches is None:
        return True
    return min_residue_matches <= int(paf_line.num_residue_matches) 

def meet_max_residue_matches_thresholds(paf_line, max_residue_matches):
    if max_residue_matches is None:
        return True
    return int(paf_line.num_residue_matches) <= max_residue_matches 

def meet_max_left_overhang_thresholds(paf_line, max_left_overhang):
    left_overhang = int(paf_line.query_start) / int(paf_line.query_length)
    if max_left_overhang is None:
        return True
    return left_overhang <= max_left_overhang

def meet_max_right_overhang_thresholds(paf_line, max_right_overhang):
    right_overhang = 1 - (int(paf_line.query_end) / int(paf_line.query_length))
    if max_right_overhang is None:
        return True
    return right_overhang <= max_right_overhang

def meet_block_length_thresholds(paf_line, min_block_length, max_block_length):
    if all([
        meet_min_block_length_thresholds(paf_line, min_block_length),
        meet_max_block_length_thresholds(paf_line, max_block_length)
    ]):
        return True
    else:
        return False
def meets_query_length_thresholds(paf_line, min_query_length, max_query_length):
    if all([
        meet_min_query_length_thresholds(paf_line, min_query_length),
        meet_max_query_length_thresholds(paf_line, max_query_length)
    ]):
        return True
    else:
        return False
def meet_target_length_thresholds(paf_line, min_target_length, max_target_length):
    if all([
        meet_min_target_length_thresholds(paf_line, min_target_length),
        meet_max_target_length_thresholds(paf_line, max_target_length)
    ]):
        return True
    else:
        return False
    
def meet_query_target_ratio_thresholds(paf_line, min_query_target_ratio, max_query_target_ratio):
    if all([
        meet_min_query_target_ratio_thresholds(paf_line, min_query_target_ratio),
        meet_max_query_target_ratio_thresholds(paf_line, max_query_target_ratio)
    ]):
        return True
    else:
        return False

def meet_mapping_quality_thresholds(paf_line, min_mapping_quality, max_mapping_quality):
    if all([
        meet_min_mapping_quality_thresholds(paf_line, min_mapping_quality),
        meet_max_mapping_quality_thresholds(paf_line, max_mapping_quality)
    ]):
        return True
    
def meet_block_query_ratio_thresholds(paf_line, min_block_query_ratio, max_block_query_ratio):
    if all([
        meet_min_block_query_ratio_thresholds(paf_line, min_block_query_ratio),
        meet_max_block_query_ratio_thresholds(paf_line, max_block_query_ratio)
    ]):
        return True
    else:
        return False
    
def meet_block_target_ratio_thresholds(paf_line, min_block_target_ratio, max_block_target_ratio):
    if all([
        meet_min_block_target_ratio_thresholds(paf_line, min_block_target_ratio),
        meet_max_block_target_ratio_thresholds(paf_line, max_block_target_ratio)
    ]):
        return True
    else:
        return False
       
def meet_residue_matches_thresholds(paf_line, min_residue_matches, max_residue_matches):
    if all([
        meet_min_residue_matches_thresholds(paf_line, min_residue_matches),
        meet_max_residue_matches_thresholds(paf_line, max_residue_matches)
    ]):
        return True
    else:
        return False    
def meet_overhang_thresholds(paf_line, max_left_overhang, max_right_overhang):
    if all([
        meet_max_left_overhang_thresholds(paf_line, max_left_overhang),
        meet_max_right_overhang_thresholds(paf_line, max_right_overhang)
    ]):
        return True
    else:
        return False    

def meet_relative_strand_thresholds(paf_line, relative_strand):
    return relative_strand == "both" or paf_line.relative_strand == relative_strand

def is_a_target(paf_line, targets):
    if targets:
        return paf_line.target_name in targets
    else:
        return True

def filter_paf(paf_lines, args, targets):
    min_nucleotide_identity = args.min_identity
    min_block_target_ratio = args.min_block_target_ratio
    min_block_query_ratio = args.min_block_query_ratio
    min_residue_matches = args.minimuim_residue_matches
    min_mapping_quality = args.minimum_mapping_quality
    min_query_length = args.min_query_length
    min_target_length = args.min_target_length
    max_query_length = args.max_query_length
    max_target_length = args.max_target_length
    max_mapping_quality = args.maximum_mapping_quality
    max_residue_matches = args.maximum_residue_matches
    max_block_target_ratio = args.max_block_target_ratio
    max_block_query_ratio = args.max_block_query_ratio
    max_left_overhang = args.max_left_overhang
    max_right_overhang = args.max_right_overhang
    max_nucleotide_identity = args.max_identity
    max_block_length = args.max_block_length
    min_block_length = args.min_block_length
    relative_strand = args.relative_strand
    min_query_target_ratio = args.min_query_target_ratio
    max_query_target_ratio = args.max_query_target_ratio
    filtered_query_names = []
    for paf_line in paf_lines:
        # to check if the paf line meets all the conditions (if defined), add paf_line to filtered_query_names.
        # if no conditions are defined, add paf_line to filtered_query_names.
        # To do so, we will use the all() function to check if all the conditions are met.
        conditions = [
            ("is_a_target", is_a_target(paf_line, targets)),
            ("meet_identity_thresholds", meet_identity_thresholds(paf_line, min_nucleotide_identity, max_nucleotide_identity)),
            ("meet_block_length_thresholds", meet_block_length_thresholds(paf_line, min_block_length, max_block_length)),
            ("meets_query_length_thresholds", meets_query_length_thresholds(paf_line, min_query_length, max_query_length)),
            ("meet_target_length_thresholds", meet_target_length_thresholds(paf_line, min_target_length, max_target_length)),
            ("meet_mapping_quality_thresholds", meet_mapping_quality_thresholds(paf_line, min_mapping_quality, max_mapping_quality)),
            ("meet_block_query_ratio_thresholds", meet_block_query_ratio_thresholds(paf_line, min_block_query_ratio, max_block_query_ratio)),
            ("meet_block_target_ratio_thresholds", meet_block_target_ratio_thresholds(paf_line, min_block_target_ratio, max_block_target_ratio)),
            ("meet_residue_matches_thresholds", meet_residue_matches_thresholds(paf_line, min_residue_matches, max_residue_matches)),
            ("meet_overhang_thresholds", meet_overhang_thresholds(paf_line, max_left_overhang, max_right_overhang)),
            ("meet_relative_strand_thresholds", meet_relative_strand_thresholds(paf_line, relative_strand)),
            ("meet_query_target_ratio_thresholds", meet_query_target_ratio_thresholds(paf_line, min_query_target_ratio, max_query_target_ratio))
        ]
        # if all conditions are "True", add paf_line to filtered_query_names
        if all([condition[1] for condition in conditions]): # condition[1] is the value of the condition. condition[0] is the name of the condition.
            filtered_query_names.append(paf_line.query_name)
    filtered_query_names = list(set(filtered_query_names))

    return filtered_query_names


def read_lines(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]


def collate_references(paf_references, target_references):
    query_ids = []
    for reference in paf_references:
        if reference in target_references:
            query_ids.append(reference)
    return query_ids


def save_query_ids(query_ids, output_file):
    with open(output_file, 'w') as file:
        for query_id in query_ids:
            file.write(query_id + '\n')


def main():
    args = arg_parser()

    paf_file = args.input
    target_file = args.target
    output_file = args.output
    if target_file:
        targets = read_lines(target_file)
    else:
        targets = None
    paf_file = read_paf(paf_file)
    
    query_ids = filter_paf(paf_file, args, targets)
    save_query_ids(query_ids, output_file)


if __name__ == '__main__':
    main()
