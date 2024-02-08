#!/usr/bin/env python
# coding: utf-8

import logging
import os
import sys
import argparse
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Setup the logging system. Configures a stream handler to output log messages to stdout.
# Default logging level is set to INFO.

logger = logging.getLogger()
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
logger.addHandler(handler)

def _get_args():
    """
    """
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument(
        '-i',
        '--input',
        help='FASTA file (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        help='output file prefix (default: basename of the file + trimmed)',
        type=str)
    parser.add_argument(
        '--fwd_name',
        help='forward primer name (optional)',
        type=str,
        default="")
    parser.add_argument(
        '--rev_name',
        help='reverse primer name (optional)',
        type=str,
        default="")
    parser.add_argument(
        '-f',
        '--fwd_seq',
        help='forward primer sequence (required)',
        type=str,
        required=True,
        default="")
    parser.add_argument(
        '-r',
        '--rev_seq',
        help='reverse primer sequence (required)',
        type=str,
        required=True,
        default="")
    parser.add_argument(
        '-m',
        '--min',
        help='Minimium output sequence length (default: 500)',
        type=int,
        default=500)
    parser.add_argument(
        '--keep_primers',
        help='Keep primer sequences from the resulting sequences (default: False).',
        action='store_true')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    args = parser.parse_args()
    return args

aligner = Align.PairwiseAligner()
aligner.mode = "local" 
aligner.match_score = 2
aligner.mismatch_score = -3
aligner.open_gap_score = -10
aligner.extend_gap_score = -2
aligner.query_end_gap_score = -1  # Slightly penalize end gaps to allow some flexibility but still discourage them
aligner.target_end_gap_score = -1  # Similar reasoning for the target sequence

def get_pairwise_identity(query: str, subject: str, aligner):
    # https://www.biostars.org/p/9553992/
    # Convert strings to Seq objects, necessary for reverse complement
    query_seq = query.seq
    query_rc_seq = query_seq.reverse_complement()    
    subject_seq = subject.seq

    
    # Perform forward alignment
    alignments_f = aligner.align(query_seq, subject_seq)
    alignment_f = alignments_f[0]
    # Perform reverse complement alignment
    alignments_r = aligner.align(query_rc_seq, subject_seq)
    alignment_r = alignments_r[0]
    # Select final alignment based on higher alignment score
    if alignment_f.score >= alignment_r.score:
        final_alignment = alignment_f
    else:
        final_alignment = alignment_r
    # Extract aligned sequences from the final alignment
    aln1, aln2 = final_alignment
    
    # Calculate identities
    identities = sum(seq1 == seq2 for seq1, seq2 in zip(aln1, aln2) if seq1 != '-' and seq2 != '-')
    
    # Calculate identity percentage
    identity = 100 * identities / len(query_seq)
    # Calculate query coverage
    query_coverage = 100 * (final_alignment.aligned[0][-1][-1] - final_alignment.aligned[0][0][0]) / len(query)

    if final_alignment == alignment_f:
        subject_start_coordinate = final_alignment.aligned[1][0][0]
        subject_end_coordinate = final_alignment.aligned[1][-1][-1]
        query_start_coordinate = final_alignment.aligned[0][0][0]
        query_end_coordinate = final_alignment.aligned[0][-1][-1]
        subject_match_sequence = subject_seq[subject_start_coordinate:subject_end_coordinate]
    else:
        subject_start_coordinate = final_alignment.aligned[1][-1][-1]
        subject_end_coordinate = final_alignment.aligned[1][0][0]
        query_start_coordinate = len(query_seq) - final_alignment.aligned[0][-1][-1]
        query_end_coordinate = len(query_seq) - final_alignment.aligned[0][0][0]
        subject_match_sequence = subject_seq[subject_end_coordinate:subject_start_coordinate].reverse_complement()
    return round(identity, 2), round(query_coverage, 2), query_start_coordinate, query_end_coordinate, subject_start_coordinate, subject_end_coordinate, subject_match_sequence


def trim_contig(contig, result_fwd, result_rev, keep_primers):
    if keep_primers == True:
        trimmed_start = result_fwd[4]
        trimmed_end = result_rev[4]
    else:
        trimmed_start = result_fwd[5]
        trimmed_end = result_rev[5]
    trimmed_contig = contig.seq[trimmed_start:trimmed_end]
    return trimmed_contig

def determine_primer_orientations(result_fwd, result_rev):
    primer_orientation = ""
    fwd_start = result_fwd[4]
    fwd_end = result_fwd[5]
    rev_start = result_rev[4]
    rev_end = result_rev[5]
    
    
    if fwd_start < fwd_end:
        fwd_is_sense = True
    else:
        fwd_is_sense = False
        
    if rev_start < rev_end:
        rev_is_sense = True
    else:
        rev_is_sense = False

    if (fwd_is_sense, rev_is_sense) == (True, False):
        if fwd_end < rev_end:
            primer_orientation = "face_to_face"
        else:
            if rev_end <= fwd_start and fwd_start < rev_start:
                primer_orientation = "dovetail"
            else:
                primer_orientation = "back_to_back"
    return primer_orientation


def determine_output_file_prefix(in_fa, out_fa):
    """
    """
    if out_fa is not None:
        return out_fa
    else:
        in_fa_base: str = os.path.splitext(in_fa)[0]
        return "{}.trimmed.fa".format(in_fa_base)

def process_contig(contig, primer_fwd, primer_rev, keep_primers, min_len):
    result_fwd = get_pairwise_identity(query=primer_fwd, subject=contig, aligner=aligner)
    result_rev = get_pairwise_identity(query=primer_rev, subject=contig, aligner=aligner)
    primer_orientation = determine_primer_orientations(result_fwd, result_rev)
    if primer_orientation == "face_to_face":
        trimmed_seq = trim_contig(contig, result_fwd, result_rev, keep_primers)  
    else:
        return None
    if len(trimmed_seq) >= min_len:
        has_valid_amplicon = True
    else:
        has_valid_amplicon = False

    if has_valid_amplicon == True:
        description = 'fwd_name:{},fwd_seq:{},rev_name:{},rev_seq:{}'.format(primer_fwd.id, str(primer_fwd.seq), primer_rev.id, str(primer_rev.seq))
        trimmed_record = SeqRecord(trimmed_seq, id=contig.id,description=description)
        return trimmed_record
    else:
        return None

def main():
    args = _get_args()
    in_fa = args.input
    out_fa = args.output
    fwd_name = args.fwd_name
    rev_name = args.rev_name
    fwd_seq = args.fwd_seq
    rev_seq = args.rev_seq
    min_len = args.min
    keep_primers = args.keep_primers
    
    out_fa = determine_output_file_prefix(in_fa, out_fa)
    
    out_records = []
    
    primer_fwd = SeqRecord(Seq(fwd_seq),id=fwd_name)
    primer_rev = SeqRecord(Seq(rev_seq),id=rev_name)
    logger.info("INFO: Forward primer name: {}".format(fwd_name))
    logger.info("INFO: Forward primer sequence: {}".format(fwd_seq))
    logger.info("INFO: Reverse primer name: {}".format(rev_name))
    logger.info("INFO: Reverse primer sequence: {}".format(rev_seq))    
    contigs = list(SeqIO.parse(in_fa, "fasta"))
    logger.info("INFO: Input FASTA file: {}".format(in_fa))
    logger.info("INFO: Number of sequences: {}".format(len(contigs)))

    for contig in contigs:
        trimmed_record = process_contig(contig, primer_fwd, primer_rev, keep_primers, min_len)
        if trimmed_record is not None:
            out_records.append(trimmed_record)
            logger.info("INFO: Valid amplicon was found for {}".format(trimmed_record.id))    
        else:
            logger.warning("WARNING: No valid amplicon was found for {}".format(contig.id))
        
    SeqIO.write(out_records, out_fa, 'fasta')
    logger.info("INFO: Trimmed sequences have been saved as {}".format(out_fa))

if __name__ == "__main__":
    main()
