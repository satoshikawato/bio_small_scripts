#!/usr/bin/env python3

import argparse
import sys
import os
from Bio import SeqIO, Seq
from Bio.Data import CodonTable
import tqdm
import logging

logging.basicConfig(level=logging.INFO, stream=sys.stderr)

def calculate_gc_content(sequence):
    """Compute the GC content of a given sequence."""
    return (sequence.count('G') + sequence.count('C')) / max(1, len(sequence))


def heuristic_back_translate(sequence, protein_to_codon, desired_gc):
    nucleotide_seq = ""

    for aa in sequence:
        if aa not in protein_to_codon:
            continue

        current_gc = calculate_gc_content(nucleotide_seq)
        codons = sorted(protein_to_codon[aa], key=calculate_gc_content)

        if current_gc < desired_gc:
            nucleotide_seq += codons[-1]
        else:
            nucleotide_seq += codons[0]

    resulting_gc = calculate_gc_content(nucleotide_seq)
    return nucleotide_seq, resulting_gc
    
def calculate_min_max_gc(sequence, protein_to_codon):
    min_gc_seq = ''.join(sorted(protein_to_codon[aa], key=calculate_gc_content)[0] for aa in sequence if aa in protein_to_codon)
    max_gc_seq = ''.join(sorted(protein_to_codon[aa], key=calculate_gc_content, reverse=True)[0] for aa in sequence if aa in protein_to_codon)

    return calculate_gc_content(min_gc_seq), calculate_gc_content(max_gc_seq)

def process_record(record, protein_to_codon, desired_gc, file_format, output_file, gc_adjusted_count, tsv_report):
    min_gc, max_gc = calculate_min_max_gc(str(record.seq), protein_to_codon)
    current_desired_gc = desired_gc
    note = ""

    if current_desired_gc < min_gc:
        current_desired_gc = min_gc
        note = "desired GC% unreachable; adjusted to theoretical min"
        gc_adjusted_count[0] += 1
    elif current_desired_gc > max_gc:
        current_desired_gc = max_gc
        note = "desired GC% unreachable; adjusted to theoretical max"
        gc_adjusted_count[0] += 1

    translated_seq, resulting_gc = heuristic_back_translate(str(record.seq), protein_to_codon, current_desired_gc)

    record.seq = Seq.Seq(translated_seq)
    record.description += f" | Resulting GC: {resulting_gc*100:.2f}%"

    with open(output_file, 'a') as outfile:
        SeqIO.write(record, outfile, file_format)

    tsv_report.append([record.id, len(record.seq), min_gc * 100, max_gc * 100, resulting_gc * 100, note])

    return min_gc, max_gc, resulting_gc

def main(input_file, output_file, file_format, genetic_table_id, desired_gc, report_file):
    table = CodonTable.unambiguous_dna_by_id[genetic_table_id]
    protein_to_codon = {value: [key for key in table.forward_table if table.forward_table[key] == value] for value in table.forward_table.values()}

    sum_min_gc = 0
    sum_max_gc = 0
    sum_actual_gc = 0
    total_proteins = 0
    gc_adjusted_count = [0]
    tsv_report = []

    records = list(SeqIO.parse(input_file, file_format))

    for record in tqdm.tqdm(records, desc="Processing records"):
        min_gc, max_gc, actual_gc = process_record(record, protein_to_codon, desired_gc/100, file_format, output_file, gc_adjusted_count, tsv_report)

        sum_min_gc += min_gc
        sum_max_gc += max_gc
        sum_actual_gc += actual_gc
        total_proteins += 1

    avg_min_gc = (sum_min_gc / total_proteins) * 100
    avg_max_gc = (sum_max_gc / total_proteins) * 100
    avg_actual_gc = (sum_actual_gc / total_proteins) * 100

    logging.info(f"Summary:\nTotal proteins back translated: {total_proteins}\nAverage theoretical minimum GC: {avg_min_gc:.2f}%\nAverage theoretical maximum GC: {avg_max_gc:.2f}%\nAverage actual GC% of the backtranslated sequences: {avg_actual_gc:.2f}%")

    if gc_adjusted_count[0] > 0:
        percentage = (gc_adjusted_count[0] / total_proteins) * 100
        logging.warning(f"{gc_adjusted_count[0]} sequences ({percentage:.2f}%) were adjusted to fit the theoretical GC bounds.")

    with open(report_file, 'w') as rpt_file:
        rpt_file.write("id\tlength\ttheoretical_min_gc\ttheoretical_max_gc\tactual_gc\tnote\n")
        for line in tsv_report:
            rpt_file.write('\t'.join(map(str, line)) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Back-translate protein sequences into nucleotide sequences.")
    
    parser.add_argument('-i', '--input', required=True, help='Input protein sequence file.')
    
    default_output_name = None
    default_report_name = None
    args, _ = parser.parse_known_args()
    if args.input:
        base_name = os.path.splitext(os.path.basename(args.input))[0]
        default_output_name = f"{base_name}.fna"
        default_report_name = f"{base_name}_report.tsv"
    
    parser.add_argument('-o', '--output', default=default_output_name, help=f'Output nucleotide sequence file. Default: {default_output_name}')
    parser.add_argument('-f', '--format', default="fasta", help='File format (default: fasta).')
    parser.add_argument('-t', '--table', type=int, default=1, help='Genetic table ID (default: 1).')
    parser.add_argument('--gc', type=float, default=0.5, help='Desired GC content (default: 50%).')
    parser.add_argument('-r', '--report', default=default_report_name, help=f'Report file in TSV format detailing backtranslation results. Default: {default_report_name}')
    
    args = parser.parse_args()

    main(args.input, args.output, args.format, args.table, args.gc, args.report)
