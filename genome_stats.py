#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
import os
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

def _get_args():
    parser = argparse.ArgumentParser(description='Produce a simple TSV summary file for microbial genome assemblies')
    parser.add_argument('-i','--input',nargs='*', help='Genbank/DDBJ flatfile (required)',type=str,required=True)
    parser.add_argument("--output","-o", "--out","--output",metavar="FILE",help="output TSV file",default="out.tsv")
    args = parser.parse_args()
    return args

def load_gbks(gbk_list):
    record_dict = dict()
    for gbk_file in gbk_list:
        record_list = []
        records = SeqIO.parse(gbk_file, 'genbank')
        for record in records:
            record_list.append(record)
        record_dict[gbk_file] = record_list
    return record_dict

def calculate_gc_percent(sequence):
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    gc_percent = 100 * ((g_count + c_count) / len(sequence))
    gc_percent = round(gc_percent, 2)
    return gc_percent

def extract_assembly_and_strain_names(first_record):
    dbxrefs_dict = dict(map(lambda s: s.split(':'), first_record.dbxrefs))
    assembly_name = dbxrefs_dict.get('Assembly', first_record.id)

    strain_name = "unknown"
    for feature in first_record.features:
        if feature.type == "source":
            if 'isolate' in feature.qualifiers:
                strain_name = feature.qualifiers['isolate'][0]
            elif 'strain' in feature.qualifiers:
                strain_name = feature.qualifiers['strain'][0]

    return assembly_name, strain_name

def count_feature(records, feature_type):
    return sum(1 for record in records for feature in record.features if feature.type == feature_type)

def count_crispr(records):
    crispr_count = 0
    for record in records:
        for feature in record.features:
            if feature.type == 'repeat_region' and 'rpt_family' in feature.qualifiers:
                if feature.qualifiers['rpt_family'][0] == 'CRISPR':
                    crispr_count += 1
    return crispr_count

def aggregate_sequence_and_calculate_gc(records):
    whole_seq = ''.join(str(record.seq) for record in records)
    gc_content = calculate_gc_percent(whole_seq)
    return len(whole_seq), gc_content

def count_unique_trnas(records):
    trna_products = set()
    for record in records:
        for feature in record.features:
            if feature.type == 'tRNA' and 'product' in feature.qualifiers:
                trna_products.add(feature.qualifiers['product'][0])
    return len(trna_products)

def get_assembly_stats(records):
    first_record = records[0]
    assembly_name, strain_name = extract_assembly_and_strain_names(first_record)
    sequence_length, gc_content = aggregate_sequence_and_calculate_gc(records)
    unique_trna_count = count_unique_trnas(records)

    out_dict = {
        'Species': first_record.annotations.get('organism', 'Unknown'),
        'strain': strain_name,
        'Accession': assembly_name,
        'Length': sequence_length,
        'num_contigs': len(records),
        'GC%': gc_content,
        'CDS': count_feature(records, 'CDS'),
        'rRNA': count_feature(records, 'rRNA'),
        'tRNA': count_feature(records, 'tRNA'),
        'Unique tRNAs': unique_trna_count,
        'CRISPR': count_crispr(records)
    }

    return out_dict

def main():
    args = _get_args()
    in_gb = args.input
    out_tsv = args.output
    record_dict = load_gbks(in_gb)
    assembly_dict = dict()
    for key in record_dict.keys():
        out_dict = get_assembly_stats(record_dict[key])
        assembly_dict[key] = out_dict
    df = pd.DataFrame.from_dict(assembly_dict, orient='columns')
    df.to_csv(out_tsv, sep='\t', header=False)

if __name__ == "__main__":
    main()
