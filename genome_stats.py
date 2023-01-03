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
    parser = argparse.ArgumentParser(description='Crop genbank file. ')
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

def get_assembly_stats(records):
    out_dict = dict()
    cds_count = 0
    rrna_count = 0
    trna_count = 0
    trna_list = []
    crispr_count = 0
    whole_seq = ''
    first_record = records[0]
    assembly_name = dict(map(lambda s : s.split(':'), first_record.dbxrefs))['Assembly']
    for feature in first_record.features:
        if feature.type == "source":
            if 'isolate' in feature.qualifiers.keys():
                strain_name = feature.qualifiers['isolate'][0]
            elif 'strain' in feature.qualifiers.keys():
                strain_name = feature.qualifiers['strain'][0]
            else:
                strain_name = "unknown"
    for record in records:
        whole_seq += record.seq
        for feature in record.features:
            if feature.type == 'CDS':
                cds_count += 1
            elif feature.type == 'rRNA':
                rrna_count += 1
            elif feature.type == 'tRNA':
                trna_count += 1
                trna_list.append(feature.qualifiers['product'][0])
            elif feature.type == 'repeat_region':
                if feature.qualifiers['rpt_family'][0] == 'CRISPR':
                    crispr_count += 1
    unique_trna = len(set(trna_list))
    out_dict['Species'] = first_record.annotations['organism']
    out_dict['strain'] = strain_name
    out_dict['Accession'] = assembly_name
    out_dict['Length'] = len(whole_seq)
    out_dict['num_contigs'] = len(records)
    out_dict['GC%'] = calculate_gc_percent(whole_seq)
    out_dict['CDS'] = cds_count
    out_dict['rRNA'] = rrna_count
    out_dict['tRNA'] = '{}({})'.format(trna_count, unique_trna)
    out_dict['CRISPR'] = crispr_count
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

