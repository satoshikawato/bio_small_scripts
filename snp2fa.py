#!/usr/bin/env python
# coding: utf-8

from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gzip
import os
import pandas as pd
import sys
import argparse
import logging
import pathlib

def _get_args():
    parser = argparse.ArgumentParser(description='Extract SNPs into a FASTA file on the basis of a VCF file')
    parser.add_argument("--vcf", "--input", "--in", "-i",  metavar="FILE", help="Input VCF file", required=True)
    parser.add_argument("-o", "--out", "--output",metavar="FILE",help="Specify the output FASTA file name")
    parser.add_argument("-b", "--bed", metavar="FILE",help="Include a set of sites on the basis of a BED file")
    parser.add_argument("-e", "--exclude_bed", metavar="FILE",type=str, help="Exclude a set of sites on the basis of a BED file")
    parser.add_argument("-m", "--max_sv_len", help="Max SV length to be considered (default: 10)", type=int, default="10")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def get_header(in_file):
    with open(in_file, "rt") as ifile:
        header = []
        for line in ifile:
            if line.startswith("#CHROM"):
                line = line.replace('\n','')
                header = [x for x in line.split('\t')]
                break
        return header

def get_strain_names(header):
    strains = header[9:]
    strains_out = [i.split('.')[0] for i in strains] 
    return strains_out

def csv_to_df(in_file, header_line):
    with open(in_file, "rt") as ifile:
        vcf = pd.read_csv(ifile, comment='#',  sep='\t', header=None, names=header_line)
        return vcf

def extract_snps(df, max_sv_len):
    lst = []
    for index,row in df.iterrows():
        len_list = []
        ref_seq = row["REF"]
        len_list.append(len(ref_seq))
        alt_seq = [x for x in row["ALT"].split(',')]
        len_list.append(len(max(alt_seq)))
        if max(len_list) <= max_sv_len:
            lst.append(row)
    df2 = pd.DataFrame(lst)
    return df2

def snp_to_seq(df_snp, strain_list):
    seq_dict = {}
    seq_records = []
    for strain in strain_list:
        seq_dict[strain] = []
    for index,row in df_snp.iterrows():
        var_dict = {}
        len_list = []
        ref_seq = row["REF"]
        var_dict["./."] = ref_seq
        len_list.append(len(ref_seq))
        alt_seq = [x for x in row["ALT"].split(',')]
        for i in range(len(alt_seq)):
            var_dict["{}/{}".format(i+1,i+1)] = alt_seq[i]
        len_list.append(len(max(alt_seq)))
        for strain in strain_list:
            column_name = "{}.bam".format(strain)
            genotype = str(row[column_name])
            strain_seq = var_dict[genotype]
            gap = "-" * (max(len_list) - len(strain_seq))
            seq_dict[strain].append(strain_seq + gap)
    
    for strain in strain_list:
        seq_record = SeqRecord(Seq(''.join(seq_dict[strain])), id=strain, description="")
        seq_records.append(seq_record)
    return seq_records

def main():
    args = _get_args()
    in_vcf = args.vcf
    output = args.out
    max_sv_len = args.max_sv_len
    header = get_header(in_vcf)
    strain_list  = get_strain_names(header)
    df = csv_to_df(in_vcf, header)
    df_snp = extract_snps(df, max_sv_len)
    out_records = snp_to_seq(df_snp, strain_list)    
    if output:
        parent_dir = pathlib.Path(output).parent
        parent_dir.mkdir(exist_ok=True)
        with open(output, "w") as f:
            SeqIO.write(out_records, output, "fasta")
    else:
        for record in out_records:
            print(record)

if __name__ == "__main__":
    main()

