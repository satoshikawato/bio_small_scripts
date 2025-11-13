#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def _get_args():
    parser = argparse.ArgumentParser(
        description='Extract the longest isoforms of protein-coding genes from a NCBI RefSeq eukaryotic genome assembly')
    parser.add_argument(
        '-i', '--input', help='GenBank flat file format of the genomic sequence(s) (required)', type=str, required=True)
    parser.add_argument(
        '-o', '--output', help='output fasta-formatted file (optional)', type=str, default="out.faa")
    parser.add_argument(
        '-t', '--table', help='output tab-separated table (optional)', type=str, default="out.tbl")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def _process_record_and_write(gb_record, out_handle, tbl_handle):
    """
    Extract the longest isoforms from a single GenBank record and write them to files immediately.
    """
    # {gene_id: (length, seq_record)}
    best_proteins = {}

    organism = gb_record.annotations.get('organism', 'unknown_organism')

    for feature in gb_record.features:
        if feature.type != 'CDS':
            continue

        qualifiers = feature.qualifiers
        
        if 'protein_id' not in qualifiers or 'gene' not in qualifiers or 'translation' not in qualifiers:
            continue

        protein_id = qualifiers['protein_id'][0]
        gene_id = qualifiers['gene'][0]
        product_name = qualifiers.get('product', ['unknown_product'])[0]
        translation_seq = qualifiers['translation'][0]
        
        seq_len = len(translation_seq)

        if gene_id not in best_proteins or seq_len > best_proteins[gene_id][0]:
            aa_seq = Seq(translation_seq)
            aa_record = SeqRecord(
                aa_seq, 
                id=protein_id,
                name=gene_id, 
                description="{} [{}]".format(product_name, organism)
            )
            best_proteins[gene_id] = (seq_len, aa_record)

    if best_proteins:
        records_to_write = [item[1] for item in best_proteins.values()]

        SeqIO.write(records_to_write, out_handle, "fasta")
        
        for record in records_to_write:
            tbl_handle.write('{}\t{}\n'.format(record.name, record.id))

def main():
    args = _get_args()
    in_gb = args.input
    out_file = args.output
    out_table = args.table

    with open(out_file, "w") as out_handle, open(out_table, "w") as tbl_handle:
        for gb_record in SeqIO.parse(in_gb, 'genbank'):
            _process_record_and_write(gb_record, out_handle, tbl_handle)

if __name__ == "__main__":
    main()
