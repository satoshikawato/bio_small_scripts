#!/usr/bin/env python
# coding: utf-8

"""
Convert DDBJ/GenBank flatfile into gff3.
Args:
     (XX): GenBank or DDBJ format sequence file with annotation.
"""
import argparse
import collections
import logging
import pathlib
import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from logging import getLogger
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def _get_args():
    parser = argparse.ArgumentParser(
        description='Convert DDBJ/GenBank flatfile into gff3')
    parser.add_argument(
        '-i',
        '--input',
        help='input GenBank/DDBJ flatfile (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        help=' output gff3 file (default:out.gff3)',
        type=str)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def load_gbks(gbk_list):
    """Load the first instance of GenBank as SerRecord
    Args:
    """
    record_list = []
    for gbk_file in gbk_list:
        logger.info("Loading GenBank...")
        records = SeqIO.parse(gbk_file, 'genbank')
        logger.info("             ... done loading GenBank")
        logger.info("Input GenBank       : {}".format(gbk_file))
        record_list.append(record)
    return record_list


def gbk_to_seqrecords(gbk_file):
    """Load the first instance of GenBank as SerRecord
    Args:
    """
    records = []
    for record in SeqIO.parse(gbk_file, 'genbank'):
        logger.info("Loading SeqRecord...")
        records.append(record)
    logger.info("             ... done loading SeqRecords")
    return records


class GeneObject:
    def __init__(self, gene_id, locations, qualifiers, gene_biotype):
        self.gene_id = gene_id
        self.locations = locations
        self.qualifiers = qualifiers
        self.gene_biotype = gene_biotype


class RepeatObject:
    def __init__(self, repeat_id, locations, qualifiers):
        self.repeat_id = repeat_id
        self.locations = locations
        self.qualifiers = qualifiers


class FeatureObject:
    def __init__(self, feature_id, locations, qualifiers):
        self.feature_id = feature_id
        self.locations = locations
        self.qualifiers = qualifiers


def create_gene_object(locus_id, feature):
    coordinates = get_coordinates(feature.location.parts)
    qualifiers = feature.qualifiers
    if feature.type == 'CDS':
        gene_biotype = "protein_coding"
    elif feature.type == 'rRNA':
        gene_biotype = "rRNA"
    elif feature.type == 'tRNA':
        gene_biotype = "tRNA"
    else:
        gene_biotype = "other"
    gene_object = GeneObject(locus_id, coordinates, qualifiers, gene_biotype)
    return gene_object


def create_repeat_object(repeat_id, feature):
    coordinates = get_coordinates(feature.location.parts)
    qualifiers = feature.qualifiers
    repeat_object = RepeatObject(repeat_id, coordinates, qualifiers)
    return repeat_object


def create_feature_object(feature_id, feature):
    coordinates = get_coordinates(feature.location.parts)
    qualifiers = feature.qualifiers
    feature_object = FeatureObject(feature_id, coordinates, qualifiers)
    return feature_object


def create_feature_dict(gb_record):
    feature_dict = {}
    cds_count = 0
    rrna_count = 0
    trna_count = 0
    repeat_count = 0
    feature_count = 0
    for feature in gb_record.features:
        if feature.type == 'CDS':
            cds_count = cds_count + 1
            locus_id = "LOCUS_" + str(cds_count).zfill(3) + "0"
            gene_object = create_gene_object(locus_id, feature)
            feature_dict[locus_id] = gene_object

        elif feature.type == 'rRNA':
            rrna_count = rrna_count + 1
            locus_id = "LOCUS_r" + str(rrna_count).zfill(3) + "0"
            gene_object = create_gene_object(locus_id, feature)
            feature_dict[locus_id] = gene_object

        elif feature.type == 'tRNA':
            trna_count = trna_count + 1
            locus_id = "LOCUS_t" + str(trna_count).zfill(3) + "0"
            gene_object = create_gene_object(locus_id, feature)
            feature_dict[locus_id] = gene_object

        elif feature.type == 'repeat_region':
            repeat_count = repeat_count + 1
            repeat_id = str(repeat_count).zfill(3) + "0"
            repeat_object = create_repeat_object(repeat_id, feature)
            feature_dict[repeat_id] = repeat_object
        elif feature.type == 'misc_feature':
            feature_count = feature_count + 1
            feature_id = str(feature_count).zfill(3) + "0"
            feature_object = create_feature_object(feature_id, feature)
            feature_dict[feature_id] = feature_object
    return feature_dict


def get_coordinates(feature_location_parts: list) -> list:
    coordinates = []
    for part in feature_location_parts:
        if part.strand == 1:
            strand = "+"
        elif part.strand == -1:
            strand = "-"
        else:
            strand = "."
        start = int(part.start)
        end = int(part.end)
        coordinate = [strand, start, end]
        coordinates.append(coordinate)
    return coordinates


def calc_phase(exon_len_of_interest):
    if exon_len_of_interest % 3 == 0:
        phase = 0
    elif exon_len_of_interest % 3 == 1:
        phase = 2
    elif exon_len_of_interest % 3 == 2:
        phase = 1
    return phase


def main():
    args = _get_args()
    gbk_file = args.input
    out_gff = args.output
    records = gbk_to_seqrecords(gbk_file)
    source = "ddbj_to_gff3.py"
    score = "."
    out_lines = []
    comments = ''
    comments += '{}\n'.format("##gff-version 3")
    for record in records:
        seqid = record.name
        feature_dict = create_feature_dict(record)
        comments += '##sequence-region {} {} {}\n'.format(
            seqid, 1, len(record.seq))
        out_features = ''
        for feature_id in feature_dict.keys():
            feature_object = feature_dict[feature_id]
            locations = feature_object.locations
            feature_strand = locations[0][0]
            feature_start = str(min([locations[i][1]
                                for i in range(len(locations))]) + 1)
            feature_end = str(max([locations[i][2]
                              for i in range(len(locations))]))
            if isinstance(feature_object, GeneObject):
                locus_id = feature_id
                gene_object = feature_object
                gene_biotype = gene_object.gene_biotype
                product = gene_object.qualifiers['product'][0]
                if 'note' in gene_object.qualifiers.keys():
                    notes = gene_object.qualifiers['note']
                else:
                    notes = ""
                new_notes = ""
                for note in notes:
                    if "source:" in note:
                        continue
                    else:
                        new_notes += "note={};".format(note)
                        break
                if gene_biotype == 'rRNA':
                    gene_qualifiers = "ID={};product={};gene_biotype={};".format(
                        locus_id, product, gene_biotype)
                    gene_line = [
                        seqid,
                        source,
                        "gene",
                        feature_start,
                        feature_end,
                        score,
                        feature_strand,
                        ".",
                        gene_qualifiers]
                    gene_line = '\t'.join(gene_line)
                    out_features += '{}\n'.format(gene_line)
                    rna_id = locus_id.replace('LOCUS_r', 'rrna-')
                    new_qualifiers = "Parent={};ID={};product={};".format(
                        locus_id, rna_id, new_notes, product)
                    for location in locations:
                        phase = "."
                        strand = str(location[0])
                        start = str(location[1] + 1)
                        end = str(location[2])
                        out_feature_line = [
                            seqid,
                            source,
                            "rRNA",
                            start,
                            end,
                            score,
                            strand,
                            phase,
                            new_qualifiers]
                        out_feature_line = '\t'.join(out_feature_line)
                        out_features += '{}\n'.format(out_feature_line)
                elif gene_biotype == 'tRNA':
                    gene_qualifiers = "ID={};product={};gene_biotype={};".format(
                        locus_id, product, gene_biotype)
                    gene_line = [
                        seqid,
                        source,
                        "gene",
                        feature_start,
                        feature_end,
                        score,
                        feature_strand,
                        ".",
                        gene_qualifiers]
                    gene_line = '\t'.join(gene_line)
                    out_features += '{}\n'.format(gene_line)
                    rna_id = locus_id.replace('LOCUS_t', 'trna-')
                    new_qualifiers = "Parent={};ID={};{}product={};".format(
                        locus_id, rna_id, new_notes, product)
                    for location in locations:
                        phase = "."
                        strand = str(location[0])
                        start = str(location[1] + 1)
                        end = str(location[2])
                        out_feature_line = [
                            seqid,
                            source,
                            "tRNA",
                            start,
                            end,
                            score,
                            strand,
                            phase,
                            new_qualifiers]
                        out_feature_line = '\t'.join(out_feature_line)
                        out_features += '{}\n'.format(out_feature_line)
                elif gene_biotype == 'protein_coding':
                    gene_qualifiers = "ID={};product={};gene_biotype={};".format(
                        locus_id, product, gene_biotype)
                    gene_line = [
                        seqid,
                        source,
                        "mRNA",
                        feature_start,
                        feature_end,
                        score,
                        feature_strand,
                        ".",
                        gene_qualifiers]
                    gene_line = '\t'.join(gene_line)
                    out_features += '{}\n'.format(gene_line)
                    cds_id = locus_id.replace('LOCUS_', 'cds-')
                    new_qualifiers = "Parent={};ID={};{}product={};".format(
                        locus_id, cds_id, new_notes, product)
                    phases = []
                    phase = int(gene_object.qualifiers['codon_start'][0]) - 1
                    phases.append(phase)
                    if len(locations) > 1:
                        for i in range(1, len(locations)):
                            exon_len_of_interest = (
                                locations[i - 1][2] - locations[i - 1][1] - phases[i - 1])
                            next_phase = calc_phase(exon_len_of_interest)
                            phases.append(next_phase)
                    for i in range(len(locations)):
                        phase = str(phases[i])
                        strand = str(locations[i][0])
                        start = str(locations[i][1] + 1)
                        end = str(locations[i][2])
                        out_feature_line = [
                            seqid,
                            source,
                            "CDS",
                            start,
                            end,
                            score,
                            strand,
                            phase,
                            new_qualifiers]
                        out_feature_line = '\t'.join(out_feature_line)
                        out_features += '{}\n'.format(out_feature_line)
            elif isinstance(feature_object, RepeatObject):
                qualifiers = feature_object.qualifiers
                new_qualifiers = ""
                for key in qualifiers.keys():
                    values = qualifiers[key]
                    for value in values:
                        new_qualifiers += "{}={};".format(key, value)
                for location in locations:
                    phase = "."
                    strand = str(location[0])
                    start = str(location[1] + 1)
                    end = str(location[2])
                    out_feature_line = [
                        seqid,
                        source,
                        "repeat_region",
                        start,
                        end,
                        score,
                        strand,
                        phase,
                        new_qualifiers]
                    out_feature_line = '\t'.join(out_feature_line)
                    out_features += '{}\n'.format(out_feature_line)
            elif isinstance(feature_object, FeatureObject):
                qualifiers = feature_object.qualifiers
                new_qualifiers = ""
                for key in qualifiers.keys():
                    values = qualifiers[key]
                    for value in values:
                        new_qualifiers += "{}={};".format(key, value)
                for location in locations:
                    phase = "."
                    strand = str(location[0])
                    start = str(location[1] + 1)
                    end = str(location[2])
                    out_feature_line = [
                        seqid,
                        source,
                        "misc_feature",
                        start,
                        end,
                        score,
                        strand,
                        phase,
                        new_qualifiers]
                    out_feature_line = '\t'.join(out_feature_line)
                    out_features += '{}\n'.format(out_feature_line)
        with open(out_gff, "w") as f:
            print(comments, file=f, end='')
            print(out_features, file=f, end='')

if __name__ == "__main__":
    main()
