#!/usr/bin/env python
# coding: utf-8

"""
Draw a linear genome diagram with introns.
Multi-track comparison (GenBank/DDBJ files + BLAST table).
"""
import argparse
import collections
import logging
import sys
import os
import re
import csv
import svgwrite
import pandas as pd

from svgwrite.container import Group
from svgwrite.shapes import Line
from svgwrite.path import Path
from svgwrite.text import Text

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIXML

from logging import getLogger

logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def _get_args():
    parser = argparse.ArgumentParser(description='Generate  plot in SVG')
    parser.add_argument(
        '-i',
        '--input',
        help='genbank (required)',
        type=str,
        required=True,
        nargs='*')
    parser.add_argument(
        '-b',
        '--blast',
        help="input BLAST result file in tab-separated format (-outfmt 6 or 7) (optional)",
        type=str,
        nargs='*')
    parser.add_argument(
        '-t',
        '--table',
        help='color table (optional)',
        type=str)
    parser.add_argument(
        '-o',
        '--output',
        help='output prefix (default: diagram)',
        type=str,
        default="diagram")
    parser.add_argument(
        '-n',
        '--nt',
        help='dinucleotide (default: GC). ',
        type=str,
        default="GC")
    parser.add_argument(
        '-w',
        '--window',
        help='window size (default: 1000) ',
        type=int,
        default="1000")
    parser.add_argument(
        '-s',
        '--step',
        help='step size (default: 100) ',
        type=int,
        default="100")
    parser.add_argument(
        '--separate_strands',
        help='separate forward and reverse strands (default: False). Features of undefined strands are shown on the forward strand. ',
        action='store_true')
    parser.add_argument(
        '--show_gc',
        help='plot GC content below genome (default: False). ',
        action='store_true')
    parser.add_argument(
        '--align_center',
        help='Align genomes to the center (default: False). ',
        action='store_true')
    parser.add_argument(
        '--evalue',
        help='evalue threshold (default=1e-2)',
        type=float,
        default="10")
    parser.add_argument(
        '--bitscore',
        help='bitscore threshold (default=50)',
        type=float,
        default="0")
    parser.add_argument(
        '--identity',
        help='identity threshold (default=0)',
        type=float,
        default="0")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args


class GeneObject:
    def __init__(self, gene_id, location, product, color, gene_biotype, note):
        self.gene_id = gene_id
        self.location = location
        self.product = product
        self.color = color
        self.gene_biotype = gene_biotype
        self.note = note


class RepeatObject:
    def __init__(self, repeat_id, location, rpt_family, color, rpt_type, note):
        self.repeat_id = repeat_id
        self.location = location
        self.rpt_family = rpt_family
        self.color = color
        self.rpt_type = rpt_type
        self.note = note


class FeatureObject:
    def __init__(self, feature_id, location, color, note):
        self.feature_id = feature_id
        self.location = location
        self.color = color
        self.note = note


def load_gbks(gbk_list):
    """Load the first instance of GenBank as SerRecord
    Args:
    """
    record_list = []
    for gbk_file in gbk_list:
        logger.info("Loading GenBank...")
        records = SeqIO.parse(gbk_file, 'genbank')
        record = next(records)
        logger.info("             ... done loading GenBank")
        logger.info("Input GenBank       : {}".format(gbk_file))
        record_list.append(record)
        print(record)
    return record_list


def calculate_gc_percent(sequence):
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    gc_percent = 100 * ((g_count + c_count) / len(sequence))
    gc_percent = round(gc_percent, 2)
    return gc_percent


def sliding_window(seq, window, step):
    for start in range(0, len(seq), step):
        end = start + window
        if end > len(seq):
            overhang_length = (end - len(seq))
            # assuming circular sequence
            out_seq = seq[start:len(seq)] + seq[0:overhang_length]
        else:
            out_seq = seq[start:end]
        yield start, out_seq


def calculate_dinucleotide_skew(seq, base1, base2):
    "Return dinucleotide skew in a given sequence"
    base1_count = seq.count(base1)
    base2_count = seq.count(base2)
    skew = ((base1_count - base2_count) / (base1_count + base2_count))
    return skew


def skew_df(record, window, step, nt):
    nt_list = list(nt)
    nt_1 = nt_list[0]
    nt_2 = nt_list[1]
    skew_sum = 0
    skew_dict = {}
    content_dict = {}
    skew_cumulative_dict = {}
    seq = record.seq.upper()
    for start, seq_part in sliding_window(seq, window, step):
        dinucleotide_content = (seq_part.count(
            nt_1) + seq_part.count(nt_1)) / len(seq_part)
        content_dict[start] = dinucleotide_content
        if dinucleotide_content == 0:
            skew = 0
        else:
            skew = calculate_dinucleotide_skew(seq_part, nt_1, nt_2)
        skew_dict[start] = skew
        skew_sum = (skew_sum + skew)
        skew_cumulative_dict[start] = (skew_sum)
    max_skew_abs = abs(max(skew_dict.values(), key=abs))
    max_skew_cumulative_abs = abs(max(skew_cumulative_dict.values(), key=abs))
    factor = (max_skew_abs / max_skew_cumulative_abs)
    skew_cumulative_dict.update((x, (y * factor))
                                for x, y in skew_cumulative_dict.items())
    content_legend = "{} content".format(nt)
    skew_legend = "{} skew".format(nt)
    cumulative_skew_legend = "Cumulative {} skew, normalized".format(nt)
    df = pd.DataFrame({content_legend: pd.Series(content_dict), skew_legend: pd.Series(
        skew_dict), cumulative_skew_legend: pd.Series(skew_cumulative_dict)})
    return df


def load_comparisons(
        comparison_files,
        evalue_threshold,
        bitscore_threshold,
        identity_threshold):
    """
    Load tab-separated comparsion files
    Args:
    """
    logger.info(
        "e-value threshold: {}; bitscore threshold: {}; identity threshold: {}".format(
            evalue_threshold,
            bitscore_threshold,
            identity_threshold))
    comparison_list = []
    for comparison_file in comparison_files:
        logger.info("Loading comparison file...")
        df = pd.read_csv(
            comparison_file,
            sep='\t',
            comment='#',
            names=(
                "query",
                "subject",
                "identity",
                "alignment_length",
                "mismatches",
                "gap_opens",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore"))
        df = df[(df['evalue'] <= evalue_threshold) & (df['bitscore'] >=
                                                      bitscore_threshold) & (df['identity'] >= identity_threshold)]
        logger.info("             ... done loading comparison file")
        logger.info("Input comparison file       : {}".format(comparison_file))
        comparison_list.append(df)
    return comparison_list


def gbk_to_seqrecords(gbk_list):
    """Load the first instance of GenBank as SerRecord
    Args:
    """
    record_list = []
    for gbk_file in gbk_list:
        logger.info("Loading GenBank...")
        records = SeqIO.parse(gbk_file, 'genbank')
        record = next(records)
        logger.info("             ... done loading GenBank")
        logger.info("Input GenBank       : {}".format(gbk_file))
        record_list.append(record)
    return record_list


def get_gene_biotype(feature):
    if feature.type == 'CDS':
        gene_biotype = "protein_coding"
    elif feature.type == 'rRNA':
        gene_biotype = "rRNA"
    elif feature.type == 'tRNA':
        gene_biotype = "tRNA"
    else:
        gene_biotype = "undefined"
    return gene_biotype


def get_color(feature, color_table):
    # Set default color
    gene_biotype = get_gene_biotype(feature)
    if gene_biotype == "protein_coding":
        # color = "#47b8f8"
        color = "lightgray"
    elif gene_biotype == "rRNA":
        color = "#009e73"
    elif gene_biotype == "tRNA":
        color = "#e69f00"
    else:
        color = "gray"
    if 'note' in feature.qualifiers.keys():
        note = feature.qualifiers['note'][0]
    else:
        note = "none"
    # Apply feature-specific color
    target_row = color_table[(color_table['feature_type'] == feature.type) & (color_table['qualifier_key'] == "note") & (
        color_table['value'].str.contains(note)) & (color_table['color'].notna())]
    if (len(target_row) == 1):
        color = target_row['color'].tolist()[0]
    return color


def create_gene_object(feature_id, feature, color_table):
    exon_coordinates = feature.location.parts
    if 'note' in feature.qualifiers.keys():
        note = feature.qualifiers['note'][0]
    else:
        note = ""
    product = feature.qualifiers['product'][0]
    gene_biotype = get_gene_biotype(feature)
    color = get_color(feature, color_table)
    exon_list = get_exon_and_intron_coordinates(exon_coordinates)
    gene_object = GeneObject(
        feature_id,
        exon_list,
        product,
        color,
        gene_biotype,
        note)
    return gene_object


def create_repeat_object(repeat_id, feature):
    coordinates = feature.location.parts
    if feature.type == 'repeat_region':
        gene_biotype = "repeat_region"
        if 'rpt_family' in feature.qualifiers.keys():
            rpt_family = feature.qualifiers['rpt_family'][0]
        else:
            rpt_family = "undefined"
        color = "#d3d3d3"
        if 'rpt_type' in feature.qualifiers.keys():
            rpt_type = feature.qualifiers['rpt_type'][0]
        else:
            rpt_type = "undefined"
        if 'note' in feature.qualifiers.keys():
            note = feature.qualifiers['note'][0]
        else:
            note = ""
        location = get_exon_and_intron_coordinates(coordinates)
        repeat_object = RepeatObject(
            repeat_id, location, rpt_family, color, rpt_type, note)
    else:
        raise ValueError("feature not repeat")
    return repeat_object


def create_feature_object(feature_id, feature):
    coordinates = feature.location.parts
    if feature.type == 'misc_feature':
        if 'note' in feature.qualifiers.keys():
            note = list(feature.qualifiers['note'])
        else:
            note = ""
        color = "#d3d3d3"
        location = get_exon_and_intron_coordinates(coordinates)

        feature_object = FeatureObject(feature_id, location, note, color, note)
    else:
        raise ValueError("feature not misc_feature")
    return feature_object


def create_feature_dict(gb_record, color_table):
    feature_dict = {}
    locus_count = 0
    repeat_count = 0
    feature_count = 0
    for feature in gb_record.features:
        if (feature.type == 'CDS') or (
                feature.type == 'rRNA') or (feature.type == 'tRNA'):
            locus_count = locus_count + 1
            locus_id = "gene_" + str(locus_count).zfill(3)
            gene_object = create_gene_object(locus_id, feature, color_table)
            feature_dict[locus_id] = gene_object
        elif feature.type == 'repeat_region':
            repeat_count = repeat_count + 1
            repeat_id = "crt_" + str(repeat_count).zfill(3)
            repeat_object = create_repeat_object(repeat_id, feature)
            feature_dict[repeat_id] = repeat_object
        elif feature.type == 'misc_feature':
            feature_count = feature_count + 1
            feature_id = "feature_" + str(feature_count).zfill(3)
            feature_object = create_feature_object(feature_id, feature)
            feature_dict[feature_id] = feature_object
    return feature_dict


def get_product_color(notes):
    """
    Not implemented yet.
    """
    product_color = "lightgray"
    return product_color


def get_strand(strand_value: int) -> str:
    """
    Convert strand value to string.
    Args:
        strand_value (int): strand value (1: positive, -1: negative).
    """
    if strand_value == 1:
        strand = "positive"
    elif strand_value == -1:
        strand = "negative"
    else:
        strand = "undefined"
    return strand


def get_exon_coordinate(exon_line, exon_count, last_or_not):
    """
    Return a given exon's id, strand, starting coordinate, ending coordinate, and boolean whether it is the final exon or not (which is important to decide if the exon should be drawn as an arrow or a rectangle)
    Args:
        exon_line (list):
    """
    exon_strand = get_strand(exon_line.strand)
    exon_count = exon_count + 1
    exon_id = str(exon_count).zfill(3)
    exon_start = int(exon_line.start)
    exon_end = int(exon_line.end)
    exon_coordinate = [
        "exon",
        exon_id,
        exon_strand,
        exon_start,
        exon_end,
        last_or_not]
    return exon_count, exon_coordinate


def get_intron_coordinate(previous_exon, current_exon, intron_count):
    intron_count = intron_count + 1
    intron_id = str(intron_count).zfill(3)
    intron_strand = previous_exon[2]
    previous_exon_start = previous_exon[3]
    previous_exon_end = previous_exon[4]
    current_exon_start = current_exon[3]
    current_exon_end = current_exon[4]
    if intron_strand == "positive":
        intron_start = int(previous_exon_end) + 1
        intron_end = int(current_exon_start) - 1
    elif intron_strand == "negative":
        intron_start = int(current_exon_end) + 1
        intron_end = int(previous_exon_start) - 1
    intron_coordinate = [
        "intron",
        intron_id,
        intron_strand,
        intron_start,
        intron_end,
        False]
    return intron_count, intron_coordinate


def get_exon_and_intron_coordinates(exons: list) -> list:
    exon_list = []
    exon_count = 0
    intron_start = ''
    intron_end = ''
    intron_count = 0
    num_exons = len(exons)
    intron_strand = ''
    if num_exons == 1:
        exon_count, exon_coord = get_exon_coordinate(
            exons[0], exon_count, True)
        exon_list.append(exon_coord)
    elif num_exons != 1:
        for exon in exons:
            if exon_count == 0:
                exon_count, exon_coord = get_exon_coordinate(
                    exon, exon_count, False)
                exon_list.append(exon_coord)
            elif exon_count != 0:
                if exon_count < (num_exons - 1):
                    exon_count, exon_coord = get_exon_coordinate(
                        exon, exon_count, False)
                    if exon_list[-1]:
                        intron_count, intron_coord = get_intron_coordinate(
                            exon_list[-1], exon_coord, intron_count)
                        exon_list.append(intron_coord)
                    else:
                        pass
                    exon_list.append(exon_coord)
                elif exon_count == (num_exons - 1):
                    exon_count, exon_coord = get_exon_coordinate(
                        exon, exon_count, True)
                    intron_count, intron_coord = get_intron_coordinate(
                        exon_list[-1], exon_coord, intron_count)
                    exon_list.append(intron_coord)
                    exon_list.append(exon_coord)
    return exon_list


def set_arrow_shoulder(feat_strand, arrow_end, cds_arrow_length) -> int:
    if feat_strand == "positive":
        shoulder = int(arrow_end - cds_arrow_length)
    else:
        shoulder = int(arrow_end + cds_arrow_length)
    return shoulder


def get_coordinate(coord):
    feat_type = coord[0]
    feat_strand = coord[2]
    feat_start = coord[3]
    feat_end = coord[4]
    return feat_type, feat_strand, feat_start, feat_end


def get_gc_path(fig_horizontal_length, df, total_len,
                genome_size_normalization_factor, track_height):
    coodinates_list = []
    skew_start_x = 0
    skew_start_y = 0
    skew_start_position = "M{} {}".format(skew_start_x, skew_start_y)
    coodinates_list.append(skew_start_position)
    column = 'GC content'
    mean = float(df[column].mean())
    max_diff = float((df[column] - mean).abs().max())
    for index, row in df.iterrows():
        value = float(row[column])
        diff = (value - mean)
        height_of_coordinate = - (0.5 * track_height * (diff / max_diff))
        x_corrdinate = normalize_position_to_track(
            index, total_len, fig_horizontal_length, genome_size_normalization_factor)
        y_corrdinate = height_of_coordinate
        corrdinate = "L{} {}".format(str(x_corrdinate), str(y_corrdinate))
        coodinates_list.append(corrdinate)
    end_coordinate = "L{} {}".format(str(skew_start_x), str(skew_start_y))
    coodinates_list.append(end_coordinate)
    gc_desc = "{}".format(''.join(coodinates_list))
    gc_desc += "z"
    return gc_desc


def gc_content(gb_record, fig_horizontal_length,
               longest_genome, df, track_height):
    total_len = len(gb_record.seq)
    genome_size_normalization_factor = total_len / longest_genome
    gc_group = Group(id="gc_content")
    gc_path_desc = get_gc_path(
        fig_horizontal_length,
        df,
        total_len,
        genome_size_normalization_factor,
        track_height)
    gc_path = Path(
        d=gc_path_desc,
        fill="#808080",
        stroke='none',
        fill_opacity=1,
        fill_rule="evenodd")
    gc_group.add(gc_path)
    return gc_group


def normalize_position_to_track(
        position,
        genome_length,
        fig_horizontal_length,
        genome_size_normalization_factor):
    normalized_position = fig_horizontal_length * \
        (position / genome_length) * genome_size_normalization_factor
    return normalized_position


def compute_strand(strandedness, strand):
    if strandedness:
        offset = 0.1
        height = 0.45
        factors_positive = [-(height + offset), -
                            ((height + 2 * offset) / 2), -offset]
        factors_negative = [
            offset,
            ((height + 2 * offset) / 2),
            (height + offset)]
        if strand == "positive":
            factors = factors_positive
        else:
            factors = factors_negative
    else:
        factors = [-0.5, 0, 0.5]
    return factors


def create_intron_path_linear(
        cds_height,
        coord_dict,
        genome_length,
        fig_horizontal_length,
        genome_size_normalization_factor,
        strandedness):
    coord_path = []
    coord_path.append("intron")
    feat_start = coord_dict['feat_start']
    feat_end = coord_dict['feat_end']
    feat_strand = coord_dict["feat_strand"]
    factors = compute_strand(strandedness, feat_strand)
    normalized_start = normalize_position_to_track(
        feat_start,
        genome_length,
        fig_horizontal_length,
        genome_size_normalization_factor)
    normalized_end = normalize_position_to_track(
        feat_end,
        genome_length,
        fig_horizontal_length,
        genome_size_normalization_factor)
    start_x = normalized_start
    start_y = cds_height * factors[1]
    end_x = normalized_end
    end_y = cds_height * factors[1]
    feature_path = "M " + str(start_x) + "," + str(start_y) + \
        "L" + str(end_x) + "," + str(end_y) + " z"
    coord_path.append(feature_path)
    return coord_path


def set_arrow_shoulder_linear(
        feat_strand, normalized_end, normalized_arrow_length):
    if feat_strand == "positive":
        shoulder = int(normalized_end - normalized_arrow_length)
    else:
        shoulder = int(normalized_end + normalized_arrow_length)
    return shoulder


def create_arrowhead_path_linear(
        coord_dict,
        genome_length,
        cds_height,
        fig_horizontal_length,
        cds_arrow_length,
        genome_size_normalization_factor,
        strandedness):
    coord_path = []
    coord_path.append("exon")
    feat_start = int(coord_dict["feat_start"])
    feat_end = int(coord_dict["feat_end"])
    feat_strand = coord_dict["feat_strand"]
    factors = compute_strand(strandedness, feat_strand)
    normalized_start = normalize_position_to_track(
        feat_start,
        genome_length,
        fig_horizontal_length,
        genome_size_normalization_factor)
    normalized_end = normalize_position_to_track(
        feat_end,
        genome_length,
        fig_horizontal_length,
        genome_size_normalization_factor)
    normalized_feat_len = normalized_end - normalized_start
    normalized_arrow_length = fig_horizontal_length * \
        (cds_arrow_length / genome_length) * genome_size_normalization_factor
    arrow_strand_dict = {
        "positive": [
            normalized_start,
            normalized_end],
        "negative": [
            normalized_end,
            normalized_start]}
    arrow_start, arrow_end = arrow_strand_dict[feat_strand]
    shoulder = set_arrow_shoulder_linear(
        feat_strand, arrow_end, normalized_arrow_length)
    if abs(normalized_feat_len) < normalized_arrow_length:
        point_x = arrow_end
        point_y = cds_height * factors[1]
        start_x_1 = arrow_start
        start_y_1 = cds_height * factors[0]
        start_x_2 = arrow_start
        start_y_2 = cds_height * factors[2]
        feature_path = "M " + str(start_x_1) + "," + str(start_y_1) + " L " + str(
            point_x) + "," + str(point_y) + " L " + str(start_x_2) + "," + str(start_y_2) + " z"
    else:
        point_x = arrow_end
        point_y = cds_height * factors[1]
        start_x_1 = arrow_start
        start_y_1 = cds_height * factors[0]
        start_x_2 = shoulder
        start_y_2 = cds_height * factors[2]
        end_x_1 = shoulder
        end_y_1 = cds_height * factors[0]
        end_x_2 = arrow_start
        end_y_2 = cds_height * factors[2]
        feature_path = "M " + str(start_x_1) + "," + str(start_y_1) + "L " + str(end_x_1) + "," + str(end_y_1) + " L " + str(
            point_x) + "," + str(point_y) + "L " + str(start_x_2) + "," + str(start_y_2) + " L " + str(end_x_2) + "," + str(end_y_2) + " z"
    coord_path.append(feature_path)
    return coord_path


def create_rectangle_path_linear(
        coord_dict,
        genome_length,
        cds_height,
        fig_horizontal_length,
        genome_size_normalization_factor,
        strandedness):
    coord_path = []
    coord_path.append("exon")
    feat_start = coord_dict['feat_start']
    feat_end = coord_dict['feat_end']
    feat_strand = coord_dict["feat_strand"]
    factors = compute_strand(strandedness, feat_strand)
    normalized_start = fig_horizontal_length * \
        (coord_dict['feat_start'] / genome_length) * genome_size_normalization_factor
    normalized_end = fig_horizontal_length * \
        (coord_dict['feat_end'] / genome_length) * genome_size_normalization_factor
    start_x_1 = normalized_start
    start_y_1 = cds_height * factors[0]
    start_x_2 = normalized_start
    start_y_2 = cds_height * factors[2]
    end_x_1 = normalized_end
    end_y_1 = cds_height * factors[0]
    end_x_2 = normalized_end
    end_y_2 = cds_height * factors[2]
    feature_path = "M " + str(start_x_1) + "," + str(start_y_1) + "L" + str(end_x_1) + "," + str(
        end_y_1) + " L" + str(end_x_2) + "," + str(end_y_2) + "L" + str(start_x_2) + "," + str(start_y_2) + " z"
    coord_path.append(feature_path)
    return coord_path


def create_canvas(
        file_name,
        width,
        canvas_padding,
        num_of_entries,
        cds_height,
        comparison_height,
        vertical_padding,
        vertical_offset):
    total_width = str(width + 2 * canvas_padding)
    total_height = str(2 *
                       vertical_offset +
                       (2 *
                        cds_height +
                        2 *
                        vertical_padding) +
                       (cds_height +
                        comparison_height +
                        2 *
                        vertical_padding) *
                       (num_of_entries -
                           1))
    dwg = svgwrite.Drawing(
        filename=file_name +
        ".svg",
        size=(
            total_width +
            'px',
            total_height +
            'px'),
        viewBox=(
            '0 0 ' +
            str(total_width) +
            ' ' +
            str(total_height)))
    return dwg


def axis_path_linear(
        genome_length,
        fig_horizontal_length,
        genome_size_normalization_factor):
    bar_length = fig_horizontal_length * genome_size_normalization_factor
    start_x = 0
    start_y = 0
    end_x = bar_length
    end_y = 0
    axis_path = Line(
        start=(
            start_x,
            start_y),
        end=(
            end_x,
            end_y),
        stroke='lightgray',
        stroke_width=2,
        fill='none')
    return axis_path


def gene_to_path_linear(
        gene_object,
        genome_length: int,
        cds_height: int,
        fig_horizontal_length: int,
        genome_size_normalization_factor,
        strandedness,
        cds_arrow_length):
    coords = gene_object.location
    gene_body_color = gene_object.color
    coordinates_paths = []
    for coord in coords:
        coord_path = []
        coord_dict = {
            'feat_type': coord[0],
            'feat_strand': coord[2],
            'feat_start': coord[3],
            'feat_end': coord[4]}
        feat_type = coord_dict['feat_type']
        if feat_type == "intron":
            coord_path = create_intron_path_linear(
                cds_height,
                coord_dict,
                genome_length,
                fig_horizontal_length,
                genome_size_normalization_factor,
                strandedness)
        elif feat_type == "exon":
            if coord[5]:
                coord_path = create_arrowhead_path_linear(
                    coord_dict,
                    genome_length,
                    cds_height,
                    fig_horizontal_length,
                    cds_arrow_length,
                    genome_size_normalization_factor,
                    strandedness)
            elif coord[5] == False:
                coord_path = create_rectangle_path_linear(
                    coord_dict,
                    genome_length,
                    cds_height,
                    fig_horizontal_length,
                    genome_size_normalization_factor,
                    strandedness)
        else:
            continue
        coordinates_paths.append(coord_path)
    return coordinates_paths


def draw_paths_linear(
        feature_dict,
        cds_height,
        genome_length,
        group: Group,
        fig_horizontal_length: int,
        fig_vertical_length: int,
        genome_size_normalization_factor,
        strandedness,
        cds_arrow_length) -> svgwrite.container.Group:
    total_len = genome_length
    record_group = Group(id="record")
    axis_path = axis_path_linear(
        genome_length,
        fig_horizontal_length,
        genome_size_normalization_factor)
    group.add(axis_path)
    for key in feature_dict:
        feature_object = feature_dict[key]
        if isinstance(feature_object, GeneObject):
            gene_object = feature_object
            gene_paths = gene_to_path_linear(
                gene_object,
                genome_length,
                cds_height,
                fig_horizontal_length,
                genome_size_normalization_factor,
                strandedness, cds_arrow_length)
            for gene_path in gene_paths:
                feat_type = gene_path[0]
                if feat_type == "exon":
                    exon_path = Path(
                        d=gene_path[1],
                        fill=gene_object.color,
                        stroke='black',
                        stroke_width=0.5)
                    group.add(exon_path)
                elif feat_type == "intron":
                    intron_path = Path(
                        d=gene_path[1], stroke='gray', stroke_width=0.5)
                    group.add(intron_path)
        elif isinstance(feature_object, RepeatObject):
            coords = feature_object.location
            color = feature_object.color
            coordinates_paths = []
            repeat_paths = []
            for coord in coords:
                coord_path = []
                coord_dict = {
                    'feat_type': coord[0],
                    'feat_strand': coord[2],
                    'feat_start': coord[3],
                    'feat_end': coord[4]}
                repeat_path = create_rectangle_path_linear(
                    coord_dict,
                    genome_length,
                    cds_height,
                    fig_horizontal_length,
                    genome_size_normalization_factor,
                    strandedness)
                unit_path = Path(
                    d=repeat_path[1],
                    fill=color,
                    stroke='none',
                    stroke_width=0.5)
                group.add(unit_path)
        elif isinstance(feature_object, FeatureObject):
            coords = feature_object.location
            color = feature_object.color
            coordinates_paths = []
            feature_paths = []
            for coord in coords:
                coord_path = []
                coord_dict = {
                    'feat_type': coord[0],
                    'feat_strand': coord[2],
                    'feat_start': coord[3],
                    'feat_end': coord[4]}
                feature_path = create_rectangle_path_linear(
                    coord_dict,
                    genome_length,
                    cds_height,
                    fig_horizontal_length,
                    genome_size_normalization_factor,
                    strandedness)
                feature_path = Path(
                    d=feature_path[1],
                    fill=color,
                    stroke='none',
                    stroke_width=0.5)
                group.add(feature_path)
    return group


def normalize(position, longest_genome, fig_horizontal_length):
    return fig_horizontal_length * (position / longest_genome)


def add_match_path(
        row,
        mutual_hit_group,
        longest_genome,
        fig_horizontal_length,
        comparison_height,
        record_dict,
        align_center):
    query_id = row.query
    subject_id = row.subject
    if align_center:
        query_offset_x = (longest_genome - record_dict[query_id]) / 2
        subject_offset_x = (longest_genome - record_dict[subject_id]) / 2
    else:
        query_offset_x = 0
        subject_offset_x = 0
    query_start = row.qstart + query_offset_x
    query_end = row.qend + query_offset_x
    subject_start = row.sstart + subject_offset_x
    subject_end = row.send + subject_offset_x
    percent_ident = row.identity
    query_start_x = normalize(
        query_start,
        longest_genome,
        fig_horizontal_length)
    query_start_y = 0
    query_end_x = normalize(query_end, longest_genome, fig_horizontal_length)
    query_end_y = 0
    subject_start_x = normalize(
        subject_start,
        longest_genome,
        fig_horizontal_length)
    subject_start_y = comparison_height
    subject_end_x = normalize(
        subject_end,
        longest_genome,
        fig_horizontal_length)
    subject_end_y = comparison_height
    '''
    if percent_ident < 60:
        opacity = 0.2
    elif 60 <= percent_ident < 70:
        opacity = 0.3
    elif 70 <= percent_ident < 80:
        opacity = 0.4
    elif 80 <= percent_ident < 90:
        opacity = 0.5
    elif 90 <= percent_ident <= 100:
        opacity = 0.6
    '''
    opacity = 0.7
    feature_path_desc = "M " + str(query_start_x) + "," + str(query_start_y) + "L" + str(query_end_x) + "," + str(
        query_end_y) + " L" + str(subject_end_x) + "," + str(subject_end_y) + "L" + str(subject_start_x) + "," + str(subject_start_y) + " z"
    feature_path = Path(
        d=feature_path_desc,
        fill='#92daff',
        fill_opacity=opacity,
        stroke='none',
        stroke_width=0.2)
    mutual_hit_group.add(feature_path)
    # fill='#ff99ff'
    return mutual_hit_group


def length_bar_to_track(fig_width: int, cds_height: int, longest_genome: int):
    if longest_genome < 2000:
        tick = 100
        label_text = str(tick) + " bp"
    elif 2000 <= longest_genome < 20000:
        tick = 1000
        label_text = str(int(tick / 1000)) + " kbp"
    elif 20000 <= longest_genome < 50000:
        tick = 5000
        label_text = str(int(tick / 1000)) + " kbp"
    elif 50000 <= longest_genome < 150000:
        tick = 10000
        label_text = str(int(tick / 1000)) + " kbp"
    elif 150000 <= longest_genome < 250000:
        tick = 50000
        label_text = str(int(tick / 1000)) + " kbp"
    elif 250000 <= longest_genome < 1000000:
        tick = 100000
        label_text = str(int(tick / 1000)) + " kbp"
    elif 1000000 <= longest_genome < 2000000:
        tick = 200000
        label_text = str(int(tick / 1000)) + " kbp"
    elif 2000000 <= longest_genome:
        tick = 500000
        label_text = str(int(tick / 1000)) + " kbp"
    record_group = Group(id="length_bar")
    bar_length = fig_width * (tick / longest_genome)
    start_x = (0.98 * fig_width - bar_length + 1)
    start_y = 0
    end_x = 0.98 * fig_width
    end_y = 0
    length_bar = Line(
        start=(
            start_x,
            start_y),
        end=(
            end_x,
            end_y),
        stroke='black',
        stroke_width=3,
        fill='none')
    label_path = Text(
        label_text,
        insert=(
            start_x - 10,
            start_y),
        stroke='none',
        fill='black',
        font_size='10pt',
        font_weight="normal",
        font_family="Arial",
        text_anchor="end",
        dominant_baseline="middle")
    record_group.add(length_bar)
    record_group.add(label_path)
    return record_group


def record_to_track(record, fig_width, cds_height,
                    longest_genome, strandedness, color_table):
    track_id = str(record.annotations["accessions"][0])
    record_group = Group(id=track_id)
    record_name = record.annotations["source"]
    recod_length = len(record.seq)
    genome_size_normalization_factor = recod_length / longest_genome
    cds_arrow_length = 0.0012 * longest_genome
    track_height = cds_height
    feature_dict = create_feature_dict(record, color_table)
    record_group = draw_paths_linear(
        feature_dict,
        cds_height,
        recod_length,
        record_group,
        fig_width,
        track_height,
        genome_size_normalization_factor,
        strandedness, cds_arrow_length)
    return record_group


def record_definition_to_track(record):
    cds_list = []
    strain_name = "unknown"
    for feature in record.features:
        if feature.type == "source":
            if 'isolate' in feature.qualifiers.keys():
                strain_name = feature.qualifiers['isolate'][0]
            elif 'strain' in feature.qualifiers.keys():
                strain_name = feature.qualifiers['strain'][0]
            else:
                strain_name = "unknown"
        if feature.type == "CDS":
            cds_list.append(feature)
    track_id = str(strain_name)
    record_group = Group(id=track_id)
    record_name = strain_name
    recoed_length = len(record.seq)
    num_cds = len(cds_list)
    length_label = "{:,} bp".format(recoed_length)
    cds_label = "{:,} CDS".format(num_cds)
    start_x = 0
    start_y = 0
    name_path = Text(
        record_name,
        insert=(
            start_x,
            start_y - 17),
        stroke='none',
        fill='black',
        font_size='10pt',
        font_weight="normal",
        font_family="Arial",
        text_anchor="middle",
        dominant_baseline="middle")
    length_path = Text(
        length_label,
        insert=(
            start_x,
            start_y),
        stroke='none',
        fill='black',
        font_size='10pt',
        font_weight="normal",
        font_family="Arial",
        text_anchor="middle",
        dominant_baseline="middle")
    num_cds_path = Text(
        cds_label,
        insert=(
            start_x,
            start_y + 17),
        stroke='none',
        fill='black',
        font_size='10pt',
        font_weight="normal",
        font_family="Arial",
        text_anchor="middle",
        dominant_baseline="middle")
    record_group.add(name_path)
    record_group.add(length_path)
    record_group.add(num_cds_path)
    return record_group


def add_records_on_canvas(
        canvas,
        records,
        longest_genome,
        cds_height,
        alignment_width,
        horizontal_offset,
        vertical_padding,
        comparison_height,
        strandedness,
        vertical_offset,
        show_gc,
        align_center,
        color_table):
    count = 0
    for record in records:
        count += 1
        record_definition_group = record_definition_to_track(record)
        record_group = record_to_track(
            record,
            alignment_width,
            cds_height,
            longest_genome,
            strandedness, color_table)
        offset = vertical_offset + \
            (cds_height + comparison_height + 2 * vertical_padding) * (count - 1)
        if align_center:
            offset_x = alignment_width * \
                ((longest_genome - len(record.seq)) / longest_genome) / 2
        else:
            offset_x = 0
        record_group.translate(offset_x + horizontal_offset, offset)
        canvas.add(record_group)
        if show_gc:
            df = skew_df(record, 1000, 100, "GC")
            if strandedness:
                gc_height = cds_height
            else:
                gc_height = 2 * cds_height
            gc_content_group = gc_content(
                record, alignment_width, longest_genome, df, gc_height)
            gc_content_group.translate(
                offset_x + horizontal_offset,
                offset + gc_height)
            canvas.add(gc_content_group)
        record_definition_group.translate(
            offset_x + horizontal_offset * 0.5, offset)
        canvas.add(record_definition_group)
    return canvas


def comparison_to_track(
        comparison_df,
        comparison_count,
        longest_genome,
        alignment_width,
        comparison_height,
        record_dict,
        align_center):
    track_id = "comparison" + str(comparison_count)
    comparison_group = Group(id=track_id)
    for row in comparison_df.itertuples():
        comparison_group = add_match_path(
            row,
            comparison_group,
            longest_genome,
            alignment_width,
            comparison_height,
            record_dict,
            align_center)
    return comparison_group


def add_comparison_on_canvas(
        canvas,
        comparisons,
        longest_genome,
        cds_height,
        alignment_width,
        horizontal_offset,
        vertical_padding,
        comparison_height,
        vertical_offset,
        record_dict,
        align_center):
    comparison_count = 0
    for comparison in comparisons:
        comparison_count += 1
        comparison_group = comparison_to_track(
            comparison,
            comparison_count,
            longest_genome,
            alignment_width,
            comparison_height, record_dict, align_center)
        offset = vertical_offset + (0.5 * cds_height + vertical_padding) + (
            (cds_height + comparison_height + 2 * vertical_padding) * (comparison_count - 1))
        comparison_group.translate(horizontal_offset, offset)
        canvas.add(comparison_group)
    return canvas


def main():
    args = _get_args()
    comps = []
    genbank_files = args.input
    color_table = args.table
    color_table = pd.read_csv(
        color_table,
        sep='\t',
        names=(
            'feature_type',
            'qualifier_key',
            'value',
            'color'))
    blast_files = args.blast
    records = load_gbks(genbank_files)
    out_file_prefix = args.output
    fig_width = 2000
    strandedness = args.separate_strands
    show_gc = args.show_gc
    align_center = args.align_center
    evalue = args.evalue
    bitscore = args.bitscore
    identity = args.identity
    if strandedness:
        canvas_padding = 20
        cds_height = 20
        comparison_height = 50
    else:
        canvas_padding = 20
        cds_height = 10
        comparison_height = 50
    comparison_height = 75
    vertical_padding = 5
    record_dict = {}
    for record in records:
        if record.id == ".":
            record.id = record.name
        print(record.id)
        record_dict[record.id] = len(record.seq)
    longest_genome = max([len(record.seq) for record in records])
    num_of_entries = len(records)
    vertical_offset = 40
    canvas = create_canvas(
        out_file_prefix,
        fig_width,
        canvas_padding,
        num_of_entries,
        cds_height,
        comparison_height,
        vertical_padding,
        vertical_offset)
    horizontal_offset = 80

    alignment_width = fig_width - horizontal_offset
    canvas = add_records_on_canvas(
        canvas,
        records,
        longest_genome,
        cds_height,
        alignment_width,
        horizontal_offset,
        vertical_padding,
        comparison_height,
        strandedness,
        vertical_offset,
        show_gc, align_center, color_table)
    if blast_files is not None:
        comparisons = load_comparisons(blast_files, evalue, bitscore, identity)
        canvas = add_comparison_on_canvas(
            canvas,
            comparisons,
            longest_genome,
            cds_height,
            alignment_width,
            horizontal_offset,
            vertical_padding,
            comparison_height,
            vertical_offset, record_dict, align_center)
    if (show_gc) and not (strandedness):
        add_margin = 2 * cds_height
    else:
        add_margin = 0
    length_bar_group = length_bar_to_track(
        alignment_width, cds_height, longest_genome)
    offset_for_length_bar = (add_margin +
                             vertical_offset +
                             2 *
                             (cds_height +
                              vertical_padding) +
                             (cds_height +
                              comparison_height +
                              2 *
                              vertical_padding) *
                             (len(records) -
                                 1))
    length_bar_group.translate(horizontal_offset, offset_for_length_bar)
    canvas.add(length_bar_group)
    canvas.filename = "{}.svg".format(out_file_prefix)
    canvas.save()


if __name__ == "__main__":
    main()
