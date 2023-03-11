#!/usr/bin/env python
# coding: utf-8

import sys
import os
import re
import argparse
import pandas as pd
import math
from math import pi
import svgwrite
from svgwrite.shapes import Circle
from svgwrite.path import Path
from svgwrite.text import Text
from svgwrite.masking import ClipPath, Mask
from svgwrite.container import Group
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation


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


def _get_args():
    parser = argparse.ArgumentParser(
        description='Generate genome diagrams in SVG. Diagrams for multiple entries are saved separately (hence the lack of output file name option).'
    parser.add_argument(
        '-i',
        '--input',
        help='Genbank/DDBJ flatfile (required)',
        type=str,
        required=True)
    parser.add_argument(
        '-t',
        '--table',
        help='color table (optional)',
        type=str,
        default="")
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
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    #args = parser.parse_args(args=['-i', 'CN01.gb', '-t', 'color_table.txt'])
    return args


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
        skew = calculate_dinucleotide_skew(seq_part, nt_1, nt_2)
        dinucleotide_content = (seq_part.count(
            nt_1) + seq_part.count(nt_1)) / len(seq_part)
        content_dict[start] = dinucleotide_content
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


def circle_path_circular(radius, total_len, norm_factor):
    norm_radius = radius * norm_factor
    circle_radius = norm_radius
    circle_start_x = circle_radius * math.cos(math.radians(360.0 * 0 - 90))
    circle_start_y = circle_radius * math.sin(math.radians(360.0 * 0 - 90))
    circle_desc = ''
    start_coordinate = "M {} {}".format(circle_start_x, circle_start_y)
    circle_desc += start_coordinate
    n_list = []
    for n in range(0, 10, 1):
        deg = (n + 1) * 0.1
        n_x = circle_radius * math.cos(math.radians(360.0 * deg - 90))
        n_y = circle_radius * math.sin(math.radians(360.0 * deg - 90))
        n_list.append(
            "A {} {} 0 0 1 {} {}".format(
                circle_radius,
                circle_radius,
                n_x,
                n_y))
    circle_desc += ''.join(n_list)
    circle_desc += "z"
    return circle_desc


def skew_path_circular(radius, df, total_len, norm_factor, track_width):
    norm_radius = radius * norm_factor
    skew_desc_list = []
    skew_start_x = norm_radius * \
        math.cos(math.radians(360.0 * (0 / total_len) - 90))
    skew_start_y = norm_radius * \
        math.sin(math.radians(360.0 * (0 / total_len) - 90))
    skew_start_position = "M{} {}".format(skew_start_x, skew_start_y)
    skew_desc_list.append(skew_start_position)
    column = 'GC skew'
    mean = float(df[column].mean())
    max_diff = float((df[column] - mean).abs().max())
    for index, row in df.iterrows():
        value = float(row[column])
        diff = (value - mean)
        radius_of_coordinate = (
            norm_radius + (0.5 * track_width * (diff / max_diff)))
        x_corrdinate = radius_of_coordinate * \
            math.cos(math.radians(360.0 * (index / total_len) - 90))
        y_corrdinate = radius_of_coordinate * \
            math.sin(math.radians(360.0 * (index / total_len) - 90))
        corrdinate = "L{} {}".format(str(x_corrdinate), str(y_corrdinate))
        skew_desc_list.append(corrdinate)
    skew_desc = "{}".format(''.join(skew_desc_list))
    skew_desc += "z"
    circle_desc = circle_path_circular(radius, total_len, norm_factor)
    skew_desc += circle_desc
    return skew_desc


def gc_path_circular(radius, df, total_len, norm_factor, track_width):
    norm_radius = radius * norm_factor
    coodinates_list = []
    skew_start_x = norm_radius * \
        math.cos(math.radians(360.0 * (0 / total_len) - 90))
    skew_start_y = norm_radius * \
        math.sin(math.radians(360.0 * (0 / total_len) - 90))
    skew_start_position = "M{} {}".format(skew_start_x, skew_start_y)
    coodinates_list.append(skew_start_position)
    column = 'GC content'
    mean = float(df[column].mean())
    max_diff = float((df[column] - mean).abs().max())
    for index, row in df.iterrows():
        value = float(row[column])
        diff = (value - mean)
        radius_of_coordinate = (
            norm_radius + (0.5 * track_width * (diff / max_diff)))
        x_corrdinate = radius_of_coordinate * \
            math.cos(math.radians(360.0 * (index / total_len) - 90))
        y_corrdinate = radius_of_coordinate * \
            math.sin(math.radians(360.0 * (index / total_len) - 90))
        corrdinate = "L{} {}".format(str(x_corrdinate), str(y_corrdinate))
        coodinates_list.append(corrdinate)
    gc_desc = "{}".format(''.join(coodinates_list))
    gc_desc += "z"
    circle_desc = circle_path_circular(radius, total_len, norm_factor)
    gc_desc += circle_desc
    return gc_desc


def get_strand(strand_value):
    if strand_value == 1:
        strand = "positive"
    elif strand_value == -1:
        strand = "negative"
    else:
        strand = "undefined"
    return strand


def get_exon_coordinate(exon_line, exon_count, last_or_not):
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


def get_exon_and_intron_coordinates(exons):
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


def set_arrow_shoulder(feat_strand, arrow_end, cds_arrow_length):
    if feat_strand == "positive":
        shoulder = int(arrow_end - cds_arrow_length)
    else:
        shoulder = int(arrow_end + cds_arrow_length)
    return shoulder


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
        color = "#47b8f8"
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
    if isinstance(color_table, pd.DataFrame):
        target_row = color_table[(color_table['feature_type'] == feature.type) & (color_table['qualifier_key'] == "note") & (
            color_table['value'].str.contains(note,case=False,na=False, regex=False)) & (color_table['color'].notna())]
        if (len(target_row) == 1):
            color = target_row['color'].tolist()[0]
    else:
        pass
    return color


def create_gene_object(feature_id, feature, color_table):
    exon_coordinates = feature.location.parts
    if 'note' in feature.qualifiers.keys():
        note = feature.qualifiers['note'][0]
    else:
        note = "none"
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
        rpt_type = feature.qualifiers['rpt_type'][0]
        if 'note' in feature.qualifiers.keys():
            note = feature.qualifiers['note'][0]
        else:
            note = "none"
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
            note = feature.qualifiers['note'][0]
        else:
            note = "none"
        color = "#d3d3d3"
        location = get_exon_and_intron_coordinates(coordinates)
        note = list(feature.qualifiers['note'])
        feature_object = FeatureObject(feature_id, location, color, note)
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


def get_coordinate(coord):
    feat_type = coord[0]
    feat_strand = coord[2]
    feat_start = coord[3]
    feat_end = coord[4]
    return feat_type, feat_strand, feat_start, feat_end


def compute_factors(strand, track_ratio):
    offset = 0.005
    cds_ratio = track_ratio * 0.25
    base = 1.0
    factors_positive = [base, base + cds_ratio * 0.5, base + cds_ratio]
    factors_negatie = [base - cds_ratio, base - cds_ratio * 0.5, base]
    if strand == "positive":
        factors = [x + offset for x in factors_positive]
    else:
        factors = [x - offset for x in factors_negatie]
    return factors


def create_intron_path_circular(radius, coord_dict, total_length, track_ratio):
    coord_path = []
    coord_path.append("intron")
    feat_strand = coord_dict['feat_strand']
    feat_start = coord_dict['feat_start']
    feat_end = coord_dict['feat_end']
    strand_dict = {"positive": " 0 0 0 ", "negative": " 0 0 1 "}
    param = strand_dict[feat_strand]
    factors = compute_factors(feat_strand, track_ratio)
    start_x_1 = (radius * factors[1]) * math.cos(
        math.radians(360.0 * ((feat_start) / total_length) - 90))
    start_y_1 = (radius * factors[1]) * math.sin(
        math.radians(360.0 * ((feat_start) / total_length) - 90))
    end_x_1 = (radius * factors[1]) * \
        math.cos(math.radians(360.0 * ((feat_end) / total_length) - 90))
    end_y_1 = (radius * factors[1]) * \
        math.sin(math.radians(360.0 * ((feat_end) / total_length) - 90))
    feature_path = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(
        radius) + "," + str(radius) + param + str(end_x_1) + "," + str(end_y_1)  # +" z"
    coord_path.append(feature_path)
    return coord_path


def create_arrowhead_path_circular(
        radius, coord_dict, total_length, cds_arrow_length, track_ratio):
    coord_path = []
    coord_path.append("exon")
    feat_strand = coord_dict['feat_strand']
    factors = compute_factors(feat_strand, track_ratio)
    feat_start = int(coord_dict["feat_start"])
    feat_end = int(coord_dict["feat_end"])
    feat_len = abs(coord_dict["feat_end"] - coord_dict["feat_start"])
    feat_strand = coord_dict["feat_strand"]
    arrow_strand_dict = {
        "positive": [
            feat_start,
            feat_end,
            " 0 0 1 ",
            " 0 0 0 "],
        "negative": [
            feat_end,
            feat_start,
            " 0 0 0 ",
            " 0 0 1 "]}
    arrow_start, arrow_end, param_1, param_2 = arrow_strand_dict[feat_strand]
    shoulder = set_arrow_shoulder(feat_strand, arrow_end, cds_arrow_length)
    if abs(feat_len) < cds_arrow_length:
        point_x = (
            radius * factors[1]) * math.cos(math.radians(360.0 * (arrow_end / total_length) - 90))
        point_y = (
            radius * factors[1]) * math.sin(math.radians(360.0 * (arrow_end / total_length) - 90))
        start_x_1 = (radius * factors[0]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_y_1 = (radius * factors[0]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_x_2 = (radius * factors[2]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_y_2 = (radius * factors[2]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        feature_path = "M " + str(start_x_1) + "," + str(start_y_1) + " L" + str(
            point_x) + "," + str(point_y) + " L" + str(start_x_2) + "," + str(start_y_2) + " z"
    else:
        point_x = (
            radius * factors[1]) * math.cos(math.radians(360.0 * (arrow_end / total_length) - 90))
        point_y = (
            radius * factors[1]) * math.sin(math.radians(360.0 * (arrow_end / total_length) - 90))
        start_x_1 = (radius * factors[0]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_y_1 = (radius * factors[0]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_x_2 = (radius * factors[2]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        start_y_2 = (radius * factors[2]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90))
        end_x_1 = (
            radius * factors[0]) * math.cos(math.radians(360.0 * ((shoulder) / total_length) - 90))
        end_y_1 = (
            radius * factors[0]) * math.sin(math.radians(360.0 * ((shoulder) / total_length) - 90))
        end_x_2 = (
            radius * factors[2]) * math.cos(math.radians(360.0 * ((shoulder) / total_length) - 90))
        end_y_2 = (
            radius * factors[2]) * math.sin(math.radians(360.0 * ((shoulder) / total_length) - 90))
        feature_path = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius) + "," + str(radius) + param_1 + str(end_x_1) + "," + str(end_y_1) + " L" + str(
            point_x) + "," + str(point_y) + " L" + str(end_x_2) + "," + str(end_y_2) + "A" + str(radius) + "," + str(radius) + param_2 + str(start_x_2) + "," + str(start_y_2) + " z"
    coord_path.append(feature_path)
    return coord_path


def create_rectangle_path_circular(
        radius, coord_dict, total_length, track_ratio):
    coord_path = []
    coord_path.append("exon")
    feat_strand = coord_dict['feat_strand']
    factors = compute_factors(feat_strand, track_ratio)
    feat_start = coord_dict['feat_start']
    feat_end = coord_dict['feat_end']
    rect_strand_dict = {
        "positive": [
            " 0 0 1 ", " 0 0 0 "], "negative": [
            " 0 0 0 ", " 0 0 1 "]}
    param_1, param_2 = rect_strand_dict[feat_strand]
    start_x_1 = (
        radius * factors[0]) * math.cos(math.radians(360.0 * (feat_start / total_length) - 90))
    start_y_1 = (
        radius * factors[0]) * math.sin(math.radians(360.0 * (feat_start / total_length) - 90))
    start_x_2 = (
        radius * factors[2]) * math.cos(math.radians(360.0 * (feat_start / total_length) - 90))
    start_y_2 = (
        radius * factors[2]) * math.sin(math.radians(360.0 * (feat_start / total_length) - 90))
    end_x_1 = (radius * factors[0]) * \
        math.cos(math.radians(360.0 * ((feat_end) / total_length) - 90))
    end_y_1 = (radius * factors[0]) * \
        math.sin(math.radians(360.0 * ((feat_end) / total_length) - 90))
    end_x_2 = (radius * factors[2]) * \
        math.cos(math.radians(360.0 * ((feat_end) / total_length) - 90))
    end_y_2 = (radius * factors[2]) * \
        math.sin(math.radians(360.0 * ((feat_end) / total_length) - 90))
    feature_path = "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius) + "," + str(radius) + param_1 + str(end_x_1) + "," + str(
        end_y_1) + " L" + str(end_x_2) + "," + str(end_y_2) + "A" + str(radius) + "," + str(radius) + param_2 + str(start_x_2) + "," + str(start_y_2) + " z"
    coord_path.append(feature_path)
    return coord_path


def create_gene_path_circular(gene_object, total_length, radius, track_ratio):
    coords = gene_object.location
    gene_body_color = gene_object.color
    coordinates_paths = []
    min_arrow_length = 30
    max_arrow_length = 500
    param_a = 3
    param_b = 5
    cds_arrow_length = (min_arrow_length + (max_arrow_length - min_arrow_length)
                        * 1 / (1 + math.exp(- param_a * (math.log10(total_length) - param_b))))

    for coord in coords:
        coord_path = []
        coord_dict = {
            'feat_type': coord[0],
            'feat_strand': coord[2],
            'feat_start': coord[3],
            'feat_end': coord[4]}
        feat_type = coord_dict['feat_type']
        if feat_type == "intron":
            coord_path = create_intron_path_circular(
                radius, coord_dict, total_length, track_ratio)
        elif feat_type == "exon":
            if coord[5]:
                coord_path = create_arrowhead_path_circular(
                    radius, coord_dict, total_length, cds_arrow_length, track_ratio)
            elif coord[5] == False:
                coord_path = create_rectangle_path_circular(
                    radius, coord_dict, total_length, track_ratio)
        coordinates_paths.append(coord_path)
    return coordinates_paths


def add_circular_axis(radius):
    axis_group = Group(id="axis")
    circular_axis = Circle(
        center=(
            0,
            0),
        r=radius,
        stroke='gray',
        stroke_width=1,
        fill='none')
    axis_group.add(circular_axis)
    return axis_group


def create_tick_paths_list(ticks, tick_width, size, total_len, radius):
    tick_paths_list = []
    ratio = {'small': [1.06, 1.08], 'large': [1.06, 1.11]}
    prox, dist = ratio[size]
    for tick in ticks:
        prox_x = (radius * prox) * \
            math.cos(math.radians(360.0 * (tick / total_len) - 90))
        prox_y = (radius * prox) * \
            math.sin(math.radians(360.0 * (tick / total_len) - 90))
        dist_x = (radius * dist) * \
            math.cos(math.radians(360.0 * (tick / total_len) - 90))
        dist_y = (radius * dist) * \
            math.sin(math.radians(360.0 * (tick / total_len) - 90))
        tick_path_desc = "M " + \
            str(prox_x) + "," + str(prox_y) + " L" + str(dist_x) + "," + str(dist_y) + " z"
        tick_path = Path(
            d=tick_path_desc,
            stroke='gray',
            stroke_width=tick_width)
        tick_paths_list.append(tick_path)
    return tick_paths_list


def set_tick_label_anchor_value(tick, total_len):
    angle = (360.0 * (tick / total_len))
    if 0 <= angle < 45:
        anchor_value, baseline_value = "middle", "text-after-edge"
    elif 45 <= angle < 155:
        anchor_value, baseline_value = "start", "middle"
    elif 155 <= angle < 205:
        anchor_value, baseline_value = "middle", "hanging"
    elif 205 <= angle < 315:
        anchor_value, baseline_value = "end", "middle"
    elif 315 <= angle < 360:
        anchor_value, baseline_value = "middle", "text-after-edge"
    else:
        raise ValueError("Abnormal angle: verify the ticks and total length")
    return anchor_value, baseline_value


def create_tick_label_paths_list(ticks, size, total_len, radius):
    tick_label_paths_list = []
    def base_to_kbp(tick): return str(int(tick / 1000)) + " kbp"
    ratio = {'small': [1.10, 1.11], 'large': [1.12, 1.14]}
    prox, dist = ratio[size]
    for tick in ticks:
        anchor_value, baseline_value = set_tick_label_anchor_value(
            tick, total_len)
        label_text = base_to_kbp(tick)
        tixk_label_x = (radius * prox) * \
            math.cos(math.radians(360.0 * (tick / total_len) - 90))
        tick_label_y = (radius * prox) * \
            math.sin(math.radians(360.0 * (tick / total_len) - 90))
        label_path = Text(
            label_text,
            insert=(
                tixk_label_x,
                tick_label_y),
            stroke='none',
            fill='black',
            font_size='10pt',
            font_weight="normal",
            font_family="Arial",
            text_anchor=anchor_value,
            dominant_baseline=baseline_value)
        tick_label_paths_list.append(label_path)
    return tick_label_paths_list


def ticks_cicular(gb_record, radius):
    total_len = len(gb_record.seq)
    if total_len <= 100000:
        tick_large, tick_small = 10000, 1000
    elif 100000 < total_len <= 1000000:
        tick_large, tick_small = 50000, 10000
    elif 1000000 < total_len <= 5000000:
        tick_large, tick_small = 200000, 50000
    elif 5000000 < total_len:
        tick_large, tick_small = 1000000, 200000
    tick_group = Group(id="tick")
    ticks_large = list(range(0, total_len, tick_large))
    tick_small = list(range(0, total_len, tick_small))
    ticks_small = [x for x in tick_small if x % tick_large != 0]
    tick_paths_large = create_tick_paths_list(
        ticks_large, 2.0, "large", total_len, radius)
    tick_paths_small = create_tick_paths_list(
        ticks_small, 2.0, "small", total_len, radius)
    ticks_large_nonzero = [x for x in ticks_large if x != 0]
    tick_label_paths_large = create_tick_label_paths_list(
        ticks_large_nonzero, "large", total_len, radius)
    tick_label_paths_small = create_tick_label_paths_list(
        ticks_small, "small", total_len, radius)
    for tick_path_large in tick_paths_large:
        tick_group.add(tick_path_large)
    for tick_label_path_large in tick_label_paths_large:
        tick_group.add(tick_label_path_large)
    return tick_group


def record_definition(record, radius):
    end_x_1 = (radius) * math.cos(math.radians(360.0 * 0 - 90))
    end_x_2 = (radius) * math.cos(math.radians(360.0 * (0.5) - 90))
    end_y_1 = (radius) * math.cos(math.radians(360.0 * (0.25) - 90))
    end_y_2 = (radius) * math.cos(math.radians(360.0 * (0.75) - 90))
    title_x = ((end_x_1 + end_x_2) / 2)
    title_y = ((end_y_1 + end_y_2) / 2)
    interval = 20
    fontsize = '{}pt'.format(18)
    cds_list = []
    for feature in record.features:
        if feature.type == "source":
            if 'organism' in feature.qualifiers.keys():
                record_name = feature.qualifiers['organism'][0]
            else:
                record_name = record.annotations["organism"]
            if 'isolate' in feature.qualifiers.keys():
                strain_name = feature.qualifiers['isolate'][0]
            elif 'strain' in feature.qualifiers.keys():
                strain_name = feature.qualifiers['strain'][0]
            else:
                strain_name = "unknown"
        if feature.type == "CDS":
            cds_list.append(feature)
    track_id = str(strain_name).replace(" ", "_")
    definition_group = Group(id=track_id)
    recoed_length = len(record.seq)
    num_cds = len(cds_list)
    accession = str(record.annotations["accessions"][0])
    length_text = "{:,} bp".format(recoed_length)
    cds_label = "{:,} CDS".format(num_cds)
    gc_percent = calculate_gc_percent(record.seq)
    gc_percent_text = str(gc_percent) + "% GC"

    name_path = Text(record_name,
                     insert=(title_x,
                             (title_y + (interval * -2))),
                     stroke='none',
                     fill='black',
                     font_size=fontsize,
                     font_style="italic",
                     font_weight="bold",
                     font_family="Arial",
                     text_anchor="middle")
    strain_path = Text(strain_name,
                       insert=(title_x,
                               (title_y + (interval * -0))),
                       stroke='none',
                       fill='black',
                       font_size=fontsize,
                       font_weight="bold",
                       font_family="Arial",
                       text_anchor="middle")
    accession_path = Text(accession,
                          insert=(title_x,
                                  (title_y + (interval * 2))),
                          stroke='none',
                          fill='black',
                          font_size=fontsize,
                          font_weight="normal",
                          font_family="Arial",
                          text_anchor="middle")
    length_path = Text(length_text,
                       insert=(title_x,
                               (title_y + (interval * 4))),
                       stroke='none',
                       fill='black',
                       font_size=fontsize,
                       font_weight="normal",
                       font_family="Arial",
                       text_anchor="middle")
    gc_percent_path = Text(gc_percent_text,
                           insert=(title_x,
                                   (title_y + (interval * 6))),
                           stroke='none',
                           fill='black',
                           font_size=fontsize,
                           font_weight="normal",
                           font_family="Arial",
                           text_anchor="middle")
    definition_group.add(name_path)
    definition_group.add(strain_path)
    definition_group.add(accession_path)
    definition_group.add(length_path)
    definition_group.add(gc_percent_path)

    return definition_group


def record_circular(gb_record, radius, track_ratio, color_table):
    total_len = len(gb_record.seq)
    record_group = Group(id="record")
    feature_dict = create_feature_dict(gb_record, color_table)
    for key in feature_dict:
        feature_object = feature_dict[key]
        if isinstance(feature_object, GeneObject):
            cds_paths = create_gene_path_circular(
                feature_object, total_len, radius, track_ratio)
            for cds_path in cds_paths:
                feat_type = cds_path[0]
                if feat_type == "exon":
                    exon_path = Path(
                        d=cds_path[1],
                        fill=feature_object.color,
                        stroke='none',
                        stroke_width=0.5)
                    record_group.add(exon_path)
                elif feat_type == "intron":
                    intron_path = Path(
                        d=cds_path[1],
                        stroke='gray',
                        fill="none",
                        stroke_width=1.0)
                    record_group.add(intron_path)
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
                repeat_path = create_rectangle_path_circular(
                    radius, coord_dict, total_len, track_ratio)
                unit_path = Path(
                    d=repeat_path[1],
                    fill=color,
                    stroke='none',
                    stroke_width=0.5)
                record_group.add(unit_path)
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
                feature_path = create_rectangle_path_circular(
                    radius, coord_dict, total_len, track_ratio)
                feature_path = Path(
                    d=feature_path[1],
                    fill=color,
                    stroke='none',
                    stroke_width=0.5)
                record_group.add(feature_path)
    return record_group


def skew_circular(gb_record, radius, df, track_width):
    total_len = len(gb_record.seq)
    norm_factor = 0.85
    norm_radius = radius * norm_factor
    skew_group = Group(id="skew")
    column = "GC skew"
    skew_desc = skew_path_circular(
        radius, df, total_len, norm_factor, track_width)
    circle_desc = circle_path_circular(radius, total_len, norm_factor)
    circle_path = ClipPath(id='clipper_circle')
    circle_path.add(Path(d=circle_desc, fill="white", stroke='none'))
    skew_group.add(circle_path)
    skew_high = Path(
        d=skew_desc,
        fill="#40e0d0",
        stroke='none',
        fill_opacity=1,
        fill_rule="evenodd")
    skew_low = Path(
        d=skew_desc,
        fill="#8a2be2",
        stroke='none',
        fill_opacity=1,
        clip_path="url(#clipper_circle)",
        clip_rule="nonzero",
        fill_rule="evenodd")

    skew_group.add(skew_high)
    skew_group.add(skew_low)

    return skew_group


def gc_content_circular(gb_record, radius, df, track_width):
    total_len = len(gb_record.seq)
    norm_factor = 0.65
    norm_radius = radius * norm_factor
    gc_group = Group(id="gc_content")
    gc_path_desc = gc_path_circular(
        radius, df, total_len, norm_factor, track_width)
    gc_path = Path(
        d=gc_path_desc,
        fill="#808080",
        stroke='none',
        fill_opacity=1,
        fill_rule="evenodd")
    gc_group.add(gc_path)
    return gc_group


def create_canvas(file_name, total_width, total_height):
    dwg = svgwrite.Drawing(
        filename=file_name +
        ".svg",
        size=(
            str(total_width) +
            'px',
            str(total_height) +
            'px'),
        viewBox=(
            '0 0 ' +
            str(total_width) +
            ' ' +
            str(total_height)))
    return dwg


def plot_genome_diagram(gb_record, window, step, dinucleotide, color_table):
    record_name = gb_record.annotations["source"]
    accession = str(gb_record.annotations["accessions"][0])
    total_width = 1000
    total_height = 1000
    radius = 300
    track_ratio = 0.19
    track_width = radius * track_ratio
    df = skew_df(gb_record, window, step, dinucleotide)
    canvas = create_canvas(record_name, total_width, total_height)
    axis_group = add_circular_axis(radius)
    axis_group.translate(total_width / 2, total_height / 2)

    record_group = record_circular(gb_record, radius, track_ratio, color_table)
    definition_group = record_definition(gb_record, radius)
    tick_group = ticks_cicular(gb_record, radius)
    skew_group = skew_circular(gb_record, radius, df, track_width)
    gc_group = gc_content_circular(gb_record, radius, df, track_width)

    definition_group.translate(total_width / 2, total_height / 2)
    tick_group.translate(total_width / 2, total_height / 2)
    record_group.translate(total_width / 2, total_height / 2)
    skew_group.translate(total_width / 2, total_height / 2)
    gc_group.translate(total_width / 2, total_height / 2)

    canvas.add(definition_group)
    canvas.add(tick_group)
    canvas.add(record_group)
    canvas.add(axis_group)
    canvas.add(skew_group)
    canvas.add(gc_group)
    outfile_name = accession + ".svg"
    canvas.filename = outfile_name
    canvas.save()


def main():
    args = _get_args()
    input_file = args.input
    dinucleotide = args.nt
    window = args.window
    step = args.step
    color_table = args.table
    column_names = ['feature_type','qualifier_key','value','color']
    if color_table == "":
        pass
    else:
        color_table = pd.read_csv(color_table,sep='\t',names=(column_names))
    gb_records = SeqIO.parse(input_file, 'genbank')
    for gb_record in gb_records:
        plot_genome_diagram(gb_record, window, step, dinucleotide, color_table)


if __name__ == "__main__":
    main()
