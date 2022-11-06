#!/usr/bin/env python
# coding: utf-8


import sys
import pandas as pd
import svgwrite
import argparse
import logging
from svgwrite.container import Group
from svgwrite.shapes import Line
from svgwrite.path import Path
from svgwrite.text import Text
from Bio.Blast import NCBIXML


logging.basicConfig(level=logging.INFO, format='%(message)s')


def _get_args():
    parser = argparse.ArgumentParser(
        description='Draw a dot plot based on a pairwise BLASTN/TBLASTX result')
    parser.add_argument(
        "--input",
        "--in",
        "-i",
        metavar="FILE",
        help="input BLASTN/TBLASTX result file in XML format (-outfmt 5)",
        required=True)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args


def draw_hits(df, query_len, subject_len, total_width, total_height, max_len):
    hit_group = Group(id="hits")
    for index, row in df.iterrows():
        qstart = int(total_width * int(row['qstart']) / query_len)
        qend = int(total_width * int(row['qend']) / query_len)
        sstart = int(total_height * int(row['sstart']) / subject_len)
        send = int(total_height * int(row['send']) / subject_len)
        pident = int(row['pident'])
        if pident <= 60:
            opacity = 0.4
        elif 60 < pident <= 70:
            opacity = 0.6
        elif 70 < pident <= 80:
            opacity = 0.8
        elif 80 < pident <= 100:
            opacity = 1
        if (row['frame'][0] > 0 and row['frame'][1] > 0) or (
                row['frame'][0] < 0 and row['frame'][1] < 0):
            color = "#1f77b4"
        else:
            color = "#ff7f0e"
        hit = Line(
            start=(
                qstart,
                sstart),
            end=(
                qend,
                send),
            stroke=color,
            stroke_opacity=opacity,
            stroke_width=2,
            fill='none')
        hit_group.add(hit)
    return hit_group


def create_tick_label_paths_list(
        ticks, tick_width, size, total_span, seq_len, axis):
    tick_label_paths_list = []

    anchor_value = "middle"
    baseline_value = "middle"
    for tick in ticks:
        if axis == "horizontal":
            tick_label_x = total_span * (tick / seq_len)
            tick_label_y = 0
            angle = 'rotate(0,0, 0)'
        elif axis == "vertical":
            tick_label_x = 0
            tick_label_y = total_span * (tick / seq_len)
            angle = 'rotate({},{}, {})'.format(-90, tick_label_x, tick_label_y)
        label_text = base_to_kbp(tick, seq_len)
        label_path = Text(label_text, insert=(tick_label_x, tick_label_y),
                          stroke='none',
                          fill='black',
                          font_size='10pt',
                          font_weight="normal",
                          font_family="Arial",
                          text_anchor=anchor_value,
                          dominant_baseline=baseline_value, transform=angle)
        tick_label_paths_list.append(label_path)
    return tick_label_paths_list


def base_to_kbp(tick, seq_len):
    if seq_len < 5000:
        tick_string = str(tick) + " bp"
    elif 5000 <= seq_len < 1000000:
        tick_string = str(int(tick / 1000)) + " kbp"
    else:
        tick_string = str(int(tick / 1000000)) + " Mbp"
    return tick_string


def create_tick_paths_list(ticks, tick_width, size, total_span, seq_len, axis):
    tick_paths_list = []
    ratio = {'small': [80, 70], 'large': [80, 60]}
    prox, dist = ratio[size]
    for tick in ticks:
        if axis == "horizontal":
            x1 = total_span * (tick / seq_len)
            y1 = prox
            x2 = total_span * (tick / seq_len)
            y2 = dist
        elif axis == "vertical":
            x1 = dist
            y1 = total_span * (tick / seq_len)
            x2 = prox
            y2 = total_span * (tick / seq_len)
        tick_path_desc = "M " + str(x1) + "," + \
            str(y1) + " L" + str(x2) + "," + str(y2)
        tick_path = Path(
            d=tick_path_desc,
            stroke='black',
            stroke_width=tick_width)
        tick_paths_list.append(tick_path)
    return tick_paths_list


def create_grid_paths_list(ticks, tick_width, size,
                           total_width, total_height, seq_len, axis):
    grid_paths_list = []
    for tick in ticks:
        if axis == "horizontal":
            x1 = total_width * (tick / seq_len)
            y1 = 0
            x2 = total_width * (tick / seq_len)
            y2 = total_height
        elif axis == "vertical":
            x1 = 0
            y1 = total_height * (tick / seq_len)
            x2 = total_width
            y2 = total_height * (tick / seq_len)
        grid_path_desc = "M " + str(x1) + "," + \
            str(y1) + " L" + str(x2) + "," + str(y2)
        grid_path = Path(
            d=grid_path_desc,
            stroke='#D3D3D3',
            stroke_width=tick_width)
        grid_paths_list.append(grid_path)
    return grid_paths_list


def create_ticks(seq_len, total_span, max_len, axis):
    tick_large, tick_small = tick_size(seq_len)
    tick_group = Group(id="tick")
    ticks_large = list(range(0, seq_len, tick_large))
    ticks_small = [
        x for x in list(
            range(
                0,
                seq_len,
                tick_small)) if x %
        tick_large != 0]
    ticks_large_nonzero = [
        x for x in list(
            range(
                0,
                seq_len,
                tick_large)) if x != 0]
    tick_paths_large = create_tick_paths_list(
        ticks_large, 1, "large", total_span, seq_len, axis)
    tick_paths_small = create_tick_paths_list(
        ticks_small, 1, "small", total_span, seq_len, axis)
    for tick_path_large in tick_paths_large:
        tick_group.add(tick_path_large)
    for tick_path_small in tick_paths_small:
        tick_group.add(tick_path_small)
    return tick_group


def create_labels(seq_len, total_span, max_len, axis):
    tick_large, tick_small = tick_size(seq_len)
    label_group = Group(id="label")
    ticks_large = list(range(0, seq_len, tick_large))
    tick_small = list(range(0, seq_len, tick_small))
    ticks_small = [x for x in tick_small if x % tick_large != 0]
    ticks_large_nonzero = [x for x in ticks_large if x != 0]
    tick_label_paths_large = create_tick_label_paths_list(
        ticks_large, 1, "large", total_span, seq_len, axis)
    for tick_label_path_large in tick_label_paths_large:
        label_group.add(tick_label_path_large)
    return label_group


def tick_size(seq_len):
    if seq_len <= 1000:
        tick_large, tick_small = 100, 10
    elif 1000 < seq_len <= 5000:
        tick_large, tick_small = 500, 100
    elif 5000 < seq_len <= 10000:
        tick_large, tick_small = 1000, 500
    elif 10000 < seq_len <= 20000:
        tick_large, tick_small = 2000, 1000
    elif 20000 < seq_len <= 50000:
        tick_large, tick_small = 5000, 1000
    elif 50000 < seq_len <= 100000:
        tick_large, tick_small = 10000, 1000
    elif 100000 < seq_len <= 1000000:
        tick_large, tick_small = 50000, 10000
    elif 1000000 < seq_len <= 5000000:
        tick_large, tick_small = 1000000, 200000
    elif 5000000 < seq_len:
        tick_large, tick_small = 1000000, 500000
    return tick_large, tick_small


def create_grids(seq_len, total_width, total_height, max_len, axis):
    tick_large, tick_small = tick_size(seq_len)
    grid_group = Group(id="grid")
    grids_large = list(range(0, seq_len, tick_large))
    grid_small = list(range(0, seq_len, tick_small))
    grids_small = [x for x in grid_small if x % tick_large != 0]
    grids_large_nonzero = [x for x in grids_large if x != 0]
    grid_paths_large = create_grid_paths_list(
        grids_large_nonzero,
        1,
        "large",
        total_width,
        total_height,
        seq_len,
        axis)
    grid_paths_small = create_grid_paths_list(
        grids_small, 1, "small", total_width, total_height, seq_len, axis)
    for grid_path_large in grid_paths_large:
        grid_group.add(grid_path_large)
    for grid_path_small in grid_paths_small:
        grid_group.add(grid_path_small)
    return grid_group


def create_frame(total_width, total_height):
    frame_group = Group(id="frame")
    frame_path_desc = "M " + str(0) + "," + str(0) + " L" + str(total_width) + "," + str(0) + " L" + str(
        total_width) + "," + str(total_height) + " L" + str(0) + "," + str(total_height) + "z"
    frame_path = Path(
        d=frame_path_desc,
        stroke='black',
        stroke_width=1,
        fill="none")
    frame_group.add(frame_path)
    return frame_group


def create_canvas(file_name, total_width, total_height):
    dwg = svgwrite.Drawing(filename=file_name + ".svg", size=(str(total_width) + 'px', str(total_height) + 'px'),
                           viewBox=('0 0 ' + str(total_width) + ' ' + str(total_height)))
    return dwg


def hsps_to_df(hsps):
    header = ['qstart', 'qend', 'sstart', 'send', 'pident', 'frame']
    hsp_lines = []
    for hsp in hsps:
        pident = 100 * (hsp.identities / hsp.align_length)
        hsp_line = [
            hsp.query_start,
            hsp.query_end,
            hsp.sbjct_start,
            hsp.sbjct_end,
            pident,
            hsp.frame]
        hsp_lines.append(hsp_line)
    df = pd.DataFrame(hsp_lines, columns=header)
    return df


def create_definition(seq_len, total_span, max_len, accession, title, axis):
    definition_group = Group(id="definiton")
    anchor_value = "middle"
    baseline_value = "middle"
    if axis == "horizontal":
        text_x = total_span * 0.5
        text_y = 0
        angle = 'rotate(0,0, 0)'
    elif axis == "vertical":
        text_x = 0
        text_y = total_span * 0.5
        angle = 'rotate({},{}, {})'.format(-90, text_x, text_y)

    definition_text = "{} {}".format(accession, title)
    definition_path = Text(definition_text, insert=(text_x, text_y),
                           stroke='none',
                           fill='black',
                           font_size='15pt',
                           font_weight="normal",
                           font_family="Arial",
                           text_anchor=anchor_value,
                           dominant_baseline=baseline_value, transform=angle)

    definition_group.add(definition_path)

    return definition_group


def main():
    offset = 80
    figsize = 1000
    args = _get_args()
    in_file = args.input
    supported = ["BLASTN", "TBLASTX"]
    with open(in_file, 'rt') as infile:
        blast_records = NCBIXML.parse(infile)
        msg = ("Input file: {}".format(in_file))
        logging.info(msg)
        for blast_record in blast_records:
            application = blast_record.application
            if application in supported:
                msg = ("{} result detectted".format(application))
                logging.info(msg)
            else:
                raise TypeError(
                    "{} is not yet supported. Exiting...".format(application))

            query_id = blast_record.query_id
            query_title = blast_record.query
            query_len = blast_record.query_length
            msg = (
                "Query:\n    ID: {}\n    Title: {}\n    Length: {}".format(
                    query_id, query_title, query_len))
            logging.info(msg)
            for alignment in blast_record.alignments:

                subject_id = alignment.accession
                subject_title = alignment.hit_def
                subject_len = alignment.length
                msg = (
                    "Subject:\n    ID: {}\n    Title: {}\n    Length: {}".format(
                        query_id, query_title, query_len))
                logging.info(msg)
                filename = "{}_{}_{}".format(query_id, subject_id, application)
                max_len = max(query_len, subject_len)
                total_width = int(figsize * query_len / max_len)
                total_height = int(figsize * subject_len / max_len)

                df = hsps_to_df(alignment.hsps)

                canvas = create_canvas(
                    filename, total_width + offset, total_height + offset)

                grids_horizontal = create_grids(
                    query_len, total_width, total_height, max_len, "horizontal")
                grids_horizontal.translate(offset, offset)
                canvas.add(grids_horizontal)

                grids_vertical = create_grids(
                    subject_len, total_width, total_height, max_len, "vertical")
                grids_vertical.translate(offset, offset)
                canvas.add(grids_vertical)

                hits = draw_hits(
                    df,
                    query_len,
                    subject_len,
                    total_width,
                    total_height,
                    max_len)
                hits.translate(offset, offset)
                canvas.add(hits)

                ticks_hotizontal = create_ticks(
                    query_len, total_width, max_len, "horizontal")
                ticks_hotizontal.translate(offset, 0)
                canvas.add(ticks_hotizontal)

                ticks_vertical = create_ticks(
                    subject_len, total_height, max_len, "vertical")
                ticks_vertical.translate(0, offset)
                canvas.add(ticks_vertical)

                labels_horizontal = create_labels(
                    query_len, total_width, max_len, "horizontal")
                labels_horizontal.translate(offset, 50)
                canvas.add(labels_horizontal)

                labels_vertical = create_labels(
                    subject_len, total_height, max_len, "vertical")
                labels_vertical.translate(50, offset)
                canvas.add(labels_vertical)

                definition_horizontal = create_definition(
                    query_len, total_width, max_len, query_id, query_title, "horizontal")
                definition_horizontal.translate(offset, 20)
                canvas.add(definition_horizontal)

                definition_vertical = create_definition(
                    subject_len, total_height, max_len, subject_id, subject_title, "vertical")
                definition_vertical.translate(20, 20)
                canvas.add(definition_vertical)

                frame_group = create_frame(total_width, total_height)
                frame_group.translate(offset, offset)
                canvas.add(frame_group)
                canvas.filename = "{}.svg".format(filename)
                canvas.save()
                msg = ("Dotplot saved as: {}".format(canvas.filename))
                logging.info(msg)


if __name__ == "__main__":
    main()
