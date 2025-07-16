import streamlit as st
import pandas as pd
import base64
from Bio.Blast import NCBIXML
from svgwrite import Drawing
from svgwrite.container import Group
from svgwrite.shapes import Line
from svgwrite.text import Text


# ==== Utility functions ====

def base_to_kbp(tick, seq_len):
    if seq_len < 5_000:
        return f"{tick} bp"
    elif seq_len < 1_000_000:
        return f"{tick // 1_000} kbp"
    else:
        return f"{tick // 1_000_000} Mbp"


def tick_size(seq_len):
    if seq_len <= 1_000:
        return 100, 10
    elif seq_len <= 5_000:
        return 500, 100
    elif seq_len <= 10_000:
        return 1_000, 500
    elif seq_len <= 20_000:
        return 2_000, 1_000
    elif seq_len <= 50_000:
        return 5_000, 1_000
    elif seq_len <= 100_000:
        return 10_000, 1_000
    elif seq_len <= 1_000_000:
        return 50_000, 10_000
    else:
        return 1_000_000, 200_000


def hsps_to_df(hsps, evalue_threshold, bitscore_threshold, identity_threshold):
    data = []
    for hsp in hsps:
        pident = 100 * (hsp.identities / hsp.align_length)
        if (
            hsp.expect <= evalue_threshold
            and hsp.bits >= bitscore_threshold
            and pident >= identity_threshold
        ):
            data.append([
                hsp.query_start, hsp.query_end,
                hsp.sbjct_start, hsp.sbjct_end,
                pident, hsp.frame
            ])
    return pd.DataFrame(data, columns=['qstart', 'qend', 'sstart', 'send', 'pident', 'frame'])


# ==== Drawing functions ====

def create_canvas(file_name, total_width, total_height):
    dwg = Drawing(
        filename=file_name,
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

def draw_hits(df, query_len, subject_len, total_width, total_height, color_plus, color_minus, stroke_width):
    hit_group = Group(id="hits")
    for _, row in df.iterrows():
        qstart = int(total_width * row['qstart'] / query_len)
        qend = int(total_width * row['qend'] / query_len)
        sstart = int(total_height * row['sstart'] / subject_len)
        send = int(total_height * row['send'] / subject_len)
        color = color_plus if (row['frame'][0] > 0 and row['frame'][1] > 0) or (row['frame'][0] < 0 and row['frame'][1] < 0) else color_minus
        opacity = min(1, row['pident'] / 100)
        hit = Line(
            start=(qstart, sstart), end=(qend, send),
            stroke=color, stroke_opacity=opacity,
            stroke_width=stroke_width, fill='none'
        )
        hit_group.add(hit)
    return hit_group


def create_frame(total_width, total_height, stroke_width=1, stroke_color="black"):
    frame_group = Group(id="frame")
    frame_path_desc = f"M0,0 L{total_width},0 L{total_width},{total_height} L0,{total_height} z"
    frame_line = Drawing().path(d=frame_path_desc,
                                stroke=stroke_color,
                                stroke_width=stroke_width,
                                fill="none")
    frame_group.add(frame_line)
    return frame_group


def create_grids(
    seq_len, total_width, total_height, axis,
    grid_color_large, grid_width_large, grid_opacity_large,
    grid_color_small, grid_width_small, grid_opacity_small
):
    group = Group(id="grids")
    tick_large, tick_small = tick_size(seq_len)
    grids_large = list(range(0, seq_len, tick_large))
    grids_small = [x for x in range(0, seq_len, tick_small) if x % tick_large != 0]

    # 大グリッド線
    for tick in grids_large:
        if tick == 0:
            continue
        if axis == "horizontal":
            x1, y1, x2, y2 = total_width * tick / seq_len, 0, total_width * tick / seq_len, total_height
        else:
            x1, y1, x2, y2 = 0, total_height * tick / seq_len, total_width, total_height * tick / seq_len
        grid_line = Line(
            start=(x1, y1), end=(x2, y2),
            stroke=grid_color_large,
            stroke_width=grid_width_large,
            stroke_opacity=grid_opacity_large
        )
        group.add(grid_line)

    # 小グリッド線
    for tick in grids_small:
        if axis == "horizontal":
            x1, y1, x2, y2 = total_width * tick / seq_len, 0, total_width * tick / seq_len, total_height
        else:
            x1, y1, x2, y2 = 0, total_height * tick / seq_len, total_width, total_height * tick / seq_len
        grid_line = Line(
            start=(x1, y1), end=(x2, y2),
            stroke=grid_color_small,
            stroke_width=grid_width_small,
            stroke_opacity=grid_opacity_small
        )
        group.add(grid_line)

    return group


def create_ticks_labels(seq_len, total_span, axis, tick_color, font_size, label_color, font_family):
    group = Group(id="ticks_labels")
    tick_large, _ = tick_size(seq_len)
    ticks_large = list(range(0, seq_len + 1, tick_large))

    for tick in ticks_large:
        if axis == "horizontal":
            x1, y1, x2, y2 = total_span * tick / seq_len, -10, total_span * tick / seq_len, 0
            text_x, text_y = total_span * tick / seq_len, -15
            angle = f"rotate(0,{text_x},{text_y})"
        else:
            x1, y1, x2, y2 = -10, total_span * tick / seq_len, 0, total_span * tick / seq_len
            text_x, text_y = -15, total_span * tick / seq_len
            angle = f"rotate(-90,{text_x},{text_y})"

        # Always draw tick line
        tick_line = Line(start=(x1, y1), end=(x2, y2), stroke=tick_color, stroke_width=1)
        group.add(tick_line)

        # Skip label for tick=0
        if tick != 0:
            label_text = base_to_kbp(tick, seq_len)
            label = Text(
                label_text,
                insert=(text_x, text_y),
                fill=label_color,
                font_size=f"{font_size}px",
                font_family=font_family,
                text_anchor="middle",
                transform=angle
            )
            group.add(label)

    return group


def create_definition_label(text, total_span, axis, font_size, font_family):
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

    label = Text(
        text,
        insert=(text_x, text_y),
        fill="black",
        font_size=f"{font_size}px",
        font_family=font_family,
        text_anchor=anchor_value,
        dominant_baseline=baseline_value,
        transform=angle
    )
    return label


def generate_dotplot(
    blast_record, output_file, figsize, offset,
    evalue_threshold, bitscore_threshold, identity_threshold,
    color_plus, color_minus, stroke_width,
    grid_color_large, grid_width_large, grid_opacity_large,
    grid_color_small, grid_width_small, grid_opacity_small,
    tick_color, label_color, label_font_size,
    frame_color, frame_width, font_family
):
    if not blast_record.alignments:
        return None

    query_len = blast_record.query_length
    subject_len = blast_record.alignments[0].length
    max_len = max(query_len, subject_len)
    total_width = int(figsize * query_len / max_len)
    total_height = int(figsize * subject_len / max_len)

    canvas = Drawing(filename=output_file, size=(f"{total_width + offset}px", f"{total_height + offset}px"))

    # Grids
    grids_horizontal = create_grids(
        query_len, total_width, total_height, "horizontal",
        grid_color_large, grid_width_large, grid_opacity_large,
        grid_color_small, grid_width_small, grid_opacity_small
    )
    grids_horizontal.translate(offset, offset)
    canvas.add(grids_horizontal)

    grids_vertical = create_grids(
        subject_len, total_width, total_height, "vertical",
        grid_color_large, grid_width_large, grid_opacity_large,
        grid_color_small, grid_width_small, grid_opacity_small
    )
    grids_vertical.translate(offset, offset)
    canvas.add(grids_vertical)

    # Ticks and labels
    ticks_labels_horizontal = create_ticks_labels(query_len, total_width, "horizontal", tick_color, label_font_size, label_color, font_family)
    ticks_labels_horizontal.translate(offset, offset)
    canvas.add(ticks_labels_horizontal)

    ticks_labels_vertical = create_ticks_labels(subject_len, total_height, "vertical", tick_color, label_font_size, label_color, font_family)
    ticks_labels_vertical.translate(offset, offset)
    canvas.add(ticks_labels_vertical)

    # Definition labels
    query_label = create_definition_label(blast_record.query, total_width, "horizontal", label_font_size + 2, font_family)
    query_label.translate(offset, 20)
    canvas.add(query_label)

    subject_label = create_definition_label(blast_record.alignments[0].hit_def, total_height, "vertical", label_font_size + 2, font_family)
    subject_label.translate(20, 20)
    canvas.add(subject_label)

    # Frame
    frame_group = create_frame(total_width, total_height, stroke_width=frame_width, stroke_color=frame_color)
    frame_group.translate(offset, offset)
    canvas.add(frame_group)

    # Hits
    for alignment in blast_record.alignments:
        df = hsps_to_df(
            alignment.hsps,
            evalue_threshold, bitscore_threshold, identity_threshold
        )
        if df.empty:
            continue
        hits = draw_hits(
            df, query_len, subject_len,
            total_width, total_height,
            color_plus, color_minus, stroke_width
        )
        hits.translate(offset, offset)
        canvas.add(hits)
    return canvas


# ==== Streamlit App ====

st.title("blast2dotplot")
st.markdown(
    "This app generates a dotplot from BLASTN/TBLASTX XML results (-outfmt 5). "
    "It visualizes the alignment of sequences in a two-dimensional plot, "
    "where each dot represents a hit between the query and subject sequences."
)

# Session state for SVG
if "svg_bytes" not in st.session_state:
    st.session_state.svg_bytes = None
    st.session_state.output_file = "dotplot.svg"

uploaded_file = st.file_uploader("Upload BLAST XML (-outfmt 5)", type="xml")

output_file = st.text_input("Output SVG file name:", st.session_state.output_file)
figsize = st.sidebar.slider("Figure Size (px)", 500, 3000, 1000, step=10)
frame_width = st.sidebar.slider("Frame Line Width", 1, 5, 1)
stroke_width = st.sidebar.slider("Hit Line Width", 1, 10, 2)
label_font_size = st.sidebar.slider("Label Font Size", 8, 20, 16, step=1)
font_family = st.sidebar.selectbox("Font Family", ["Arial"], index=0)
offset = st.sidebar.slider("Offset (px)", 10, 200, 80, step=1)

frame_color = st.sidebar.color_picker("Frame Color", "#000000")
tick_color = st.sidebar.color_picker("Tick Color", "#000000")
label_color = st.sidebar.color_picker("Label Color", "#000000")
grid_color_large = st.sidebar.color_picker("Large Grid Color", "#D3D3D3")
grid_color_small = st.sidebar.color_picker("Small Grid Color", "#D3D3D3")
color_plus = st.sidebar.color_picker("Color (+ frame)", "#1f77b4")
color_minus = st.sidebar.color_picker("Color (- frame)", "#ff7f0e")

evalue_threshold = st.sidebar.number_input("Max E-value", value=1e-2, format="%.1e")
bitscore_threshold = st.sidebar.number_input("Min Bitscore", value=50.0)
identity_threshold = st.sidebar.number_input("Min Identity (%)", value=30.0)

if st.button("RUN"):
    if uploaded_file is None:
        st.error("Please upload a BLAST XML file first.")
    else:
        blast_records = list(NCBIXML.parse(uploaded_file))
        if blast_records:
            svg_canvas = generate_dotplot(
                blast_records[0], output_file, figsize, offset,
                evalue_threshold, bitscore_threshold, identity_threshold,
                color_plus, color_minus, stroke_width,
                grid_color_large, 1, 1.0,
                grid_color_small, 1, 0.5,
                tick_color, label_color, label_font_size,
                frame_color, frame_width, font_family
            )
            svg_bytes = svg_canvas.tostring().encode('utf-8')
            st.session_state.svg_bytes = svg_bytes
            st.session_state.output_file = output_file
            st.success("Dotplot generated!")

# Display existing SVG if available
if st.session_state.svg_bytes:
    b64_svg = base64.b64encode(st.session_state.svg_bytes).decode('utf-8')
    st.markdown(
        f'<div style="text-align:center;"><img src="data:image/svg+xml;base64,{b64_svg}"/></div>',
        unsafe_allow_html=True
    )
    st.download_button(
        label="Download SVG",
        data=st.session_state.svg_bytes,
        file_name=st.session_state.output_file,
        mime="image/svg+xml"
    )

st.markdown("---")
st.markdown(
    "Author: [Satoshi Kawato](https://github.com/satoshikawato)  |  "
    "Source: [bio_small_scripts](https://github.com/satoshikawato/bio_small_scripts)  |  "
    "Blog post: [Qiita](https://qiita.com/satoshi_kawato/items/0e5c13621e53bad8d9a0)",
    unsafe_allow_html=True
)
