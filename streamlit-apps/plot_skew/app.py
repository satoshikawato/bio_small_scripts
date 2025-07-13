#!/usr/bin/env python
# coding: utf-8


import sys
import argparse
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
from io import StringIO, BytesIO

def _get_args():
    parser = argparse.ArgumentParser(
        description='Generate GC skew plot in SVG. Multiple SVG files produced if multi-FASTA file is provided')
    parser.add_argument(
        '-i',
        '--input',
        help='Fasta (required)',
        type=str,
        required=True)
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
    #args = parser.parse_args(args=['--input','input.fa'])
    return args


def fasta_to_records(in_fa):
    seq_records = [record for record in SeqIO.parse(in_fa, 'fasta')]
    return seq_records


def calculate_dinucleotide_skew(seq, base1, base2):
    "Return dinucleotide skew in a given sequence"
    base1_count = seq.count(base1)
    base2_count = seq.count(base2)
    skew = ((base1_count - base2_count) / (base1_count + base2_count))
    return skew


def sliding_window(seq, window, step):
    for start in range(0, len(seq), step):
        end = start + window
        if end > len(seq):
            break
        out_seq = seq[start:end]
        yield start, out_seq


def skew_df(record, window, step, nt):
    nt_list = list(nt)
    nt_1 = nt_list[0]
    nt_2 = nt_list[1]
    skew_sum = 0
    skew_dict = {}
    skew_cumulative_dict = {}
    sequence = record.seq.upper()
    for start, seq_part in sliding_window(sequence, window, step):
        skew = calculate_dinucleotide_skew(seq_part, nt_1, nt_2)
        skew_dict[start] = skew
        skew_sum = (skew_sum + skew)
        skew_cumulative_dict[start] = (skew_sum)
    max_skew_abs = abs(max(skew_dict.values(), key=abs))
    max_skew_cumulative_abs = abs(max(skew_cumulative_dict.values(), key=abs))
    factor = (max_skew_abs / max_skew_cumulative_abs)
    skew_cumulative_dict.update((x, (y * factor))
                                for x, y in skew_cumulative_dict.items())
    skew_legend = "{} skew".format(nt)
    cumulative_skew_legend = "Cumulative {} skew, normalized".format(nt)
    df = pd.DataFrame({skew_legend: pd.Series(skew_dict),
                       cumulative_skew_legend: pd.Series(skew_cumulative_dict)})
    return df

def draw_plot(record, df, width_px=1200, height_px=400, font_size=20, dpi=100):

    fig_width_in = width_px / dpi
    fig_height_in = height_px / dpi

    fig = plt.figure(figsize=(fig_width_in, fig_height_in), dpi=dpi)
    ax = fig.add_subplot(111)
    record_id = record.id
    ax.set_title(record_id, fontsize=font_size)
    ax.set_xlim([0, len(record.seq)])
    ax.ticklabel_format(style='plain', axis='both')
    ax.tick_params(axis="both", labelsize=font_size)
    coordinate_of_max_cumulative_skew = df.iloc[:, 1].idxmax()
    coordinate_of_min_cumulative_skew = df.iloc[:, 1].idxmin()
    max_cumulative_skew = df.iloc[:, 1].max()
    min_cumulative_skew = df.iloc[:, 1].min()
    fig.tight_layout()
    ax.plot(df)
    ax.axvline(x=coordinate_of_max_cumulative_skew, color="red")
    ax.text(
        x=coordinate_of_max_cumulative_skew,
        y=max_cumulative_skew * 0.8,
        s=coordinate_of_max_cumulative_skew,
        fontsize=font_size)
    ax.axvline(x=coordinate_of_min_cumulative_skew, color="red")
    ax.text(
        x=coordinate_of_min_cumulative_skew,
        y=min_cumulative_skew * 1.4,
        s=coordinate_of_min_cumulative_skew,
        fontsize=font_size)
    ax.legend(df.columns, fontsize=font_size)
    return fig


def main():
    args = _get_args()
    in_fa = args.input
    window = args.window
    step = args.step
    nt = args.nt.upper()
    if len(nt) != 2:
        nt = "GC"
    nt_list = list(nt)
    window = 1000
    step = 100
    records = fasta_to_records(in_fa)
    for record in records:
        df = skew_df(record, window, step, nt)
        draw_plot(record, df)

# === Streamlit UI ===
st.title("SkewPlot v0.1.0")

uploaded_file = st.file_uploader("Upload FASTA file")
# Window and step size inputs
col1, col2, col3 = st.columns(3)
with col1:
    nt = st.text_input("Dinucleotide (e.g. GC, AT)", value="GC")
with col2:
    window = st.number_input("Window size", min_value=10, value=1000)
with col3:
    step = st.number_input("Step size", min_value=1, value=100)

# Width, height, and font size inputs
col4, col5, col6 = st.columns(3)
with col4:
    width_px = st.number_input("Width (px)", min_value=400, max_value=4000, value=2000)
with col5:
    height_px = st.number_input("Height (px)", min_value=200, max_value=2000, value=400)
with col6:
    font_size = st.number_input("Font size", min_value=6, max_value=40, value=20, step=1)


# Initialize session state for plot buffers
if "plot_bufs" not in st.session_state:
    st.session_state.plot_bufs = {}

# === Process and Plot ===
if uploaded_file and st.button("Plot"):
    fasta_text = uploaded_file.getvalue().decode("utf-8")

    # Check if the uploaded file is a valid FASTA format
    if not fasta_text.lstrip().startswith(">"):
        st.error("❌ The uploaded file does not appear to be a valid FASTA file.")
    else:

        seq_records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
        st.session_state.plot_bufs.clear()  # Reset previous plots

        total_records = len(seq_records)
        progress_bar = st.progress(0, text="Processing sequences...")

        for idx, record in enumerate(seq_records, 1):
            with st.spinner(f"Processing {record.id} ({idx}/{total_records})"):
                df = skew_df(record, window, step, nt)
                fig = draw_plot(record, df, width_px, height_px, font_size)

                # Export to BytesIO
                buf = BytesIO()
                fig.savefig(buf, format="png", dpi=300, bbox_inches='tight')
                buf.seek(0)

                # Satore in session state
                st.session_state.plot_bufs[record.id] = buf

            # Refresh progress bar
            progress_bar.progress(idx / total_records, text=f"Processed {idx}/{total_records}")

        progress_bar.empty()  # Clear the progress bar
        st.success("✅ All sequences processed!")

# === Preview and Download  ===
for record_id, buf in st.session_state.plot_bufs.items():
    st.subheader(f"Record: {record_id}")
    st.image(buf, use_container_width=True)
    st.download_button(
        label=f"Download {record_id}_skew.png",
        data=buf,
        file_name=f"{record_id}_skew.png",
        mime="image/png",
        key=f"download_{record_id}"
    )

st.markdown("---")
st.markdown(
    "Author: [Satoshi Kawato](https://github.com/satoshikawato)  |  "
    "Source: [bio_small_scripts](https://github.com/satoshikawato/bio_small_scripts)",
    unsafe_allow_html=True
)
