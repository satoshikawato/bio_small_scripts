#!/usr/bin/env python
# coding: utf-8
import re
import argparse
import sys
import os

def update_svg_labels(fasta_path, tree_path, output_path):
    """
    Extracts species names from a FASTA file and updates accession numbers 
    in an SVG tree file to the format: 'Species Name Accession_ID'.
    The species name is wrapped in a tspan for italic styling.
    """
    id_to_species = {}
    
    # Regex to capture ID and [Species Name] from FASTA headers
    # Example: >XP_020502845.1 cytochrome P450 2K1-like [Labrus bergylta]
    fasta_pattern = re.compile(r'>(\S+).*?\[(.*?)\]')
    
    try:
        with open(fasta_path, 'r', encoding='utf-8') as f:
            content = f.read()
            matches = fasta_pattern.findall(content)
            for seq_id, species in matches:
                id_to_species[seq_id] = species
    except FileNotFoundError:
        print(f"Error: FASTA file not found at {fasta_path}")
        return

    if not id_to_species:
        print("Warning: No species data could be extracted from the FASTA file.")
        return

    # Load the SVG tree file
    try:
        with open(tree_path, 'r', encoding='utf-8') as f:
            svg_content = f.read()
    except FileNotFoundError:
        print(f"Error: SVG file not found at {tree_path}")
        return

    # Replace IDs with 'Species Name + ID'
    updated_count = 0
    for seq_id, species in id_to_species.items():
        # Target the ID inside the SVG text tags (e.g., >XP_12345.1<)
        old_label = f">{seq_id}<"
        # Use <tspan> with italic style for the species name
        new_label = f'><tspan style="font-style:italic">{species}</tspan> {seq_id}<'
        
        if old_label in svg_content:
            svg_content = svg_content.replace(old_label, new_label)
            updated_count += 1

    # Save the modified SVG
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(svg_content)

    print(f"Success: Updated {updated_count} labels.")
    print(f"Output saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Append italicized species names from a FASTA file to tip labels in an SVG phylogenetic tree."
    )
    
    # Define flagged arguments
    parser.add_argument("-f", "--fasta", required=True, help="Path to the input FASTA/FAA file")
    parser.add_argument("-s", "--svg", required=True, help="Path to the input SVG tree file")
    parser.add_argument("-o", "--out", required=True, help="Path for the output SVG file")

    args = parser.parse_args()

    update_svg_labels(args.fasta, args.svg, args.out)

if __name__ == "__main__":
    main()