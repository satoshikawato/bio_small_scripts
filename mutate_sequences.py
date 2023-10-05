#!/usr/bin/env python
# coding: utf-8
import os
import argparse
import random
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Mutation function
def mutate_sequence(seq, seq_type, mutation_percentage, conserved_sites):
    if seq_type == "dna":
        alphabet = "ACTG"
    elif seq_type == "protein":
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
    
    variable_sites = len([i for i in range(len(seq)) if i not in conserved_sites])
    mutation_count = int(len(seq) * mutation_percentage / 100)
    mutation_count = min(variable_sites, mutation_count)
    print(f"Length of the sequence: {len(seq)}")
    print(f"Maximum mutated residues: {mutation_count} ({100* mutation_count/len(seq)}%) ")
        
    seq = list(seq)
    reserved_sites = list(conserved_sites)  # Initialize with conserved sites for every mutation round

    # Ensure the starting methionine of the original sequence is preserved
    first_residue_position = next((i for i, residue in enumerate(seq) if residue != '-'), None)  # Find the first non-gap position
    if seq[first_residue_position] == 'M':
        reserved_sites.append(first_residue_position)
    
    mutations = []  # Create a list to store mutation details
    actual_mutation_count = 0
    for _ in range(mutation_count):
        valid_positions = [i for i in range(len(seq)) if i not in reserved_sites and seq[i] != '-']
        if not valid_positions:
            print("Warning: No valid positions left for mutation!")
            break
        mutation_position = random.choice(valid_positions)
        #mutation_position = random.choice([i for i in range(len(seq)) if i not in reserved_sites and seq[i] != '-'])  # Ensure that gaps are not selected for mutation
        reserved_sites.append(mutation_position)
        original_residue = seq[mutation_position]
        new_residue = random.choice([residue for residue in alphabet if residue != original_residue])
        seq[mutation_position] = new_residue
        actual_mutation_count += 1
        mutations.append((mutation_position + 1, original_residue, new_residue))  # mutation_position + 1 because sequence positions typically start from 1
    print(f"Number of mutated residues: {actual_mutation_count} ({100* actual_mutation_count/len(seq):.3f}%)")
    mutated_seq = "".join(seq)
    mutated_seq_no_gaps = mutated_seq.replace("-", "")  # Remove gaps from the mutated sequence
    return mutated_seq_no_gaps, mutations  # Return mutations as well
    
# Identify conserved sites
def get_conserved_sites(alignment):
    conserved_sites = []

    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        if '-' not in column and len(set(column)) == 1:
            conserved_sites.append(i)

    return conserved_sites

# Determine sequence type
def determine_seq_type(seq):
    dna_alphabet = "ACTGRYSWKMBDHVN"
    protein_alphabet = "ACDEFGHIKLMNPQRSTVWYBXZ"
    for residue in seq:
        if residue not in dna_alphabet and residue in protein_alphabet:
            return "protein"
    return "dna"

# Main function
def main(input_file, output_file, entry_name, mutation_percentage, num_mutated_seqs):
    alignment = AlignIO.read(input_file, "fasta")
    conserved_sites = get_conserved_sites(alignment)

    seq_type = determine_seq_type(str(alignment[0].seq))

    found_entry = False
    mutated_records = []
    base_filename, _ = os.path.splitext(output_file)
    tsv_file = f"{base_filename}.tsv"
    with open(tsv_file, 'w') as tsv_out:
        tsv_out.write("Sequence\tResidue\tOriginal\tMutated\n")  # Header for the TSV file
        for record in alignment:
            if record.id == entry_name:
                found_entry = True
                for i in range(num_mutated_seqs):
                    mutated_seq, mutations = mutate_sequence(record.seq, seq_type, mutation_percentage, conserved_sites)  # Unpack the returned mutations
                    mutated_records.append(SeqRecord(Seq(mutated_seq), id=f"{entry_name}_mutated_{i+1}", description="Mutated sequence"))
                    
                    # Sort mutations by their position before writing to TSV
                    for mutation in sorted(mutations, key=lambda x: x[0]):
                        tsv_out.write(f"{entry_name}_mutated_{i+1}\t{mutation[0]}\t{mutation[1]}\t{mutation[2]}\n")


    if not found_entry:
        print(f"Error: Entry name '{entry_name}' not found in the input file.")
        return

    SeqIO.write(mutated_records, output_file, "fasta")
    print(f"{num_mutated_seqs} mutated sequences saved to {output_file}")
    print(f"Details of mutations saved to {tsv_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Introduce mutations in a sequence from a multiple sequence alignment.")
    parser.add_argument("-i", "--input", required=True, help="Input multiple sequence alignment file in FASTA format.")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with the mutated sequence(s).")
    parser.add_argument("-e", "--entry", required=True, help="Entry name of the sequence to be mutated.")
    parser.add_argument("-p", "--percentage", type=int, required=True, help="Percentage of mutations to be introduced.")
    parser.add_argument("-n", "--num", type=int, default=1, help="Number of mutated sequences to generate.")

    args = parser.parse_args()
    main(args.input, args.output, args.entry, args.percentage, args.num)
