#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def calculate_n_stretches(fasta_file):
    gap_sizes = {}
    
    # Parse each scaffold sequence in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        scaffold_id = record.id
        sequence = str(record.seq)
        
        # Find all contiguous stretches of "N"s
        n_sizes = []
        gap_length = 0
        
        for base in sequence:
            if base.upper() == "N":
                gap_length += 1  # Increment gap length if "N" is found
            elif gap_length > 0:
                n_sizes.append(gap_length)  # Save the gap length when "N" sequence ends
                gap_length = 0
        
        # Check for any trailing "N" gaps at the end of the scaffold
        if gap_length > 0:
            n_sizes.append(gap_length)
        
        # Save gap sizes for each scaffold
        gap_sizes[scaffold_id] = n_sizes
    
    return gap_sizes

def main():
    parser = argparse.ArgumentParser(description="Calculate 'N' gap sizes in each scaffold of a FASTA file.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file")
    args = parser.parse_args()
    
    n_stretches = calculate_n_stretches(args.fasta_file)
    
    # Display results
    for scaffold, gaps in n_stretches.items():
        total_n_gap_size = sum(gaps)
        print(f"Scaffold: {scaffold}, N-Gap Sizes: {gaps}, Total N-Gap Size: {total_n_gap_size}")

if __name__ == "__main__":
    main()
