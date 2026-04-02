#!/usr/bin/env python3

import os
import sys
import pandas as pd

def parse_attributes(attributes_str):
    """Parse the attributes field into a dictionary."""
    attributes = {}
    for field in attributes_str.strip(";").split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            attributes[key.strip()] = value.strip()
    return attributes

def process_file(input_file, output_file):
    rows = []
    with open(input_file, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("#") or not line:
            i += 1
            continue
            
        parts = line.split("\t")
        if len(parts) < 9:
            i += 1
            continue

        scaffold = parts[0]
        feature_type = parts[2]
        attributes_str = parts[8]
        attributes = parse_attributes(attributes_str)

        if feature_type == "gene":
            gene_id = attributes.get("ID", "-")
            mRNA_id = "-"
            product = "-"
            go_terms = "-"
            domains = "-"
            cog = "-"
            eggnog = "-"
            signalp = "-"
            tm_domains = "-"
            merops = "-"

            # Check if next line is mRNA
            if i + 1 < len(lines):
                next_line = lines[i + 1].strip()
                if not next_line.startswith("#"):
                    next_parts = next_line.split("\t")
                    if len(next_parts) >= 9 and next_parts[2] == "mRNA":
                        mrna_attributes = parse_attributes(next_parts[8])
                        mRNA_id = mrna_attributes.get("ID", "-")
                        product = mrna_attributes.get("product", "-")
                        go_terms = mrna_attributes.get("Ontology_term", "-")

                        dbxref = mrna_attributes.get("Dbxref", "")
                        dbxref_items = [x.strip() for x in dbxref.split(",")] if dbxref else []

                        pfam_list = []
                        interpro_list = []
                        for item in dbxref_items:
                            if item.startswith("PFAM:"):
                                pfam_list.append(item)
                            if item.startswith("InterPro:"):
                                interpro_list.append(item)
                        domains = ";".join(pfam_list + interpro_list) if pfam_list or interpro_list else "-"

                        note = mrna_attributes.get("note", "")
                        note_items = [x.strip() for x in note.split(",")] if note else []

                        cog_list = []
                        eggnog_list = []
                        for item in note_items:
                            if item.startswith("COG:"):
                                cog_list.append(item)
                            if item.startswith("EggNog:"):
                                eggnog_list.append(item)
                            if "SignalP" in item:
                                signalp = item
                            if "TransMembrane" in item:
                                tm_domains = item
                            if "MEROPS" in item:
                                merops = item
                        
                        cog = ";".join(cog_list) if cog_list else "-"
                        eggnog = ";".join(eggnog_list) if eggnog_list else "-"

                        i += 1  # Skip the mRNA line since we processed it here

            rows.append([
                scaffold, gene_id, mRNA_id, product, go_terms,
                domains, cog, eggnog, signalp, tm_domains, merops
            ])
        i += 1

    # Create DataFrame
    columns = [
        "Scaffold", "Gene ID", "mRNA", "Product", "GO Terms",
        "Domains (PFAM/InterPro)", "COG", "EggNog",
        "Signal Peptide", "TM Domains", "MEROPS"
    ]
    df = pd.DataFrame(rows, columns=columns)

    # Save as Excel
    df.to_excel(output_file, index=False)
    print(f"Created: {output_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python create_annotation_excel_from_gff.py <input_file> <output_file.xlsx>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not os.path.isfile(input_file):
        print(f"Error: {input_file} does not exist.")
        sys.exit(1)

    process_file(input_file, output_file)

if __name__ == "__main__":
    main()
