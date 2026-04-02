import sys
import requests

def get_gene_name(protein_name):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": protein_name,
        "fields": "gene_names,protein_name",
        "format": "tsv",
        "size": 1
    }
    response = requests.get(url, params=params)
    if response.ok:
        lines = response.text.strip().split("\n")
        if len(lines) > 1:
            return lines[1].split("\t")[0]
    return "Not found"

if len(sys.argv) != 3:
    print("Usage: python script.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    outfile.write("Protein Description\tGene Name\n")
    for line in infile:
        protein_desc = line.strip()
        if protein_desc:
            gene_name = get_gene_name(protein_desc)
            print(f"{protein_desc} -> {gene_name}")
            outfile.write(f"{protein_desc}\t{gene_name}\n")

