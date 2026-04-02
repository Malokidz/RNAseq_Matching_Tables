#!/bin/bash

# ==========================================
# Extract longest isoform per Trinity gene
# ==========================================

# Default values
THREADS=1

# Help function
usage() {
cat << EOF

Usage:
  $(basename "$0") -i <input_fasta> -o <output_fasta> [-t threads]

Description:
  Extracts the longest isoform per Trinity gene (DNxxx_cX_gY)
  using sequence length.

Required arguments:
  -i    Input FASTA file (e.g. transcripts.fasta)
  -o    Output FASTA file (longest isoforms)

Optional arguments:
  -t    Number of threads for seqkit (default: 1)
  -h    Show this help message

Example:
  $(basename "$0") -i transcripts.fasta -o longest_isoforms.fasta -t 8

EOF
exit 1
}

# Parse arguments
while getopts "i:o:t:h" opt; do
  case $opt in
    i) INPUT=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Check required arguments
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input exists
if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file '$INPUT' not found"
    exit 1
fi

# Temporary files
TMP_IDS=$(mktemp)

echo "🔍 Extracting longest isoform IDs..."

# Step 1: get longest isoform IDs
seqkit fx2tab -l -n -i -j "$THREADS" "$INPUT" | \
awk '{
    gene=$1
    sub(/_i[0-9]+.*/, "", gene)

    if(!(gene in max) || $2 > max[gene]) {
        max[gene]=$2
        id[gene]=$1
    }
}
END {
    for(g in id) print id[g]
}' > "$TMP_IDS"

echo "📦 Extracting sequences..."

# Step 2: extract sequences
seqkit grep -f "$TMP_IDS" -j "$THREADS" "$INPUT" > "$OUTPUT"

# Cleanup
rm -f "$TMP_IDS"

echo "✅ Done!"
echo "Output file: $OUTPUT"
