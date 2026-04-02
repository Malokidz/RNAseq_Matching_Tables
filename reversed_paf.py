import sys

def reverse_paf(input_paf, output_paf):
    with open(input_paf, 'r') as infile, open(output_paf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            
            # Mandatory fields
            query_name, query_len, query_start, query_end = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            strand = fields[4]
            target_name, target_len, target_start, target_end = fields[5], int(fields[6]), int(fields[7]), int(fields[8])
            matches, aln_len, mapq = fields[9], fields[10], fields[11]
            
            # Optional tags (if any)
            tags = fields[12:] if len(fields) > 12 else []

            # Flip start and end for reverse strand
            if strand == '-':
                new_target_start = target_len - target_end
                new_target_end = target_len - target_start
                strand = '+'
            else:
                new_target_start = target_start
                new_target_end = target_end

            # Prepare reversed line including tags
            reversed_line = [
                target_name, str(target_len), str(new_target_start), str(new_target_end),
                strand, query_name, str(query_len), str(query_start), str(query_end),
                matches, aln_len, mapq
            ] + tags

            outfile.write('\t'.join(reversed_line) + '\n')

if __name__ == "__main__":
    # Check for correct usage
    if len(sys.argv) != 3:
        print("Usage: python reverse_paf.py <input.paf> <output.paf>")
        sys.exit(1)

    input_paf = sys.argv[1]
    output_paf = sys.argv[2]
    reverse_paf(input_paf, output_paf)

