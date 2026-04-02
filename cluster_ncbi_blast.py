import re
import csv
import sys
import os
import argparse

def parse_cluster_file(cluster_filename, gene_list_filename, output_filename, locus_tag):
    """
    Parse cluster file and create table with best matches for genes in the list
    Handles both base gene IDs (FUN_000096) and transcript IDs (FUN_000096-T1)
    """
    # Check if input files exist
    if not os.path.exists(cluster_filename):
        print(f"Error: Cluster file '{cluster_filename}' not found!")
        return
    if not os.path.exists(gene_list_filename):
        print(f"Error: Gene list file '{gene_list_filename}' not found!")
        return
    
    # Read gene list from file
    with open(gene_list_filename, 'r') as f:
        gene_ids = [line.strip() for line in f if line.strip()]
    
    print(f"Processing {len(gene_ids)} genes from {gene_list_filename}...")
    print(f"Using locus tag pattern: {locus_tag}")
    
    # Parse cluster file and select best cluster for each gene
    cluster_data = {}
    current_gene = None
    current_cluster = None
    
    # Create regex pattern for the locus tag
    # Escape any special characters in the locus tag and match the gene number
    escaped_locus_tag = re.escape(locus_tag)
    gene_pattern = re.compile(f'({escaped_locus_tag}\\d+)-T\\d+')
    
    with open(cluster_filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Check if this line starts a new query
            if line.startswith('Query #'):
                match = gene_pattern.search(line)
                if match:
                    current_gene = match.group(1)  # Extract base gene ID without transcript suffix
                    current_cluster = {'bit_score': -1}  # Initialize with low score
            
            # Extract cluster information for current gene
            elif current_gene:
                if line.startswith('Cluster:'):
                    current_cluster['Cluster'] = line.split(':', 1)[1].strip()
                elif line.startswith('Highest Bit Score:'):
                    try:
                        current_bit_score = float(line.split(':', 1)[1].strip())
                        current_cluster['Bit_Score'] = current_bit_score
                    except ValueError:
                        current_cluster['Bit_Score'] = 0
                elif line.startswith('Percent Identity:'):
                    current_cluster['Identity'] = line.split(':', 1)[1].strip().replace('%', '')
                
                # When we reach the end of a cluster section, compare and store the best one
                elif line.startswith('Accession Length:') or line.startswith('---') or line.startswith('>'):
                    if current_cluster and 'Cluster' in current_cluster and 'Bit_Score' in current_cluster:
                        # Only process genes that are in our list
                        if current_gene in gene_ids:
                            if current_gene not in cluster_data:
                                cluster_data[current_gene] = current_cluster.copy()
                            else:
                                if current_cluster['Bit_Score'] > cluster_data[current_gene]['Bit_Score']:
                                    cluster_data[current_gene] = current_cluster.copy()
                        
                        # Reset for next cluster
                        current_cluster = {'bit_score': -1}
    
    # Create output table
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        
        # Write header
        writer.writerow(['GeneId', 'Cluster', 'Percent_Identity', 'Highest_Bit_Score'])
        
        # Write data for each gene
        matched_count = 0
        for gene_id in gene_ids:
            if gene_id in cluster_data:
                data = cluster_data[gene_id]
                writer.writerow([
                    gene_id,
                    data.get('Cluster', 'No cluster'),
                    data.get('Identity', 'N/A'),
                    data.get('Bit_Score', 'N/A')
                ])
                matched_count += 1
            else:
                writer.writerow([gene_id, 'No match found', 'N/A', 'N/A'])
    
    print(f"Results saved to {output_filename}")
    print(f"Found matches for {matched_count} out of {len(gene_ids)} genes")
    
    # Also print summary to console
    if matched_count > 0:
        print("\nSummary Table (first few matches):")
        print("-" * 100)
        print(f"{'GeneId':<20} {'Cluster':<50} {'Identity':<15} {'Bit_Score':<15}")
        print("-" * 100)
        
        # Show first 5 matches for preview
        count = 0
        for gene_id in gene_ids:
            if count >= 5:  # Show only first 5 for preview
                break
            if gene_id in cluster_data:
                data = cluster_data[gene_id]
                cluster_name = data.get('Cluster', 'No cluster')
                # Truncate long cluster names for display
                if len(cluster_name) > 45:
                    cluster_name = cluster_name[:42] + "..."
                print(f"{gene_id:<20} {cluster_name:<50} {data.get('Identity', 'N/A'):<15} {data.get('Bit_Score', 'N/A'):<15}")
                count += 1

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Parse cluster file and create table with best matches for genes')
    parser.add_argument('cluster_file', help='Input cluster file (e.g., AWB7Z5TJ014-Alignment.txt)')
    parser.add_argument('gene_list_file', help='Input gene list file (e.g., unique.txt)')
    parser.add_argument('output_file', help='Output file (e.g., list1.tsv)')
    parser.add_argument('--locus_tag', '-l', default='FUN_', 
                       help='Locus tag prefix (default: FUN_)')
    
    # Parse arguments
    args = parser.parse_args()
    
    cluster_filename = args.cluster_file
    gene_list_filename = args.gene_list_file
    output_filename = args.output_file
    locus_tag = args.locus_tag
    
    print(f"Input cluster file: {cluster_filename}")
    print(f"Input gene list: {gene_list_filename}")
    print(f"Output file: {output_filename}")
    print(f"Locus tag: {locus_tag}")
    print("-" * 50)
    
    parse_cluster_file(cluster_filename, gene_list_filename, output_filename, locus_tag)

if __name__ == "__main__":
    main()
