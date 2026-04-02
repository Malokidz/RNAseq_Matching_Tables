#!/usr/bin/env python3

import argparse
import os
import glob
import subprocess
import sys

def run(cmd):
    print(f"\nRunning:\n{cmd}\n")
    subprocess.run(cmd, shell=True, check=True)

def index_exists(prefix):
    """
    Check if Bowtie2 index files exist
    """
    required_files = [
        f"{prefix}.1.bt2",
        f"{prefix}.2.bt2",
        f"{prefix}.3.bt2",
        f"{prefix}.4.bt2",
        f"{prefix}.rev.1.bt2",
        f"{prefix}.rev.2.bt2"
    ]
    return all(os.path.exists(f) for f in required_files)

def main():

    parser = argparse.ArgumentParser(
        description="Bowtie2 RNA-seq mapping + flagstat summary pipeline"
    )

    parser.add_argument("--index", required=True,
                        help="FASTA file OR existing Bowtie2 index prefix")

    parser.add_argument("--reads", required=True,
                        help="Directory containing *_R1*.fastq.gz and *_R2*.fastq.gz")

    parser.add_argument("--outdir", required=True,
                        help="Output directory")

    parser.add_argument("--threads", type=int, default=16,
                        help="Number of threads (default=16)")

    args = parser.parse_args()

    index_input = args.index
    reads_path = args.reads
    outdir = args.outdir
    threads = args.threads

    os.makedirs(outdir, exist_ok=True)

    # --------------------------------------------------
    # STEP 1 — Check or Build Bowtie2 index
    # --------------------------------------------------

    if os.path.isfile(index_input) and index_input.endswith((".fa", ".fasta", ".fna")):
        # User gave FASTA → build index inside outdir
        index_prefix = os.path.join(outdir, "transcripts_index")

        if index_exists(index_prefix):
            print("Index already exists. Skipping indexing step.")
        else:
            print("Index not found. Building Bowtie2 index...")
            run(f"bowtie2-build {index_input} {index_prefix}")

    else:
        # Assume user gave an existing index prefix
        if not index_exists(index_input):
            print("ERROR: Provided index prefix does not exist or is incomplete.")
            sys.exit(1)

        print("Using existing Bowtie2 index.")
        index_prefix = index_input

    # --------------------------------------------------
    # STEP 2 — Map reads directly to sorted BAM
    # --------------------------------------------------

    r1_files = sorted(glob.glob(os.path.join(reads_path, "*_R1*.fastq.gz")))

    if not r1_files:
        print("No R1 files found. Check reads directory.")
        sys.exit(1)

    for r1 in r1_files:
        r2 = r1.replace("_R1", "_R2")

        if not os.path.exists(r2):
            print(f"Missing pair for {r1}")
            continue

        sample = os.path.basename(r1).split("_R1")[0]
        sorted_bam = os.path.join(outdir, f"{sample}_sorted.bam")

        if os.path.exists(sorted_bam):
            print(f"{sample} already mapped. Skipping.")
            continue

        run(f"""
        bowtie2 -x {index_prefix} \
            -1 {r1} -2 {r2} \
            -p {threads} | \
        samtools sort -@ {threads} -o {sorted_bam}
        """)

        run(f"samtools index {sorted_bam}")

    # --------------------------------------------------
    # STEP 3 — Flagstat
    # --------------------------------------------------

    sorted_bams = glob.glob(os.path.join(outdir, "*_sorted.bam"))

    for bam in sorted_bams:
        flagstat_out = bam.replace(".bam", ".flagstat.txt")

        if os.path.exists(flagstat_out):
            continue

        run(f"samtools flagstat {bam} > {flagstat_out}")

    # --------------------------------------------------
    # STEP 4 — Merge flagstats
    # --------------------------------------------------

    merged_file = os.path.join(outdir, "merged_flagstat.tsv")

    with open(merged_file, "w") as out:
        out.write("Sample\tTotal\tMapped\tMapped_%\tProperly_Paired\tProperly_Paired_%\n")

        flagstat_files = glob.glob(os.path.join(outdir, "*_sorted.flagstat.txt"))

        for f in flagstat_files:
            sample = os.path.basename(f).replace("_sorted.flagstat.txt", "")

            with open(f) as infile:
                lines = infile.readlines()

            total = lines[0].split()[0]
            mapped_line = [l for l in lines if " mapped (" in l][0]
            proper_line = [l for l in lines if "properly paired" in l][0]

            mapped = mapped_line.split()[0]
            mapped_pct = mapped_line.split("(")[1].split("%")[0]

            proper = proper_line.split()[0]
            proper_pct = proper_line.split("(")[1].split("%")[0]

            out.write(f"{sample}\t{total}\t{mapped}\t{mapped_pct}\t{proper}\t{proper_pct}\n")

    print("\nPipeline completed successfully.")
    print(f"Final merged table: {merged_file}")

if __name__ == "__main__":
    main()
