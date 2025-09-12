#!/bin/bash
#BSUB -P acc_oscarlr
#BSUB -q premium
#BSUB -n 10
#BSUB -W 72:00
#BSUB -J extract_hifi
#BSUB -o extract_hifi.out
#BSUB -e extract_hifi.err

set -euo pipefail

scratch="/sc/arion/scratch/arayan01/projects/nanopore/data/Pacbio"
projectdata="/sc/arion/scratch/arayan01/projects/nanopore/data/oscardir"

# Paths for BAMs inside the tars
WGS_BAM="2024-04-17_human_wgs/1_A01/hifi_reads/m84248_240410_212543_s1.hifi_reads.bam"
SEM_BAM="data/2022-08-29_CW48_7-8/demultiplex.bc1058--bc1058.bam"
RS411_BAM="data/2022-09-21_Seq_CW48_10-12/demultiplex.bc1017--bc1017.bam"

# Extract all needed BAMs in one go per tar
echo "Extracting WGS BAM..."
tar -xvzf "${projectdata}/2024-04-17_human_wgs.tar.gz" -C "${scratch}/hifi-reads" "$WGS_BAM"

echo "Extracting SEM + RS411 BAMs..."
tar -xvzf "${projectdata}/data.tar.gz" -C "${scratch}/hifi-reads" "$SEM_BAM" "$RS411_BAM"

# Convert all BAMs to FASTQ with proper labels
declare -A bam_labels=(
  ["$WGS_BAM"]="SEM-HiFi-WGS"
  ["$SEM_BAM"]="SEM-Targeted-HiFi"
  ["$RS411_BAM"]="RS411-Targeted-HiFi"
)

for bam_path in "${!bam_labels[@]}"; do
    bam_file="${scratch}/hifi-reads/${bam_path}"
    fastq_file="${scratch}/hifi-reads/${bam_labels[$bam_path]}.fastq.gz"
    echo "Converting BAM → FASTQ: $bam_file → $fastq_file"
    samtools fastq "$bam_file" | pigz -p 36 > "$fastq_file"
done

echo "All BAMs extracted and converted!"
