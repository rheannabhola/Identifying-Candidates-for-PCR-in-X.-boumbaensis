#!/bin/bash
#SBATCH --job-name=XL_CDS_boumb
#SBATCH --account=rrg-ben
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=XL_CDS_vs_boumb_genome_%j.out
#SBATCH --error=XL_CDS_vs_boumb_genome_%j.err

module load StdEnv/2023
module load blast+

# Purpose:
# BLAST X. laevis CDS sequences against the female X. boumbaensis reference genome.
# This identifies coding/CDS-like regions in the X. boumbaensis assembly.

WORKDIR="/home/rheanna/projects/rrg-ben/rheanna/genomewide_coding_coverage"

QUERY="/home/rheanna/projects/rrg-ben/rheanna/XL_CDS_only_nospaces.fasta"
DB="${WORKDIR}/boumb_wbubble_genome_db"

OUT="${WORKDIR}/XL_CDS_vs_boumb_genome.tsv"

cd "$WORKDIR"

blastn \
  -query "$QUERY" \
  -db "$DB" \
  -evalue 1e-10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
  -num_threads 8 \
  -max_target_seqs 20 \
  -out "$OUT"

echo "Finished BLAST search:"
echo "$OUT"
