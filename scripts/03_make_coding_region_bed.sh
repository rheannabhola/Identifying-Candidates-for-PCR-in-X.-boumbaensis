#!/bin/bash
#SBATCH --job-name=make_coding_bed
#SBATCH --account=rrg-ben
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=make_coding_bed_%j.out
#SBATCH --error=make_coding_bed_%j.err

module load StdEnv/2023
module load bedtools/2.31.0

# Purpose:
# Convert X. laevis CDS vs X. boumbaensis BLAST hits into BED-style intervals.
# These intervals define putative coding/CDS-like regions in the female X. boumbaensis reference.

WORKDIR="/home/rheanna/projects/rrg-ben/rheanna/genomewide_coding_coverage"

BLAST="${WORKDIR}/XL_CDS_vs_boumb_genome.tsv"

RAW_BED="${WORKDIR}/XL_CDS_vs_boumb_genome.raw.bed"
FILTERED_BED="${WORKDIR}/XL_CDS_vs_boumb_genome.filtered.bed"
MERGED_BED="${WORKDIR}/genomewide_coding_regions.filtered.merged.bed"

cd "$WORKDIR"

# BLAST outfmt columns used here:
# 1 qseqid
# 2 sseqid
# 3 pident
# 4 length
# 5 mismatch
# 6 gapopen
# 7 qstart
# 8 qend
# 9 sstart
# 10 send
# 11 evalue
# 12 bitscore
# 13 qlen

# Convert BLAST subject coordinates to BED coordinates.
# BED start is 0-based, so subtract 1 from the smaller subject coordinate.
awk 'BEGIN{OFS="\t"} {
  start = ($9 < $10 ? $9 : $10)
  end = ($9 > $10 ? $9 : $10)
  bed_start = start - 1
  if (bed_start < 0) bed_start = 0

  print $2, bed_start, end, $1, $3, $4, $11, $12, $13
}' "$BLAST" \
| sort -k1,1 -k2,2n \
> "$RAW_BED"

# Filter to higher-confidence CDS-like hits:
# percent identity >= 80
# alignment length >= 50 bp
# e-value <= 1e-10
awk 'BEGIN{OFS="\t"} $5 >= 80 && $6 >= 50 && $7 <= 1e-10 {
  print
}' "$RAW_BED" \
> "$FILTERED_BED"

# Merge overlapping coding-like intervals for genome-wide depth calculation.
cut -f1-3 "$FILTERED_BED" \
| sort -k1,1 -k2,2n \
| bedtools merge \
> "$MERGED_BED"

echo "Wrote:"
echo "$RAW_BED"
echo "$FILTERED_BED"
echo "$MERGED_BED"
