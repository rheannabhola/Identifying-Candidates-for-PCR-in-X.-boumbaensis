#!/bin/bash
#SBATCH --job-name=maleWGS_coding_depth
#SBATCH --account=rrg-ben
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=maleWGS_coding_depth_%j.out
#SBATCH --error=maleWGS_coding_depth_%j.err

module load StdEnv/2023
module load samtools

# Purpose:
# Calculate male X. boumbaensis WGS depth across genome-wide coding/CDS-like regions.
# The male WGS reads were mapped to the female X. boumbaensis reference genome.
# The -aa option is important because it reports positions with zero coverage.

WORKDIR="/home/rheanna/projects/rrg-ben/rheanna/genomewide_coding_coverage"

BAM="/home/rheanna/projects/rrg-ben/rheanna/boum_WGS_male/aligned/SRR35972782_male_vs_femaleBoumb.sorted.bam"

BED="${WORKDIR}/genomewide_coding_regions.filtered.merged.bed"

OUT="${WORKDIR}/maleWGS_depth_genomewide_coding.tsv"

cd "$WORKDIR"

samtools depth \
  -aa \
  -q 20 \
  -Q 20 \
  -b "$BED" \
  "$BAM" \
  > "$OUT"

echo "Finished male WGS depth calculation:"
echo "$OUT"
