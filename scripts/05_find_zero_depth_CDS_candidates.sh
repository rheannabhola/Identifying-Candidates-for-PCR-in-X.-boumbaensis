#!/bin/bash
#SBATCH --job-name=zeroDepth_CDS_candidates
#SBATCH --account=rrg-ben
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=zeroDepth_CDS_candidates_%j.out
#SBATCH --error=zeroDepth_CDS_candidates_%j.err

module load StdEnv/2023
module load bedtools/2.31.0

# Purpose:
# Identify CDS-like regions in the female X. boumbaensis reference genome
# that have 100% zero coverage in male X. boumbaensis WGS.
#
# Inputs:
# 1. maleWGS_depth_genomewide_coding.tsv
#    Per-base male WGS depth across coding-like regions.
#
# 2. XL_CDS_vs_boumb_genome.filtered.bed
#    Filtered X. laevis CDS BLAST hits against the X. boumbaensis reference.
#
# Output:
# A strict candidate table requiring:
# - alignment length >= 300 bp
# - percent identity >= 90%
# - fraction zero in male WGS = 1.0

WORKDIR="/home/rheanna/projects/rrg-ben/rheanna/genomewide_coding_coverage"

DEPTH="${WORKDIR}/maleWGS_depth_genomewide_coding.tsv"
CDS_BED="${WORKDIR}/XL_CDS_vs_boumb_genome.filtered.bed"

ZERO_BED="${WORKDIR}/maleWGS_zeroDepth_coding_merged.bed"
CDS_WITH_IDS="${WORKDIR}/XL_CDS_vs_boumb_genome.filtered.withIDs.bed"
OVERLAP="${WORKDIR}/XL_CDS_hits_overlap_maleZeroDepth.tsv"
SUMMARY="${WORKDIR}/maleWGS_zeroDepth_CDS_hit_summary.tsv"
STRICT="${WORKDIR}/maleWGS_zeroDepth_CDS_candidates_100pctZero_min300bp_pid90.tsv"

cd "$WORKDIR"

# 1. Extract positions where male WGS depth is exactly zero.
# samtools depth is 1-based, BED is 0-based, so position becomes position-1 to position.
awk 'BEGIN{OFS="\t"} $3 == 0 {print $1, $2-1, $2}' "$DEPTH" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> "$ZERO_BED"

# 2. Add unique IDs to each filtered CDS BLAST hit.
# Columns:
# contig, start, end, hit_id, laevis_CDS_ID, pident, align_len, evalue, bitscore, qlen/hit_length
awk 'BEGIN{OFS="\t"} {
  id = "hit_" NR
  len = $3 - $2
  print $1, $2, $3, id, $4, $5, $6, $7, $8, len
}' "$CDS_BED" \
> "$CDS_WITH_IDS"

# 3. Intersect CDS hits with male-zero regions and report overlap length.
bedtools intersect \
  -a "$CDS_WITH_IDS" \
  -b "$ZERO_BED" \
  -wo \
> "$OVERLAP"

# 4. Summarize zero-depth overlap per CDS hit.
awk 'BEGIN{OFS="\t"} {
  hit_id=$4
  contig=$1
  start=$2
  end=$3
  gene=$5
  pident=$6
  aln_len=$7
  evalue=$8
  bitscore=$9
  hit_len=$10
  overlap=$14

  key=hit_id

  if (!(key in seen)) {
    seen[key]=1
    c[key]=contig
    s[key]=start
    e[key]=end
    g[key]=gene
    pid[key]=pident
    alen[key]=aln_len
    eval[key]=evalue
    bit[key]=bitscore
    hlen[key]=hit_len
    zero_bp[key]=0
  }

  zero_bp[key] += overlap
}
END {
  print "contig","start","end","hit_id","laevis_CDS_ID","pident","align_len","evalue","bitscore","hit_length","zero_bp","frac_zero_maleWGS"

  for (key in seen) {
    frac = zero_bp[key] / hlen[key]
    print c[key],s[key],e[key],key,g[key],pid[key],alen[key],eval[key],bit[key],hlen[key],zero_bp[key],frac
  }
}' "$OVERLAP" \
> "$SUMMARY"

# 5. Strict candidate filter:
# - hit length >= 300 bp
# - percent identity >= 90%
# - 100% of hit has zero male WGS coverage
awk 'BEGIN{OFS="\t"} NR==1 || ($10 >= 300 && $6 >= 90 && $12 == 1)' "$SUMMARY" \
> "$STRICT"

echo "Wrote:"
echo "$ZERO_BED"
echo "$SUMMARY"
echo "$STRICT"
