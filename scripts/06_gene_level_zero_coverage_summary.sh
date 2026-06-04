#!/bin/bash
#SBATCH --job-name=gene_zero_summary
#SBATCH --account=rrg-ben
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=gene_zero_summary_%j.out
#SBATCH --error=gene_zero_summary_%j.err

module load StdEnv/2023
module load bedtools/2.31.0

# Purpose:
# Summarize male-zero coverage at the gene level.
#
# This script asks:
# For each non-LOC gene in the regions of interest, do all detected high-quality
# CDS-like intervals have 100% zero coverage in male X. boumbaensis WGS?
#
# Regions of interest:
# - all tig00003177
# - tig00004591 SNP-associated block: 5,175,617-10,105,961
#
# Important interpretation:
# These are CDS-like intervals based on X. laevis CDS BLAST hits,
# not confirmed X. boumbaensis exon annotations.

WORKDIR="/home/rheanna/projects/rrg-ben/rheanna/genomewide_coding_coverage"

CDS_BED="${WORKDIR}/XL_CDS_vs_boumb_genome.filtered.bed"
ZERO_BED="${WORKDIR}/maleWGS_zeroDepth_coding_merged.bed"

ROI_CDS="${WORKDIR}/ROI_CDS_hits_3177_4591SNPblock_pid90_min300.bed"
ROI_ZERO="${WORKDIR}/ROI_CDS_hits_100pct_maleZero.bed"
AWK_SCRIPT="${WORKDIR}/summarize_gene_zero_coverage.awk"
OUT="${WORKDIR}/gene_level_CDS_zeroCoverage_summary_noLOC.tsv"

cd "$WORKDIR"

# 1. Extract all high-quality CDS-like intervals in the regions of interest.
# Input CDS_BED columns:
# 1 contig
# 2 start
# 3 end
# 4 laevis_CDS_ID
# 5 pident
# 6 align_len
# 7 evalue
# 8 bitscore
# 9 qlen, if present from BLAST output
#
# Gene symbol is extracted from the X. laevis CDS ID by splitting on "_"
# and taking field 3, matching the naming pattern used in this project.

awk 'BEGIN{FS=OFS="\t"}
{
  contig=$1
  start=$2
  end=$3
  cds_id=$4
  pident=$5
  align_len=$6
  evalue=$7
  bitscore=$8

  split(cds_id, parts, "_")
  gene=parts[3]

  if (pident < 90) next
  if (align_len < 300) next

  if (contig=="tig00003177") {
    print contig,start,end,gene,cds_id,pident,align_len,evalue,bitscore
  }

  else if (contig=="tig00004591" && start >= 5175617 && end <= 10105961) {
    print contig,start,end,gene,cds_id,pident,align_len,evalue,bitscore
  }
}' "$CDS_BED" \
> "$ROI_CDS"

# 2. Identify which ROI CDS-like intervals are fully contained in male-zero regions.
# -f 1.0 requires 100% of the CDS-like interval to overlap male-zero coverage.

bedtools intersect \
  -a "$ROI_CDS" \
  -b "$ZERO_BED" \
  -f 1.0 \
  -u \
> "$ROI_ZERO"

# 3. Make an awk script to summarize by gene.

cat > "$AWK_SCRIPT" << 'EOF'
BEGIN {
  FS = OFS = "\t"
}

FNR == NR {
  key = $1 ":" $2 "-" $3 ":" $4 ":" $5
  zero[key] = 1
  next
}

{
  contig = $1
  start = $2
  end = $3
  gene = $4
  cds_id = $5
  pident = $6
  align_len = $7

  if (gene ~ /LOC/) next

  key = contig ":" start "-" end ":" gene ":" cds_id

  total[gene]++

  if (key in zero) {
    zero_count[gene]++
  }

  if (!(gene in seen)) {
    seen[gene] = 1
    contigs[gene] = contig
    min_start[gene] = start
    max_end[gene] = end
    max_pident[gene] = pident
    max_align[gene] = align_len
    intervals[gene] = contig ":" start "-" end
  } else {
    if (start < min_start[gene]) min_start[gene] = start
    if (end > max_end[gene]) max_end[gene] = end
    if (pident > max_pident[gene]) max_pident[gene] = pident
    if (align_len > max_align[gene]) max_align[gene] = align_len

    interval = contig ":" start "-" end
    if (intervals[gene] !~ interval) {
      intervals[gene] = intervals[gene] ";" interval
    }
  }
}

END {
  print "gene", "contig", "coordinate_span", "num_CDS_intervals", "num_100pct_maleZero_intervals", "all_CDS_intervals_100pct_maleZero", "max_pident", "max_align_len", "intervals"

  for (gene in seen) {
    if (total[gene] == zero_count[gene]) {
      all_zero = "yes"
    } else {
      all_zero = "no"
    }

    coord = min_start[gene] "-" max_end[gene]

    print gene, contigs[gene], coord, total[gene], zero_count[gene], all_zero, max_pident[gene], max_align[gene], intervals[gene]
  }
}
EOF

# 4. Run gene-level summary.
# Sort by the all-zero column and gene name.

awk -f "$AWK_SCRIPT" "$ROI_ZERO" "$ROI_CDS" \
| sort -k6,6 -k1,1 \
> "$OUT"

echo "Wrote:"
echo "$ROI_CDS"
echo "$ROI_ZERO"
echo "$OUT"
