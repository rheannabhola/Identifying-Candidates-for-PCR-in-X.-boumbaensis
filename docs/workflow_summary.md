# Workflow Summary: Male-zero CDS-like Candidate Gene Screen

## Purpose

This workflow identifies candidate PCR marker genes in Xenopus boumbaensis by finding CDS-like regions in the female reference genome that show 100% zero coverage in male X. boumbaensis WGS.

The candidates are intended as putative sex-linked marker regions, not confirmed sex-determining genes.

## Biological logic

The X. boumbaensis reference genome used here is female. If female-specific or W-linked sequence is present in the female reference, then male WGS reads mapped to the female reference may show zero coverage across those regions.

The candidate pattern is:

    present in female reference
    coding/CDS-like based on X. laevis annotation
    zero coverage in male X. boumbaensis WGS
    located on candidate sex-linked contigs
    suitable for PCR marker testing

## Main input files

Large input files are not included in the GitHub repository.

Important inputs used in the workflow:

- Female X. boumbaensis reference genome:
  boumb.wbubble.fa

- Male X. boumbaensis WGS BAM mapped to the female reference:
  SRR35972782_male_vs_femaleBoumb.sorted.bam

- X. laevis CDS FASTA:
  XL_CDS_only_nospaces.fasta

- Candidate contigs:
  tig00003177
  tig00004591

## Main tools

- BLASTN
- samtools depth
- bedtools intersect
- bedtools merge
- awk
- R

## Step-by-step summary

### 1. Build BLAST database

A BLAST nucleotide database was made from the female X. boumbaensis reference genome.

This allowed X. laevis CDS sequences to be searched against the X. boumbaensis assembly.

### 2. BLAST X. laevis CDS against X. boumbaensis

X. laevis CDS sequences were BLASTed against the female X. boumbaensis reference.

The BLAST output was used to define CDS-like intervals in the X. boumbaensis assembly and assign putative gene names.

### 3. Convert BLAST hits to BED intervals

BLAST subject coordinates were converted to BED format.

The resulting intervals represented putative coding/CDS-like regions in the X. boumbaensis reference.

### 4. Calculate male WGS depth

Male X. boumbaensis WGS depth was calculated across the genome-wide CDS-like intervals using samtools depth.

The -aa option was used so that zero-coverage positions were retained.

### 5. Identify male-zero coding intervals

Positions with male WGS depth equal to zero were extracted and merged into continuous intervals.

These intervals represent CDS-like regions in the female reference with no coverage from male WGS reads.

### 6. Intersect zero-depth intervals with CDS annotations

Male-zero intervals were intersected back with the X. laevis CDS BLAST hit annotations.

This recovered putative gene names and CDS-hit information for each male-zero region.

### 7. Apply strict filters

Strict candidate filters were applied:

- alignment length >= 300 bp
- percent identity >= 90%
- fraction zero in male WGS = 1.0

This produced a high-confidence set of CDS-like intervals that are 100% zero coverage in the male WGS.

### 8. Focus on regions of interest

Candidates were prioritized on:

- tig00003177
- tig00004591 SNP-associated block from approximately 5.18 Mb to 10.11 Mb

These regions were selected because they were supported by earlier sex-association and coverage analyses.

### 9. Remove LOC-only annotations

Generic LOC-only annotations were removed from the final candidate gene table to produce a more interpretable list of named candidate genes.

LOC genes may still be biologically relevant, but they were excluded from the clean summary table for readability and prioritization.

### 10. Gene-level summary

For each non-LOC gene in the regions of interest, all detected high-quality CDS-like intervals were counted.

The key question was:

    Do all detected CDS-like intervals for this gene have 100% zero male WGS coverage?

Genes where the answer was yes were retained in the clean gene-level candidate table.

## Final interpretation

The final table contains genes where all detected high-quality X. laevis CDS-like intervals in the region of interest had 100% zero coverage in male X. boumbaensis WGS.

These are strong PCR marker candidates, but they are not confirmed W-linked genes or sex-determining genes.

Because the annotation is based on X. laevis CDS BLAST hits, these regions should be described as CDS-like or exon-like intervals rather than confirmed X. boumbaensis exons.

## Recommended next steps

1. Extract candidate sequences from the female X. boumbaensis reference genome.
2. Check primer specificity against the X. boumbaensis reference.
3. Design PCR primers for selected candidate intervals.
4. Test primers on known-sex male and female adult samples.
5. Use validated markers to genotype RNA-seq tadpoles if needed.
