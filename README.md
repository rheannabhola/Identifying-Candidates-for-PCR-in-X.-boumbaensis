# Xenopus boumbaensis Candidate Gene Identification Workflow

This repository contains scripts used to identify candidate PCR marker genes for putative sex-linked regions in Xenopus boumbaensis.

## Project goal

The goal of this workflow is to identify coding-like regions in the female X. boumbaensis reference genome that show 100% zero coverage in male X. boumbaensis whole-genome sequencing data.

These regions may represent female-specific, W-linked, or highly sex-diverged candidate marker loci.

## Biological rationale

The X. boumbaensis reference genome used here is from a female assembly. If female-specific or W-linked sequences are present in this reference, male WGS reads mapped to the female reference may show no coverage across those regions.

Candidate logic:

    female reference contains sequence
    + sequence is coding/CDS-like based on X. laevis BLAST annotation
    + male WGS mapped to female reference has zero coverage
    + region occurs on a candidate sex-linked contig
    = potential PCR marker candidate

These candidates should be interpreted as candidate marker regions, not confirmed sex-determining genes.

## Main data inputs

Large input files are not included in this repository.

Key input data used in the workflow:

- Female X. boumbaensis reference genome: boumb.wbubble.fa
- Male X. boumbaensis WGS BAM mapped to the female reference:
  SRR35972782_male_vs_femaleBoumb.sorted.bam
- X. laevis CDS FASTA file
- Candidate sex-linked contigs from previous analyses:
  - tig00003177
  - tig00004591

## Main tools used

- BLASTN
- samtools
- bedtools
- awk
- R

## Workflow overview

1. Build a BLAST database from the female X. boumbaensis reference genome.
2. BLAST X. laevis CDS sequences against the X. boumbaensis reference.
3. Convert CDS BLAST hits into BED-style coding-like intervals.
4. Use samtools depth to calculate male WGS coverage across coding-like intervals.
5. Extract coding-like intervals with 100% zero male WGS coverage.
6. Intersect male-zero intervals with CDS BLAST annotations.
7. Filter for high-confidence hits:
   - alignment length >= 300 bp
   - percent identity >= 90%
   - male WGS zero-depth fraction = 1.0
8. Prioritize candidates on:
   - tig00003177
   - tig00004591 SNP-associated block
9. Remove generic LOC-only annotations for a cleaner candidate list.
10. Summarize genes where all detected CDS-like intervals have 100% zero male coverage.

## Final output

The final gene-level output identifies non-LOC genes in the regions of interest where all detected high-quality CDS-like intervals have 100% zero coverage in male X. boumbaensis WGS.

Example final table columns:

- Gene
- Contig
- Coordinate span
- Number of CDS-like intervals
- Number of 100% male-zero CDS-like intervals
- Whether all CDS-like intervals are 100% male-zero
- Max identity percentage
- Max alignment length

## Interpretation

The final candidates are coding-like, male-zero regions on candidate sex-linked contigs. They are promising PCR marker candidates, but require experimental validation across known male and female individuals.

## Notes and limitations

- Candidate regions are based on X. laevis CDS BLAST annotations, not a finalized X. boumbaensis gene annotation.
- Therefore, these are described as CDS-like or exon-like intervals rather than confirmed X. boumbaensis exons.
- A gene-level “yes” result means that all detected high-quality CDS-like intervals for that gene in the region of interest had 100% zero male WGS coverage.
- This does not prove the gene is sex-determining; it only identifies strong PCR marker candidates.

## Repository structure

    scripts/
        Shell and R scripts used for candidate discovery and table generation.

    docs/
        Workflow notes and documentation.

    example_outputs/
        Small example output tables. Large intermediate files are excluded.
