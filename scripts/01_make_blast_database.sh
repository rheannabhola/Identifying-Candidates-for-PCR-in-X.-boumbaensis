#!/bin/bash
#SBATCH --job-name=make_boumb_blastdb
#SBATCH --account=rrg-ben
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=make_boumb_blastdb_%j.out
#SBATCH --error=make_boumb_blastdb_%j.err

module load StdEnv/2023
module load blast+

# Purpose:
# Build a BLAST nucleotide database from the female X. boumbaensis reference genome.
# This database is used to map X. laevis CDS sequences onto the X. boumbaensis assembly.

# Edit this path if running in a new location.
REF="/home/rheanna/projects/rrg-ben/rheanna/boumbaensis_Adam_genome_assembly/boumb.wbubble.fa"

OUTDIR="/home/rheanna/projects/rrg-ben/rheanna/genomewide_coding_coverage"
DB_PREFIX="${OUTDIR}/boumb_wbubble_genome_db"

mkdir -p "$OUTDIR"

makeblastdb \
  -in "$REF" \
  -dbtype nucl \
  -out "$DB_PREFIX"

echo "Finished building BLAST database:"
echo "$DB_PREFIX"
