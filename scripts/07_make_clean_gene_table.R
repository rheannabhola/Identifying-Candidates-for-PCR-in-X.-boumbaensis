# 07_make_clean_gene_table.R

# Purpose:
# Make a clean gene-level candidate table from the gene-level CDS zero-coverage summary.
#
# Input:
# gene_level_CDS_zeroCoverage_summary_noLOC.tsv
#
# Output:
# gene_level_CDS_zeroCoverage_summary_clean.csv
#
# This final table keeps genes where all detected high-quality CDS-like intervals
# in the regions of interest have 100% zero coverage in male X. boumbaensis WGS.

workdir <- "/home/rheanna/projects/rrg-ben/rheanna/genomewide_coding_coverage"

input_file <- file.path(workdir, "gene_level_CDS_zeroCoverage_summary_noLOC.tsv")
out_csv <- file.path(workdir, "gene_level_CDS_zeroCoverage_summary_clean.csv")

# Read without trusting the header, because the shell summary can sometimes be sorted
# in a way that moves the header.
df <- read.delim(
  input_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

colnames(df) <- c(
  "gene",
  "contig",
  "coordinate_span",
  "num_CDS_intervals",
  "num_100pct_maleZero_intervals",
  "all_CDS_intervals_100pct_maleZero",
  "max_pident",
  "max_align_len",
  "intervals"
)

# Remove accidental header row if present
df <- df[df$gene != "gene", ]

# Clean whitespace
df$all_CDS_intervals_100pct_maleZero <- trimws(df$all_CDS_intervals_100pct_maleZero)

# Keep only genes where all detected CDS-like intervals are 100% male-zero
df_yes <- df[df$all_CDS_intervals_100pct_maleZero == "yes", ]

notes_text <- paste(
  "All detected high-quality X. laevis CDS-like intervals for this gene",
  "in the region of interest have 100% zero coverage in male X. boumbaensis WGS",
  "mapped to the female reference."
)

final <- data.frame(
  Gene = df_yes$gene,
  Contig = df_yes$contig,
  `Coordinate span` = df_yes$coordinate_span,
  `Number of CDS-like intervals` = df_yes$num_CDS_intervals,
  `Number of 100% male-zero CDS-like intervals` = df_yes$num_100pct_maleZero_intervals,
  `All CDS-like intervals 100% male-zero?` = df_yes$all_CDS_intervals_100pct_maleZero,
  `Max identity (%)` = df_yes$max_pident,
  `Max alignment length` = df_yes$max_align_len,
  Notes = rep(notes_text, nrow(df_yes)),
  check.names = FALSE
)

# Sort by contig, coordinate start, then gene
get_start <- function(x) as.numeric(sub("-.*", "", x))
final$sort_start <- get_start(final$`Coordinate span`)
final <- final[order(final$Contig, final$sort_start, final$Gene), ]
final$sort_start <- NULL

write.csv(final, out_csv, row.names = FALSE)

message("Wrote CSV file: ", out_csv)
message("Number of genes retained: ", nrow(final))
