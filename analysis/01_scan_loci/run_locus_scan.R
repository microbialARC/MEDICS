# Stage 01 - Scan the reference genome for candidate barcoding loci.
#
# Reads the entropy profile in data/reference/ and enumerates every window that
# concentrates informative sites while keeping its flanks free of them, so the
# flanks can serve as conserved primer-binding sites. Writes the candidate table
# and one figure per candidate.
#
# Runs on the files distributed with this repository - no external data needed.
# Run from the repository root:
#
#   Rscript analysis/01_scan_loci/run_locus_scan.R
#
# Outputs
#   results/tables/locus_candidates.xlsx
#   results/figures/locus_candidate_site_distribution.pdf
#   results/figures/locus_candidates/locus_<i>.pdf
#
# The run behind the proposal produced 78 candidates, committed under
# results/. Re-running overwrites them in place.

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(Biostrings)
  library(IRanges)
  library(ggplot2)
  library(patchwork)
  library(openxlsx2)
})

source(Sys.getenv("MEDICS_CONFIG", unset = "analysis/config.R"))
source(repo_path("analysis", "01_scan_loci", "locus_scan_functions.R"))

## Inputs ----
site_entropy  <- readRDS(reference_path("entropy.RDS"))
new_pos       <- readRDS(reference_path("new_pos.RDS"))
mges_df       <- readRDS(reference_path("mges.RDS"))
chr_bins      <- readRDS(reference_path("chr_bins.RDS"))
concat_genome <- Biostrings::readDNAStringSet(reference_path("concat.fasta"))

candidate_dir <- file.path(FIGURE_DIR, "locus_candidates")

## Scan ----
# Settings as used for the grant analysis. min_informative_sites = 3 keeps the
# search permissive; candidates are ranked afterwards, and n_top = "all" keeps
# every one of them for the downstream typing stages.
locus_candidates <- entropy_barcoding(
  site_entropy          = site_entropy,
  new_pos               = new_pos,
  mges_df               = mges_df,
  concat_genome         = concat_genome,
  chr_bins              = chr_bins,
  min_locus_size        = 100,
  max_locus_size        = 1000,
  flank_size            = 100,
  sliding_step          = 50,
  min_informative_sites = 3,
  n_top                 = "all",
  output_dir            = candidate_dir
)

## Collect outputs ----
# entropy_barcoding() writes the per-locus figures, the workbook and the summary
# figure into one directory; the repository keeps tables and figures apart, so
# the two non-figure outputs are moved up.
collected_outputs <- c(
  "locus_candidates.xlsx"                 = file.path(TABLE_DIR, "locus_candidates.xlsx"),
  "locus_candidate_site_distribution.pdf" = file.path(FIGURE_DIR, "locus_candidate_site_distribution.pdf")
)
for (base_name in names(collected_outputs)) {
  from_path <- file.path(candidate_dir, base_name)
  if (file.exists(from_path)) file.rename(from_path, collected_outputs[[base_name]])
}

message("Candidate loci found: ", nrow(locus_candidates))
