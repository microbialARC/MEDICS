# Shared paths and settings for the MEDICS proof-of-concept analysis.
#
# Every script in analysis/ begins by sourcing this file, so no script carries an
# absolute path. Two kinds of path are resolved here:
#
#   repo_path()     files distributed with this repository (reference profile,
#                   figures, aggregate tables)
#   external_path() the CHOP NICU source datasets, which are NOT distributed
#                   (see docs/data-availability.md). Point the environment
#                   variables below at a local copy to re-run those stages.

## Repository root ----
# Walks up from the working directory until the marker file is found, so scripts
# can be run from anywhere (RStudio, Rscript, an HPC job) without setwd().
find_repo_root <- function(start_dir = getwd()) {
  current_dir <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current_dir, "analysis", "config.R"))) return(current_dir)
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not locate the MEDICS repository root above '", start_dir,
           "'. Set MEDICS_ROOT, or run from inside the repository.")
    }
    current_dir <- parent_dir
  }
}

MEDICS_ROOT <- Sys.getenv("MEDICS_ROOT", unset = NA_character_)
if (is.na(MEDICS_ROOT) || !nzchar(MEDICS_ROOT)) MEDICS_ROOT <- find_repo_root()

repo_path <- function(...) file.path(MEDICS_ROOT, ...)

## Reference genome profile (distributed) ----
# SILOSim profiler output for S. aureus NCTC 8325 (GCF_000013425.1), produced by
# analysis/00_profile_reference.sh. The profiler names every file after the input
# assembly, so the prefix is kept exactly as the tool wrote it.
REFERENCE_PREFIX <- "GCF_000013425.1_ASM1342v1_genomic"

reference_path <- function(suffix) {
  repo_path("data", "reference", paste0(REFERENCE_PREFIX, "_", suffix))
}

## Result directories (distributed) ----
FIGURE_DIR <- repo_path("results", "figures")
TABLE_DIR  <- repo_path("results", "tables")

## External source datasets (not distributed) ----
# Each entry is an environment variable holding a path on the analyst's machine.
EXTERNAL_PATH_VARS <- c(
  # THRESHER SNP summary: isolate -> strain -> transmission cluster.
  # Contains patient-identifying fields and is never copied into this repository.
  thresher_snp_summary = "MEDICS_THRESHER_SNP_SUMMARY",
  # MLST / clonal-complex assignment per genome.
  mlst                 = "MEDICS_MLST",
  # Per-locus haplotype calls from stage 02 (candidate_aln_haplotype.RDS).
  locus_haplotypes     = "MEDICS_LOCUS_HAPLOTYPES",
  # Stage 03 working directory (locus selection intermediates and best_combo_df.RDS).
  locus_selection_dir  = "MEDICS_LOCUS_SELECTION_DIR",
  # Directory of QC-passed genome assemblies, and the table of their paths.
  genome_assembly_dir  = "MEDICS_GENOME_ASSEMBLY_DIR",
  genome_assembly_df   = "MEDICS_GENOME_ASSEMBLY_DF"
)

external_path <- function(key, must_exist = TRUE) {
  if (!key %in% names(EXTERNAL_PATH_VARS)) {
    stop("Unknown external path '", key, "'. Known keys: ",
         paste(names(EXTERNAL_PATH_VARS), collapse = ", "))
  }
  env_var <- EXTERNAL_PATH_VARS[[key]]
  path_value <- Sys.getenv(env_var, unset = "")
  if (!nzchar(path_value)) {
    stop("This stage needs the '", key, "' dataset, which is not distributed with ",
         "the repository.\nSet ", env_var, " to its location on your machine. ",
         "See docs/data-availability.md.", call. = FALSE)
  }
  path_value <- path.expand(path_value)
  if (must_exist && !file.exists(path_value)) {
    stop(env_var, " is set to '", path_value, "', which does not exist.", call. = FALSE)
  }
  path_value
}

## Benchmark window ----
# The grant reports the stretch of the collection with sustained, even sampling.
# Isolates outside it are excluded from every benchmark number.
BENCHMARK_START <- as.Date("2022-09-01")
BENCHMARK_END   <- as.Date("2025-03-01")

## Locus-typing thresholds ----
# A locus counts as present in a genome when its BLAST hit clears both, which is
# why the stage 02/03 working directories are named "80_80".
LOCUS_MIN_IDENTITY <- 80
LOCUS_MIN_COVERAGE <- 80

# A haplotype is said to match a transmission cluster when it holds at least this
# share of the cluster's typed genomes and no more than this share of any other
# unit's genomes
CLUSTER_MATCH_MIN_OWN_FRACTION   <- 0.8
CLUSTER_MATCH_MAX_OTHER_FRACTION <- 0.2
