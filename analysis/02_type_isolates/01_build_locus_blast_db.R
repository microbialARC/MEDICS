# Stage 02.1 - Build a BLAST database of the candidate loci.
#
# For each candidate from stage 01, cuts the reference sequence spanning the
# locus and both flanks, and writes all of them to one FASTA. The flanks are
# included because they are the conserved regions a PCR primer pair would anneal
# to, so a genome only counts as typeable when the whole amplicon is present.
#
# The FASTA is the BLAST subject: each isolate assembly is queried against it in
# stage 02.2, which is the reverse of the usual arrangement and is what lets one
# database serve every isolate.
#
# Run from the repository root:
#   Rscript analysis/02_type_isolates/01_build_locus_blast_db.R <output_dir>

suppressPackageStartupMessages({
  library(Biostrings)
  library(openxlsx2)
})

source(Sys.getenv("MEDICS_CONFIG", unset = "analysis/config.R"))

output_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(output_dir)) stop("Usage: Rscript 01_build_locus_blast_db.R <output_dir>")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

candidate_df <- openxlsx2::read_xlsx(file.path(TABLE_DIR, "locus_candidates.xlsx"))
reference_genome <- Biostrings::readDNAStringSet(reference_path("concat.fasta"))
reference_contig_ids <- sub("\\s.*$", "", names(reference_genome))

locus_seq_list <- vector("list", nrow(candidate_df))

for (candidate_idx in seq_len(nrow(candidate_df))) {

  contig_name <- candidate_df$contig[candidate_idx]
  contig_index <- match(contig_name, reference_contig_ids)
  if (is.na(contig_index)) stop("Contig not found in reference: ", contig_name)

  # Locus plus both flanks, clamped to the contig.
  contig_length <- Biostrings::width(reference_genome)[contig_index]
  region_start <- max(1L, candidate_df$left_flank_start[candidate_idx])
  region_end   <- min(contig_length, candidate_df$right_flank_end[candidate_idx])

  region_seq <- Biostrings::DNAStringSet(
    Biostrings::subseq(reference_genome[[contig_index]], start = region_start, end = region_end)
  )
  # Candidates are named by row order in the stage 01 table, which is how they
  # are referred to throughout stages 02-04 (cand_13 is the benchmark locus).
  names(region_seq) <- paste0("cand_", candidate_idx)
  locus_seq_list[[candidate_idx]] <- region_seq
}

locus_db_fasta <- file.path(output_dir, "all_nicu_candidate_seqs.fasta")
Biostrings::writeXStringSet(do.call(c, locus_seq_list), filepath = locus_db_fasta)

message("Wrote ", nrow(candidate_df), " candidate loci to ", locus_db_fasta)
message("Now run: makeblastdb -in ", basename(locus_db_fasta), " -dbtype nucl -title all_nicu_candidates")
