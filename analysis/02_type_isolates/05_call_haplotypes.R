# Stage 02.5 - Call haplotypes from the per-locus alignments.
#
# A haplotype is one distinct aligned sequence at one locus, and an isolate's
# haplotype is the sequence it carries there. This is the whole typing scheme:
# no distance threshold and no clustering, just exact sequence identity over the
# aligned locus, which is what makes a haplotype call reproducible between
# laboratories reading the same amplicon.
#
# Writes candidate_aln_haplotype.RDS - a list with one entry per locus, each
# holding the alignment and its haplotypes with the genomes carrying them. That
# file is the input to stages 03 and 04 (MEDICS_LOCUS_HAPLOTYPES in config.R).
#
# Run from the repository root:
#   Rscript analysis/02_type_isolates/05_call_haplotypes.R <alignment_dir> <output_dir>

suppressPackageStartupMessages({
  library(Biostrings)
  library(parallel)
})

command_args <- commandArgs(trailingOnly = TRUE)
alignment_dir <- command_args[1]
output_dir    <- command_args[2]
if (is.na(alignment_dir) || is.na(output_dir)) {
  stop("Usage: Rscript 05_call_haplotypes.R <alignment_dir> <output_dir>")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

alignment_files <- list.files(path = alignment_dir, pattern = "*aln.fasta", full.names = TRUE)
if (length(alignment_files) == 0) stop("No *_aln.fasta files in ", alignment_dir)

aln_haplotype <- mclapply(alignment_files, function(alignment_file) {

  locus_name <- gsub("\\_aln.fasta$", "", basename(alignment_file))
  aln_seq <- Biostrings::readDNAMultipleAlignment(alignment_file)
  # Identical aligned sequences collapse into one haplotype
  seq_unique_haplotype <- unique(as.character(aln_seq))

  seq_haplotype_list <- vector("list", length(seq_unique_haplotype))

  for (haplotype_idx in seq_along(seq_unique_haplotype)) {

    haplotype_seq <- seq_unique_haplotype[haplotype_idx]
    # Alignment names are "<genome>_<locus>"; strip the locus to recover the genome
    haplotype_genome_entry <- gsub(paste0("_", locus_name), "",
                                   names(which(as.character(aln_seq) == haplotype_seq)))

    seq_haplotype_list[[haplotype_idx]] <- list(
      candidate       = locus_name,
      haplotype_idx   = haplotype_idx,
      haplotype_seq   = haplotype_seq,
      haplotype_genome = haplotype_genome_entry
    )
  }

  list(
    candidate = locus_name,
    aln_seq = aln_seq,
    seq_haplotype_list = seq_haplotype_list
  )
}, mc.cores = detectCores())

saveRDS(aln_haplotype, file.path(output_dir, "candidate_aln_haplotype.RDS"))

message("Called haplotypes for ", length(aln_haplotype), " loci -> ",
        file.path(output_dir, "candidate_aln_haplotype.RDS"))
