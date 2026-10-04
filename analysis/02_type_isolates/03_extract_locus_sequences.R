# Stage 02.3 - Cut each locus out of every isolate that carries it.
#
# Reads the BLAST hit tables from stage 02.2, calls a locus present in an isolate
# when its best hit clears both thresholds in config.R (>= 80% identity and
# >= 80% coverage - the "80_80" naming in the working directories), and writes
# one FASTA per locus holding that locus' sequence from every isolate carrying
# it. Those per-locus FASTAs are what stage 02.4 aligns.
#
# Isolates whose hit falls below either threshold are left untyped at that locus
# rather than forced into a haplotype; this is why the benchmark counts 1,060
# typed isolates out of 1,594 in the window.
#
# Expects the per-locus presence table assembled from stage 02.2:
#   blastn_output_sum_df  one row per locus x isolate, with columns
#                         query (locus), genome, contig, start, end, presence
#   local_genome_assembly_df  genome -> assembly path
#
# Run as a SLURM job; parallel workers come from --cpus-per-task.

suppressPackageStartupMessages(library(parallel))

source(Sys.getenv("MEDICS_CONFIG", unset = "analysis/config.R"))

n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

# Load Biostrings once here so every forked worker inherits it instead of loading it itself
invisible(loadNamespace("Biostrings"))

# Extract one locus' hit region from every genome where it is present and write
# them to a single FASTA. Returns TRUE so failed loci can be identified after
# the parallel run.
extract_locus_seq <- function(locus_entry) {

  export_file <- file.path(extract_locus_seq_dir, paste0(locus_entry, "_extract_seq.fasta"))

  if (file.exists(export_file)) {
    message(locus_entry, " found. Skip")
    return(TRUE)
  }

  message(locus_entry, " not found. Generating")
  locus_blastn_res <- blastn_output_sum_df[blastn_output_sum_df$query == locus_entry &
                                             blastn_output_sum_df$presence, ]

  locus_query_seqlist <- vector("list", nrow(locus_blastn_res))

  for (query_idx in seq_len(nrow(locus_blastn_res))) {

    query_name <- locus_blastn_res$genome[query_idx]
    query_genome <- Biostrings::readDNAStringSet(
      local_genome_assembly_df$path[local_genome_assembly_df$genome == query_name]
    )

    hit_contig <- locus_blastn_res$contig[query_idx]
    hit_start  <- locus_blastn_res$start[query_idx]
    hit_end    <- locus_blastn_res$end[query_idx]

    query_contig_ids <- sub("\\s.*$", "", names(query_genome))
    hit_contig_idx <- match(hit_contig, query_contig_ids)

    hit_seqset <- Biostrings::DNAStringSet(
      Biostrings::subseq(query_genome[[hit_contig_idx]], start = hit_start, end = hit_end)
    )
    names(hit_seqset) <- paste0(query_name, "_", locus_entry)
    locus_query_seqlist[[query_idx]] <- hit_seqset
  }

  # Write to a temp file and rename it, so a job killed mid-write cannot leave
  # a partial FASTA that the file.exists() check would skip on the next run
  tmp_file <- paste0(export_file, ".tmp")
  Biostrings::writeXStringSet(do.call(c, locus_query_seqlist), filepath = tmp_file)
  file.rename(tmp_file, export_file)

  TRUE
}

# Run loci in parallel. mclapply forks worker processes, which see the data
# frames and paths above without exporting them. Errors are caught per locus, so
# one failure does not stop the remaining loci.
locus_entries <- unique(blastn_output_sum_df$query)
extract_status <- parallel::mclapply(locus_entries, function(locus_entry) {
  tryCatch(extract_locus_seq(locus_entry), error = function(e) {
    message(locus_entry, " failed: ", conditionMessage(e))
    FALSE
  })
}, mc.cores = n_threads, mc.preschedule = FALSE)

# A worker that was killed (e.g. out of memory) returns NULL rather than FALSE,
# so check for TRUE explicitly
failed_loci <- locus_entries[!vapply(extract_status, isTRUE, logical(1))]
if (length(failed_loci) > 0) {
  stop("Sequence extraction failed for ", length(failed_loci), " locus/loci: ",
       paste(failed_loci, collapse = ", "))
}
