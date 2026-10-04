# Stage 04 - Benchmark the selected locus against WGS and against MLST.
#
# Answers three questions about the benchmark locus (cand_13) on the CHOP NICU
# S. aureus collection, and writes the aggregate tables the grant quotes:
#
#   1. How many isolates does the locus type, into how many haplotype types, and
#      how many WGS-defined transmission clusters do those types reach?
#   2. How faithfully does a single haplotype correspond to a cluster?
#   3. Does the locus resolve genomes below the MLST sequence-type level?
#
# Needs the external datasets (see docs/data-availability.md):
#   MEDICS_LOCUS_HAPLOTYPES      stage 02.5 haplotype calls
#   MEDICS_THRESHER_SNP_SUMMARY  WGS strain and cluster assignment
#   MEDICS_MLST                  per-genome ST / CC
#
# Run from the repository root:
#   Rscript analysis/04_benchmark/benchmark_locus_vs_wgs.R
#
# Outputs (all aggregate - no row identifies an isolate):
#   results/tables/benchmark_summary.csv
#   results/tables/per_st_haplotype_counts.csv
#   results/tables/haplotype_cluster_match.csv
#   results/tables/monthly_typing_capacity.csv
#   results/figures/locus13_monthly_typing_capacity.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source(Sys.getenv("MEDICS_CONFIG", unset = "analysis/config.R"))

# The locus this benchmark reports on, as named by its row in the stage 01 table.
BENCHMARK_LOCUS <- "cand_13"

## Inputs ----
snp_summary   <- readRDS(external_path("thresher_snp_summary"))
strain_compo  <- snp_summary[["genomes"]]
cluster_compo <- snp_summary[["clusters"]]
mlst_df       <- readRDS(external_path("mlst"))[["genomes_mlst_cc"]]

cand_haplotype <- readRDS(external_path("locus_haplotypes"))
names(cand_haplotype) <- vapply(
  cand_haplotype, function(cand_entry) cand_entry$seq_haplotype_list[[1]]$candidate, character(1)
)

cluster2strain_df <- do.call(rbind, lapply(cluster_compo, function(cluster_entry) {
  data.frame(cluster = cluster_entry$cluster, strain = cluster_entry$strain)
}))

# Each genome's epidemiological unit: its transmission cluster, or NC_<strain>
# when its strain was never placed in one. Non-cluster strains are kept because
# a haplotype that bleeds into them is not cluster-specific.
strain_compo$cluster <- vapply(strain_compo$strain, function(strain_entry) {
  cluster_id <- unique(cluster2strain_df$cluster[cluster2strain_df$strain == strain_entry])
  if (length(cluster_id)) cluster_id else paste0("NC_", strain_entry)
}, character(1))

## Benchmark window ----
# Sampling before Sep 2022 and after Mar 2025 is sparse and uneven, so the
# window is fixed in config.R and every number below is computed inside it.
strain_compo$collection_date <- as.Date(strain_compo$collection_date)
window_compo <- strain_compo[
  !is.na(strain_compo$collection_date) &
    strain_compo$collection_date >= BENCHMARK_START &
    strain_compo$collection_date <= BENCHMARK_END, ]

window_genome  <- window_compo$genome
window_cluster <- unique(cluster2strain_df$cluster[cluster2strain_df$strain %in% window_compo$strain])
window_years   <- as.numeric(BENCHMARK_END - BENCHMARK_START) / 365

## Typing at the benchmark locus ----
# One row per typed genome. A genome is typed when it carries a haplotype at the
# locus, i.e. its amplicon cleared the stage 02.3 presence thresholds. Type names
# follow the multi-locus convention ("L1H3" = haplotype 3 at the first locus), so
# the single-locus case reads the same way as a combination; a genome matching
# more than one haplotype gets them joined with "/".
locus_haplotype_list <- cand_haplotype[[BENCHMARK_LOCUS]]$seq_haplotype_list
names(locus_haplotype_list) <- vapply(
  locus_haplotype_list, function(hap_entry) as.character(hap_entry$haplotype_idx), character(1)
)

typed_genome <- unique(unlist(
  lapply(locus_haplotype_list, function(hap_entry) hap_entry$haplotype_genome), use.names = FALSE
))

genome_type <- vapply(typed_genome, function(genome_entry) {
  genome_hap <- names(locus_haplotype_list)[vapply(
    locus_haplotype_list,
    function(hap_entry) genome_entry %in% hap_entry$haplotype_genome,
    logical(1)
  )]
  paste0("L1H", paste(sort(genome_hap), collapse = "/"))
}, character(1))

genome_cluster <- setNames(strain_compo$cluster, strain_compo$genome)

typing_df <- data.frame(
  genome  = typed_genome,
  type    = unname(genome_type),
  cluster = unname(genome_cluster[typed_genome]),
  stringsAsFactors = FALSE
)
typing_df$is_cluster <- typing_df$cluster %in% cluster2strain_df$cluster
# Every non-cluster strain is pooled into one unit for the matching below
typing_df$unit <- ifelse(typing_df$is_cluster, typing_df$cluster, "Non-cluster")

window_typing_df <- typing_df[typing_df$genome %in% window_genome, ]

n_typed_isolate  <- nrow(window_typing_df)
n_haplotype_type <- length(unique(window_typing_df$type))
# Clusters the locus reaches: those with at least one typed isolate in the window.
# This is coverage of the cluster set, not a claim that a haplotype reproduces a
# cluster one-to-one - that stricter question is the matching section below.
n_cluster_reached <- length(unique(window_typing_df$cluster[window_typing_df$is_cluster]))
pct_cluster_reached <- 100 * n_cluster_reached / length(window_cluster)

## Matching haplotypes to clusters ----
# How faithfully a type corresponds to a cluster, on the typed genomes in the
# window. A type is credited to the cluster holding the largest share of it when:
#   condition 1  every genome of the type is in that cluster, and the type holds
#                the whole cluster          (the type reproduces the cluster)
#   condition 2  every genome of the type is in that cluster, but the type holds
#                only part of it            (a pure subset)
#   condition 3  the type holds >= 80% of that cluster's typed genomes and
#                <= 20% of every other unit's
type_unit_n <- unclass(table(window_typing_df$type, window_typing_df$unit))
type_unit_frac <- sweep(type_unit_n, 2, colSums(type_unit_n), "/")

match_name      <- setNames(rep(NA_character_, nrow(type_unit_n)), rownames(type_unit_n))
match_condition <- setNames(rep(NA_integer_,   nrow(type_unit_n)), rownames(type_unit_n))

for (type_entry in rownames(type_unit_n)) {

  unit_frac <- type_unit_frac[type_entry, ]
  unit_frac <- unit_frac[unit_frac > 0]
  cluster_frac <- unit_frac[names(unit_frac) != "Non-cluster"]
  if (length(cluster_frac) == 0) next

  top_cluster <- names(cluster_frac)[which.max(cluster_frac)]
  other_unit_frac <- unit_frac[names(unit_frac) != top_cluster]

  if (length(other_unit_frac) == 0) {
    match_name[type_entry] <- top_cluster
    match_condition[type_entry] <- if (cluster_frac[[top_cluster]] == 1) 1L else 2L
  } else if (cluster_frac[[top_cluster]] >= CLUSTER_MATCH_MIN_OWN_FRACTION &&
             all(other_unit_frac <= CLUSTER_MATCH_MAX_OTHER_FRACTION)) {
    match_name[type_entry] <- top_cluster
    match_condition[type_entry] <- 3L
  }
}

match_df <- data.frame(
  condition = c(1L, 2L, 3L),
  description = c(
    "type reproduces the cluster exactly",
    "type is a pure subset of one cluster",
    sprintf("type holds >=%.0f%% of the cluster and <=%.0f%% of every other unit",
            100 * CLUSTER_MATCH_MIN_OWN_FRACTION, 100 * CLUSTER_MATCH_MAX_OTHER_FRACTION)
  ),
  n_types = vapply(c(1L, 2L, 3L), function(k) sum(match_condition == k, na.rm = TRUE), integer(1)),
  n_clusters = vapply(c(1L, 2L, 3L), function(k)
    length(unique(match_name[which(match_condition == k)])), integer(1)),
  stringsAsFactors = FALSE
)

## Resolution below the sequence-type level ----
# An ST-specific haplotype is a type seen in exactly one ST. Counting those per
# ST says how far the locus subdivides each ST.
window_mlst_df <- merge(window_typing_df, mlst_df, by = "genome", all.x = TRUE)
window_st <- setdiff(unique(window_mlst_df$ST), c(NA, "Unassigned"))

type_st_df <- as.data.frame(table(window_mlst_df[, c("type", "ST")])) %>%
  filter(Freq != 0 & ST != "Unassigned")

type_st_summary <- do.call(rbind, lapply(as.character(unique(type_st_df$type)), function(type_entry) {
  type_st <- unique(as.character(type_st_df$ST[type_st_df$type == type_entry]))
  data.frame(type = type_entry, st = paste(type_st, collapse = ", "),
             n_st = length(type_st), is_st_specific = length(type_st) == 1,
             stringsAsFactors = FALSE)
}))

st_specific_counts <- type_st_summary %>%
  filter(is_st_specific) %>%
  count(st, name = "n_st_specific_haplotypes") %>%
  arrange(desc(n_st_specific_haplotypes))

st_count_quantiles <- quantile(st_specific_counts$n_st_specific_haplotypes)

## Monthly typing capacity ----
# Per month: isolates collected, STs and strains seen, and how many of those
# isolates the locus could type. The remainder - isolates the locus cannot type -
# is the sequencing that would still be necessary if the locus were used for
# first-pass surveillance.
window_compo <- window_compo %>%
  mutate(collection_month = as.Date(format(collection_date, "%Y-%m-01")))

monthly_df <- do.call(rbind, lapply(sort(unique(window_compo$collection_month)), function(month_entry) {

  month_genome <- window_compo$genome[window_compo$collection_month == month_entry]
  month_typed  <- intersect(month_genome, typing_df$genome)
  month_strain <- unique(window_compo$strain[window_compo$genome %in% month_genome])

  data.frame(
    month = month_entry,
    n_genome = length(month_genome),
    n_strain = length(month_strain),
    n_ST = length(setdiff(unique(mlst_df$ST[mlst_df$genome %in% month_genome]), c("Unassigned", NA))),
    n_cluster = length(unique(cluster2strain_df$cluster[cluster2strain_df$strain %in% month_strain])),
    n_haplotype = length(unique(typing_df$type[typing_df$genome %in% month_typed])),
    n_typed_genome = length(month_typed),
    stringsAsFactors = FALSE
  )
}))
monthly_df$n_untypeable_genome <- monthly_df$n_genome - monthly_df$n_typed_genome

## Summary table ----
summary_df <- data.frame(
  metric = c(
    "benchmark_window_start", "benchmark_window_end", "benchmark_window_years",
    "collection_isolates", "collection_strains", "collection_clusters",
    "benchmark_locus",
    "isolates_typed", "haplotype_types",
    "clusters_reached", "pct_clusters_reached",
    "types_reproducing_a_cluster", "clusters_reproduced_exactly",
    "sequence_types", "sequence_types_with_specific_haplotypes",
    "st_specific_haplotypes_min", "st_specific_haplotypes_q1",
    "st_specific_haplotypes_median", "st_specific_haplotypes_q3",
    "st_specific_haplotypes_max", "st_specific_haplotypes_mean",
    "st_specific_haplotypes_sd"
  ),
  value = c(
    format(BENCHMARK_START), format(BENCHMARK_END), sprintf("%.2f", window_years),
    nrow(window_compo), length(unique(window_compo$strain)), length(window_cluster),
    BENCHMARK_LOCUS,
    n_typed_isolate, n_haplotype_type,
    n_cluster_reached, sprintf("%.2f", pct_cluster_reached),
    sum(!is.na(match_condition)), sum(match_condition == 1L, na.rm = TRUE),
    length(window_st), nrow(st_specific_counts),
    st_count_quantiles[["0%"]], st_count_quantiles[["25%"]],
    st_count_quantiles[["50%"]], st_count_quantiles[["75%"]],
    st_count_quantiles[["100%"]],
    sprintf("%.2f", mean(st_specific_counts$n_st_specific_haplotypes)),
    sprintf("%.2f", sd(st_specific_counts$n_st_specific_haplotypes))
  ),
  definition = c(
    "first collection date included", "last collection date included",
    "window length in years",
    "isolates collected in the window",
    "WGS-defined strains among them",
    "WGS-defined transmission clusters among them (THRESHER)",
    "locus reported on, by its row in results/tables/locus_candidates.xlsx",
    "isolates in the window carrying a haplotype at the locus",
    "distinct haplotype types among those isolates",
    "clusters with at least one typed isolate",
    "clusters_reached as a percentage of collection_clusters",
    "types credited to a cluster under condition 1, 2 or 3",
    "clusters a single type reproduces exactly (condition 1)",
    "STs among the typed isolates",
    "STs having at least one haplotype seen in no other ST",
    rep("ST-specific haplotypes per ST, across those STs", 7)
  ),
  stringsAsFactors = FALSE
)

## Write outputs ----
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, file.path(TABLE_DIR, "benchmark_summary.csv"), row.names = FALSE)
write.csv(st_specific_counts, file.path(TABLE_DIR, "per_st_haplotype_counts.csv"), row.names = FALSE)
write.csv(match_df, file.path(TABLE_DIR, "haplotype_cluster_match.csv"), row.names = FALSE)
write.csv(monthly_df, file.path(TABLE_DIR, "monthly_typing_capacity.csv"), row.names = FALSE)

## Monthly figure ----
line_colors <- c("Sequenced Genomes"    = "#808285",
                 "Strains"              = "#786452",
                 "Sequence Type"         = "#FFDE21",
                 "Necessary Sequencing"  = "#91a01e")
fill_colors <- c("Unnecessary Sequencing" = "#CC79A7",
                 "Haplotypes"             = "#005587")

monthly_plot <- ggplot(monthly_df, aes(x = month)) +
  # The gap between isolates collected and isolates the locus cannot type is the
  # sequencing the locus would make unnecessary
  geom_ribbon(aes(ymin = n_untypeable_genome, ymax = n_genome, fill = "Unnecessary Sequencing"),
              alpha = 0.35) +
  geom_point(aes(y = n_haplotype, fill = "Haplotypes"), size = 1, color = "#005587") +
  geom_line(aes(y = n_genome, color = "Sequenced Genomes"), linewidth = 1) +
  geom_line(aes(y = n_strain, color = "Strains"), linewidth = 1) +
  geom_line(aes(y = n_ST, color = "Sequence Type"), linewidth = 1) +
  geom_line(aes(y = n_untypeable_genome, color = "Necessary Sequencing"), linewidth = 1) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y", expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = line_colors, breaks = names(line_colors)) +
  scale_fill_manual(values = fill_colors, breaks = names(fill_colors)) +
  labs(x = "Month", y = "Count", color = NULL, fill = NULL) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        axis.line = element_blank(),
        panel.border = element_rect(fill = "transparent", linewidth = 1, color = "black"))

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(FIGURE_DIR, "locus13_monthly_typing_capacity.pdf"),
       monthly_plot, width = 6.5, height = 3.5)

## Report ----
message("Benchmark window: ", BENCHMARK_START, " to ", BENCHMARK_END,
        " (", sprintf("%.2f", window_years), " years)")
message("Collection in window: ", nrow(window_compo), " isolates, ",
        length(window_cluster), " clusters")
message("Locus ", BENCHMARK_LOCUS, ": ", n_typed_isolate, " isolates typed into ",
        n_haplotype_type, " types, reaching ", n_cluster_reached, " clusters (",
        sprintf("%.2f", pct_cluster_reached), "%)")
message("MLST: ", length(window_st), " STs, ", nrow(st_specific_counts),
        " with ST-specific haplotypes (median ",
        st_count_quantiles[["50%"]], ", max ", st_count_quantiles[["100%"]], ")")
