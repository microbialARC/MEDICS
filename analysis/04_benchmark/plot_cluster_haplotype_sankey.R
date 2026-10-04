# Stage 04 - Sankey of WGS transmission clusters against locus haplotype types.
#
# One figure per locus combination selected in stage 03. The left axis is the
# epidemiological truth (THRESHER transmission clusters, plus one pooled
# Non-cluster block); the right axis is the haplotype types the loci assign. A
# type is drawn in its cluster's color when it is credited to that cluster, so
# the figure shows at a glance how much of the cluster structure the locus
# recovers and where it splits or merges clusters.
#
# Needs the external datasets (see docs/data-availability.md):
#   MEDICS_LOCUS_HAPLOTYPES      stage 02.5 haplotype calls
#   MEDICS_THRESHER_SNP_SUMMARY  WGS strain and cluster assignment
#   MEDICS_LOCUS_SELECTION_DIR   stage 03 output (best_combo_df.RDS)
#
# Run from the repository root:
#   Rscript analysis/04_benchmark/plot_cluster_haplotype_sankey.R
#
# Outputs
#   results/figures/cluster_haplotype_sankey_<tag>.pdf
#   per-combination typing detail, written to the stage 03 working directory
#   rather than into the repository because it is isolate-level

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggalluvial)
  library(colorspace)
})

source(Sys.getenv("MEDICS_CONFIG", unset = "analysis/config.R"))

# TRUE: the Haplotype side shows only the colored types (those meeting Condition 1, 2 or 3),
# in the middle. Genomes of the other types are hidden in two blocks above and below them, so
# the Transmission Cluster side keeps its full sizes. FALSE: every type is shown
show_matched_type_only <- TRUE

# TRUE: restrict to the benchmark window, so the figure covers exactly the
# isolates behind the numbers in benchmark_locus_vs_wgs.R. FALSE: use the whole
# collection, which draws more isolates but no longer matches those numbers.
restrict_to_window <- TRUE

# By default every combination selected in stage 03 is drawn. Set
# MEDICS_SANKEY_LOCI to draw a specific locus or combination instead, as it is
# named in the stage 01 table - "cand_13" for the benchmark locus, or
# "cand_1 + cand_13" for a pair. This is how the committed locus 13 figure was
# produced, since the stage 03 objective selects cand_1 as its best single locus.
sankey_loci_override <- Sys.getenv("MEDICS_SANKEY_LOCI", unset = "")

# Condition 3 thresholds, as shares of each unit's typed genomes
own_cluster_min_frac <- CLUSTER_MATCH_MIN_OWN_FRACTION
other_unit_max_frac <- CLUSTER_MATCH_MAX_OTHER_FRACTION

## Inputs ----
selection_dir <- external_path("locus_selection_dir")
snp_summary   <- readRDS(external_path("thresher_snp_summary"))
strain_compo  <- snp_summary[["genomes"]]
cluster_compo <- snp_summary[["clusters"]]

cand_haplotype <- readRDS(external_path("locus_haplotypes"))
names(cand_haplotype) <- vapply(
  cand_haplotype, function(cand_entry) cand_entry$seq_haplotype_list[[1]]$candidate, character(1)
)

if (nzchar(sankey_loci_override)) {
  best_combo_df <- data.frame(
    n_cand = length(strsplit(sankey_loci_override, split = " + ", fixed = TRUE)[[1]]),
    tie_idx = 1L,
    candidates = sankey_loci_override,
    stringsAsFactors = FALSE
  )
} else {
  best_combo_df <- readRDS(file.path(selection_dir, "best_combo_df.RDS"))
}

cluster2strain_df <- do.call(rbind, lapply(cluster_compo, function(cluster_entry) {
  data.frame(cluster = cluster_entry$cluster, strain = cluster_entry$strain)
}))

# group of each genome: its THRESHER cluster, or NC_<strain> when its strain is in no cluster
strain_compo$cluster <- vapply(strain_compo$strain, function(strain_entry) {
  cluster_id <- unique(cluster2strain_df$cluster[cluster2strain_df$strain == strain_entry])
  if (length(cluster_id)) cluster_id else paste0("NC_", strain_entry)
}, character(1))

strain_compo$collection_date <- as.Date(strain_compo$collection_date)
if (restrict_to_window) {
  strain_compo <- strain_compo[
    !is.na(strain_compo$collection_date) &
      strain_compo$collection_date >= BENCHMARK_START &
      strain_compo$collection_date <= BENCHMARK_END, ]
}

genome_cluster <- setNames(strain_compo$cluster, strain_compo$genome)
# TRUE for THRESHER clusters, FALSE for NC_ strains; thresher_cluster lists the clusters only
total_cluster <- sort(unique(strain_compo$cluster))
thresher_cluster <- total_cluster[total_cluster %in% cluster2strain_df$cluster]

for(row_idx in seq_len(nrow(best_combo_df))){

  ## Haplotypes of each candidate in the combination ----
  row_cand <- strsplit(best_combo_df$candidates[row_idx], split = " + ", fixed = TRUE)[[1]]
  # interchangeable candidates are listed as "cand_a|cand_b"; the first one is used
  row_cand <- sapply(strsplit(row_cand, split = "|", fixed = TRUE), function(member_entry) member_entry[1])

  cand_hap_list <- vector("list", length(row_cand))
  names(cand_hap_list) <- row_cand

  for(cand_idx in seq_along(row_cand)){
    cand_entry <- row_cand[cand_idx]
    cand_hap_list[[cand_idx]] <- lapply(cand_haplotype[[cand_entry]]$seq_haplotype_list, function(haplotype_entry){
      haplotype_entry$haplotype_genome[haplotype_entry$haplotype_genome %in% strain_compo$genome]
    })
    names(cand_hap_list[[cand_idx]]) <- sapply(cand_haplotype[[cand_entry]]$seq_haplotype_list, function(haplotype_entry) haplotype_entry$haplotype_idx)
  }

  # genomes typed at one or more candidates of the combination
  total_cand_genomes <- unique(unlist(cand_hap_list, use.names = FALSE))

  ## Multi-locus type of each genome ----
  # one row per genome with its haplotype at each candidate:
  # NA when the genome is untyped there, "a/b" when it carries more than one haplotype
  locus_haplo_genome_df <- do.call(rbind, lapply(total_cand_genomes, function(genome_entry){

    genome_cand_hap <- vapply(cand_hap_list, function(cand_entry){
      genome_hap <- names(cand_entry)[vapply(cand_entry, function(hap_genomes) genome_entry %in% hap_genomes, logical(1))]
      if(length(genome_hap) == 0) NA_character_ else paste(sort(genome_hap), collapse = "/")
    }, character(1))

    data.frame(genome = genome_entry, as.list(genome_cand_hap), check.names = FALSE)
  }))

  # the type gives each candidate's haplotype by locus position in the combination,
  # e.g. "L1H3-L2H1-L3H12" for cand_1 + cand_13 + cand_75 ("L2HNA" when untyped there)
  locus_haplo_genome_df$type <- apply(locus_haplo_genome_df[, row_cand, drop = FALSE], 1, function(genome_hap){
    paste0("L", seq_along(genome_hap), "H", genome_hap, collapse = "-")
  })

  n_genome_untyped <- sum(!complete.cases(locus_haplo_genome_df[, row_cand, drop = FALSE]))
  if(n_genome_untyped > 0){
    message(n_genome_untyped, " genomes are untyped at one or more of ", paste(row_cand, collapse = " + "), " (NA in their type)")
  }

  locus_isolate_df <- data.frame(
    genome = locus_haplo_genome_df$genome,
    cluster = unname(genome_cluster[locus_haplo_genome_df$genome]),
    type = locus_haplo_genome_df$type
  )
  locus_isolate_df$is_cluster <- locus_isolate_df$cluster %in% thresher_cluster
  # NC_ strains are lumped into one Non-cluster unit, for display and for the matching below
  locus_isolate_df$cluster_display <- ifelse(locus_isolate_df$is_cluster, locus_isolate_df$cluster, "Non-cluster")

  ## Define haplo_match_name ----
  # A type takes the name of the THRESHER cluster it matches. A unit is a THRESHER cluster
  # or Non-cluster (all NC_ strains together), and its genomes are those typed at one or more candidates.
  # Condition 1: Haplotype fully reproduces the unit
  # Condition 2: Haplotype contains purely part of the unit (no other unit contained)
  # Condition 3: Haplotype contains >= 80% of the unit's genomes and <= 20% of the genomes of each other unit
  type_genome_list <- split(locus_isolate_df$genome, locus_isolate_df$type)

  # share of each unit's typed genomes that falls in each type (type x unit)
  type_unit_n <- unclass(table(locus_isolate_df$type, locus_isolate_df$cluster_display))
  type_unit_frac <- sweep(type_unit_n, 2, colSums(type_unit_n), "/")

  haplo_match_name <- setNames(rep(NA_character_, length(type_genome_list)), names(type_genome_list))
  haplo_match_condition <- setNames(rep(NA_integer_, length(type_genome_list)), names(type_genome_list))

  for(type_entry in names(type_genome_list)){

    # units present in this type, and the THRESHER cluster with the largest share among them
    unit_frac <- setNames(type_unit_frac[type_entry, ], colnames(type_unit_frac))
    unit_frac <- unit_frac[unit_frac > 0]
    cluster_frac <- unit_frac[names(unit_frac) %in% thresher_cluster]
    if(length(cluster_frac) == 0) next

    top_cluster <- names(cluster_frac)[which.max(cluster_frac)]
    other_unit_frac <- unit_frac[names(unit_frac) != top_cluster]

    if(length(other_unit_frac) == 0){
      # Conditions 1 and 2: every genome of the type belongs to the cluster
      haplo_match_name[type_entry] <- top_cluster
      haplo_match_condition[type_entry] <- if(cluster_frac[top_cluster] == 1) 1L else 2L
    }else if(cluster_frac[top_cluster] >= own_cluster_min_frac && all(other_unit_frac <= other_unit_max_frac)){
      # Condition 3
      haplo_match_name[type_entry] <- top_cluster
      haplo_match_condition[type_entry] <- 3L
    }
  }

  # clusters that a type reproduces exactly
  cluster_exact_name <- unique(haplo_match_name[which(haplo_match_condition == 1L)])

  ## Build Sankey plot data ----
  # both axes ranked by genome count, most genomes at the top; Non-cluster is ranked
  # with the clusters, and ties keep their previous order
  cluster_present <- thresher_cluster[thresher_cluster %in% locus_isolate_df$cluster_display]
  cluster_size <- table(locus_isolate_df$cluster_display)
  cluster_levels <- c(cluster_present, "Non-cluster")
  cluster_levels <- cluster_levels[order(-cluster_size[cluster_levels])]

  # with show_matched_type_only, genomes of types matching no cluster are hidden in two equal
  # blocks, one above and one below the colored types, so those sit in the middle of the
  # Haplotype side; hidden genomes of the higher-ranked clusters fill the top block
  hidden_type <- c(top = "Unmatched_top", bottom = "Unmatched_bottom")
  type_display <- locus_isolate_df$type
  if(show_matched_type_only){
    hidden_genome_idx <- which(is.na(haplo_match_name[type_display]))
    hidden_genome_idx <- hidden_genome_idx[order(match(locus_isolate_df$cluster_display[hidden_genome_idx], cluster_levels))]
    type_display[hidden_genome_idx] <- hidden_type["bottom"]
    type_display[hidden_genome_idx[seq_len(floor(length(hidden_genome_idx) / 2))]] <- hidden_type["top"]
  }

  # one flow per cluster x type cell, width = isolates, in ggalluvial lodes form:
  # two rows per flow, one for each axis
  flow_df <- as.data.frame(table(
    cluster_display = locus_isolate_df$cluster_display,
    type = type_display
  ), stringsAsFactors = FALSE)
  colnames(flow_df)[3] <- "n_isolate"
  flow_df <- flow_df[flow_df$n_isolate > 0, ]
  flow_df$alluvium <- seq_len(nrow(flow_df))
  # flows into the hidden blocks are fully transparent
  flow_df$flow_alpha <- ifelse(flow_df$type %in% hidden_type, 0, 0.65)

  sankey_plot_df <- rbind(
    data.frame(alluvium = flow_df$alluvium, axis = "cluster", stratum = flow_df$cluster_display,
               cluster_display = flow_df$cluster_display, n_isolate = flow_df$n_isolate, flow_alpha = flow_df$flow_alpha),
    data.frame(alluvium = flow_df$alluvium, axis = "type", stratum = flow_df$type,
               cluster_display = flow_df$cluster_display, n_isolate = flow_df$n_isolate, flow_alpha = flow_df$flow_alpha)
  )
  # the hidden blocks get no outline either
  sankey_plot_df$stratum_border <- ifelse(sankey_plot_df$stratum %in% hidden_type, "transparent", "grey")

  type_size <- table(type_display)
  type_levels <- names(type_size)[order(-type_size)]
  # the hidden blocks go above and below the colored types
  type_levels <- c(
    intersect(hidden_type["top"], type_levels),
    setdiff(type_levels, hidden_type),
    intersect(hidden_type["bottom"], type_levels)
  )
  sankey_plot_df$axis <- factor(sankey_plot_df$axis, levels = c("cluster", "type"))
  # ggalluvial draws the first level at the top, so no rev() here
  sankey_plot_df$stratum <- factor(sankey_plot_df$stratum, levels = c(cluster_levels, type_levels))
  sankey_plot_df$cluster_display <- factor(sankey_plot_df$cluster_display, levels = cluster_levels)
  rownames(sankey_plot_df) <- NULL

  # only the clusters are labelled; the Haplotype strata get no text (matches still show as colors)
  stratum_label <- c(
    setNames(cluster_levels, cluster_levels),
    setNames(rep("", length(type_levels)), type_levels)
  )

  ## Colors ----
  # one scale covers both axes: clusters keep their color, the background and the types
  # stay neutral, and a type matching a cluster takes that cluster's color.
  # Colors go by cluster order, not rank, so the ranking doesn't reshuffle them
  cluster_color <- qualitative_hcl(length(cluster_present), palette = "Dark 3")
  names(cluster_color) <- cluster_present
  cluster_color["Non-cluster"] <- "#E3E0DC"
  stratum_color <- c(cluster_color, setNames(rep("#E3E0DC", length(type_levels)), type_levels))
  stratum_color[names(haplo_match_name)[!is.na(haplo_match_name)]] <- cluster_color[haplo_match_name[!is.na(haplo_match_name)]]
  stratum_color[intersect(hidden_type, type_levels)] <- "transparent"

  ## Plot ----
  plot <- ggplot(sankey_plot_df,
                 aes(x = axis,
                     y = n_isolate,
                     stratum = stratum,
                     alluvium = alluvium)) +
    geom_flow(aes(fill = cluster_display, alpha = flow_alpha),
              color = "transparent") +
    geom_stratum(aes(fill = stratum, color = "transparent"),
                 linewidth = 0.0001,
                 width = 0.35) +
    geom_text(stat = "stratum",
              aes(label = stratum_label[as.character(after_stat(stratum))]),
              size = 2.5) +
    scale_fill_manual(values = stratum_color, guide = "none") +
    # flow_alpha and stratum_border hold the values to draw with
    scale_alpha_identity() +
    scale_color_identity() +
    scale_x_discrete(labels = c(cluster = "Transmission Cluster", type = "Haplotype"),
                     expand = c(0.15, 0.15)) +
    labs(y = "Genomes",
         title = paste0(ifelse(length(row_cand) == 1, "Candidate ", "Candidates "),
                        paste(gsub("cand_", "", row_cand), collapse = " + "))) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_text(face = "bold",
                                     angle = 0,
                                     size = 27.5),
          axis.text.y = element_text(face = "bold",
                                     size = 25),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold", size = 30),
          plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
          plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(linewidth = 1,fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_blank(),
          axis.ticks.length = unit(0.3, "cm"))

  # file names carry n, the tie index and the candidates, e.g. ..._n3_tie1_1_13_75
  combo_file_tag <- paste0("n", best_combo_df$n_cand[row_idx], "_tie", best_combo_df$tie_idx[row_idx], "_",
                           paste(gsub("cand_", "", row_cand), collapse = "_"))

  # the matched-only version gets its own file, so it doesn't overwrite the full plot
  ggsave(
    filename = file.path(FIGURE_DIR, paste0("cluster_haplotype_sankey_", combo_file_tag,
                                            if(show_matched_type_only) "_matched_only", ".pdf")),
    plot = plot,
    width = 15,
    height = 25,
    limitsize = FALSE
  )

  ## Save ----
  # Isolate-level typing detail stays outside the repository (see docs/data-availability.md)
  combo_rds_dir <- file.path(selection_dir, "RDS")
  if(!dir.exists(combo_rds_dir)){
    dir.create(combo_rds_dir, recursive = TRUE)
  }
  saveRDS(
    list(
      row_cand = row_cand,
      locus_haplo_genome_df = locus_haplo_genome_df,
      locus_isolate_df = locus_isolate_df,
      cluster_exact_name = cluster_exact_name,
      haplo_match_name = haplo_match_name,
      haplo_match_condition = haplo_match_condition
    ),
    file.path(combo_rds_dir, paste0("Cluster_combo_", combo_file_tag, ".RDS"))
  )
}
