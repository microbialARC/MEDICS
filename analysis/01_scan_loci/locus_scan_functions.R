# Helper Functions -----
## Core detection and parallel driver ----

detect_available_cores <- function() {
  slurm_allocated_cores <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (nzchar(slurm_allocated_cores)) {
    return(max(1L, as.integer(slurm_allocated_cores)))
  }
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}

## Parallel  ----
# If windows no parallel used
run_tasks_in_parallel <- function(task_list,
                                  worker_function,
                                  n_cores) {
  
  if (length(task_list) == 0) return(list())
  
  if (n_cores <= 1L || .Platform$OS.type == "windows") {
    return(lapply(task_list, worker_function))
  }
  
  task_results <- parallel::mclapply(
    task_list,
    worker_function,
    mc.cores       = min(n_cores, length(task_list)),
    mc.preschedule = FALSE
  )
  
  worker_failed <- vapply(task_results, inherits, logical(1), "try-error")
  if (any(worker_failed)) {
    stop("Parallel worker failed on ", sum(worker_failed), " task(s). ",
         "First error: ", as.character(task_results[[which(worker_failed)[1]]]))
  }
  
  task_results
}

## Prefix-sum lookup tables ----

build_entropy_lookups <- function(site_entropy,
                                  genome_length,
                                  mge_ranges) {
  
  entropy_vector <- rep(0, genome_length)
  entropy_vector[site_entropy$position] <- site_entropy$entropy
  
  is_informative_site <- entropy_vector > 0
  
  # Every informative site in ascending genomic order. A window's site set is
  # always a contiguous index range into these two vectors, which is what lets
  # a window be identified by a single (first, last) index pair.
  informative_positions <- which(is_informative_site)
  informative_entropy_values <- entropy_vector[informative_positions]
  
  # coverage() is much cheaper than expanding MGE ranges to integer positions.
  mge_coverage <- IRanges::coverage(mge_ranges, width = genome_length)
  non_mge_mask <- as.vector(mge_coverage) == 0
  
  # Element k + 1 holds the total over positions 1..k, so any closed interval
  # [a, b] is table[b + 1] - table[a] in constant time. The degenerate empty
  # flank (b == a - 1) falls out of the same arithmetic as 0.
  list(
    entropy_vector = entropy_vector,
    informative_positions = informative_positions,
    informative_entropy_values = informative_entropy_values,
    cumulative_site_count = c(0L, cumsum(is_informative_site)),
    cumulative_entropy = c(0,  cumsum(entropy_vector)),
    cumulative_nonmge_site_count = c(0L, cumsum(is_informative_site & non_mge_mask))
  )
}

## Task list ----

build_scan_tasks <- function(contig_boundaries,
                             mge_ranges,
                             lookups,
                             min_locus_size,
                             min_informative_sites) {
  
  cumulative_site_count <- lookups$cumulative_site_count
  # Create a list of tasks, each of which is a non-MGE interval on a contig that is large
  # for subsequent parallel scanning. Each task is a list with the contig name, contig start and end, and the interval start and end.
  scan_tasks <- list()
  
  for (contig_idx in seq_len(nrow(contig_boundaries))) {
    
    contig_start <- contig_boundaries$new_start[contig_idx]
    contig_end   <- contig_boundaries$new_end[contig_idx]
    contig_name  <- contig_boundaries$ctg[contig_idx]
    
    contig_range <- IRanges::IRanges(start = contig_start, end = contig_end)
    non_mge_ranges_in_contig <- IRanges::setdiff(contig_range, mge_ranges)
    
    for (interval_idx in seq_along(non_mge_ranges_in_contig)) {
      
      interval_start <- IRanges::start(non_mge_ranges_in_contig)[interval_idx]
      interval_end <- IRanges::end(non_mge_ranges_in_contig)[interval_idx]
      # First check 
      # If interval smaller than locus size
      if (interval_end - interval_start + 1 < min_locus_size) next
      # Second check
      # An interval that cannot possibly hold enough informative sites never
      # reaches a worker.
      sites_in_interval <- cumulative_site_count[interval_end + 1] -
        cumulative_site_count[interval_start]
      if (sites_in_interval < min_informative_sites) next
      
      scan_tasks[[length(scan_tasks) + 1]] <- list(
        contig_name = contig_name,
        contig_start = contig_start,
        contig_end = contig_end,
        interval_idx = interval_idx,
        interval_start = interval_start,
        interval_end = interval_end
      )
    }
  }
  scan_tasks
}

## De-duplication ----

collapse_duplicate_site_sets <- function(window_df,
                                         n_informative_total) {
  
  if (is.null(window_df) || nrow(window_df) == 0) return(window_df)
  
  site_set_key <- as.numeric(window_df$first_site_idx) * (n_informative_total + 1) +
    window_df$last_site_idx
  
  window_df[!duplicated(site_set_key), , drop = FALSE]
}

## Interval scan (one parallel task) ----
#
# Returns a data frame of clean candidate windows for one non-MGE interval,
# or NULL if the interval yields none.

scan_interval_windows <- function(task,
                                  lookups,
                                  min_locus_size,
                                  max_locus_size,
                                  flank_size,
                                  sliding_step,
                                  min_informative_sites,
                                  n_informative_total) {
  
  cumulative_site_count <- lookups$cumulative_site_count
  cumulative_entropy <- lookups$cumulative_entropy
  cumulative_nonmge_site_count <- lookups$cumulative_nonmge_site_count
  informative_positions <- lookups$informative_positions
  
  interval_length <- task$interval_end - task$interval_start + 1
  effective_max_size <- min(max_locus_size, interval_length)
  candidate_window_sizes <- seq(min_locus_size, effective_max_size, by = sliding_step)
  
  clean_windows_per_size <- vector("list", length(candidate_window_sizes))
  
  ## PLACEHOLDER (dirty layer): re-enable alongside clean_windows_per_size.
  # dirty_windows_per_size <- vector("list", length(candidate_window_sizes))
  
  for (size_idx in seq_along(candidate_window_sizes)) {
    
    current_window_size <- candidate_window_sizes[size_idx]
    
    # Every grid start for this window size, evaluated at once.
    window_starts <- seq(task$interval_start,
                         task$interval_end - current_window_size + 1,
                         by = sliding_step)
    if (length(window_starts) == 0) next
    window_ends <- window_starts + current_window_size - 1
    
    # O(1) per window: the min_informative_sites filter costs two lookups.
    internal_site_counts <- cumulative_site_count[window_ends + 1] -
      cumulative_site_count[window_starts]
    
    passing_windows <- which(internal_site_counts >= min_informative_sites)
    if (length(passing_windows) == 0) next
    
    window_starts <- window_starts[passing_windows]
    window_ends <- window_ends[passing_windows]
    internal_site_counts <- internal_site_counts[passing_windows]
    
    left_flank_starts <- pmax(task$contig_start, window_starts - flank_size)
    left_flank_ends <- window_starts - 1
    right_flank_starts <- window_ends + 1
    right_flank_ends <- pmin(task$contig_end, window_ends + flank_size)
    
    # Replaces IRanges::intersect() + extract_informative_sites(): the non-MGE
    # mask is already baked into the prefix sum, and only the count is needed.
    flank_site_counts <-
      (cumulative_nonmge_site_count[left_flank_ends + 1] -
         cumulative_nonmge_site_count[left_flank_starts]) +
      (cumulative_nonmge_site_count[right_flank_ends + 1] -
         cumulative_nonmge_site_count[right_flank_starts])
    
    # Clean-only: discard everything with informative flanks before any row is
    # built. Dirty windows are the large majority, so this is where most of the
    # remaining per-window cost disappears.
    clean_windows <- which(flank_site_counts == 0)
    if (length(clean_windows) == 0) next
    
    window_starts <- window_starts[clean_windows]
    window_ends <- window_ends[clean_windows]
    internal_site_counts <- internal_site_counts[clean_windows]
    left_flank_starts <- left_flank_starts[clean_windows]
    left_flank_ends <- left_flank_ends[clean_windows]
    right_flank_starts <- right_flank_starts[clean_windows]
    right_flank_ends <- right_flank_ends[clean_windows]
    
    ## PLACEHOLDER (dirty layer): capture the complement before the subset
    ## above overwrites the vectors, then build a dirty frame the same way the
    ## clean one is built below, with flank_info_count taken from
    ## flank_site_counts[dirty_windows] and stored in dirty_windows_per_size.
    # dirty_windows <- which(flank_site_counts > 0)
    
    # Index range into informative_positions, by binary search. No per-window
    # subsetting of entropy_vector.
    first_site_idx <- findInterval(window_starts - 1L, informative_positions) + 1L
    last_site_idx  <- findInterval(window_ends,        informative_positions)
    
    total_entropy <- cumulative_entropy[window_ends + 1] -
      cumulative_entropy[window_starts]
    
    clean_windows_per_size[[size_idx]] <- data.frame(
      contig  = task$contig_name,
      interval_idx = task$interval_idx,
      window_start = window_starts,
      window_end = window_ends,
      window_size = current_window_size,
      total_entropy = total_entropy,
      n_informative_sites = as.integer(internal_site_counts),
      site_density = internal_site_counts / current_window_size,
      first_site_idx = first_site_idx,
      last_site_idx = last_site_idx,
      flank_info_count = 0L,
      left_flank_start = left_flank_starts,
      left_flank_end = left_flank_ends,
      right_flank_start = right_flank_starts,
      right_flank_end = right_flank_ends,
      stringsAsFactors  = FALSE
    )
  }
  
  clean_windows_in_interval <- dplyr::bind_rows(clean_windows_per_size)
  if (nrow(clean_windows_in_interval) == 0) return(NULL)
  
  # Collapsing inside the worker keeps what crosses the fork boundary small.
  collapse_duplicate_site_sets(clean_windows_in_interval, n_informative_total)
  
  ## PLACEHOLDER (dirty layer): return both layers instead of the clean frame,
  ## and update the bind step in entropy_barcoding() to match.
  # list(
  #   clean = collapse_duplicate_site_sets(clean_windows_in_interval, n_informative_total),
  #   dirty = collapse_duplicate_site_sets(dplyr::bind_rows(dirty_windows_per_size),
  #                                        n_informative_total)
  # )
}

## Materialise the string columns for surviving rows only  ----

finalise_candidate_table <- function(candidate_df, lookups, layer_label = "clean") {
  
  if (is.null(candidate_df) || nrow(candidate_df) == 0) return(candidate_df)
  
  informative_positions      <- lookups$informative_positions
  informative_entropy_values <- lookups$informative_entropy_values
  
  site_index_ranges <- Map(`:`, candidate_df$first_site_idx, candidate_df$last_site_idx)
  
  candidate_df$info_site <- vapply(
    site_index_ranges,
    function(site_index_range) paste(informative_positions[site_index_range], collapse = ", "),
    character(1)
  )
  candidate_df$site_entropy <- vapply(
    site_index_ranges,
    function(site_index_range) paste(informative_entropy_values[site_index_range], collapse = ", "),
    character(1)
  )
  
  candidate_df$first_info_site    <- informative_positions[candidate_df$first_site_idx]
  candidate_df$last_info_site     <- informative_positions[candidate_df$last_site_idx]
  candidate_df$contrast_score     <- NA_real_
  candidate_df$search_layer_label <- layer_label
  
  candidate_df[, c(
    "contig", "interval_idx", "window_start", "window_end", "window_size",
    "total_entropy", "n_informative_sites", "site_density",
    "info_site", "site_entropy", "first_info_site", "last_info_site",
    "flank_info_count",
    "left_flank_start", "left_flank_end", "right_flank_start", "right_flank_end",
    "contrast_score", "search_layer_label",
    "first_site_idx", "last_site_idx"
  ), drop = FALSE]
}

## Scoring and plotting (one parallel task per candidate)  ----


#  Arrow polygons for the gene track
build_gene_arrow_polygons <- function(coding_bins,
                                      track_start,
                                      track_end,
                                      body_half_height = 0.2,
                                      head_half_height = 0.3,
                                      head_length_fraction = 0.03) {
  
  head_length_bp <- head_length_fraction * (track_end - track_start)
  
  polygon_per_gene <- lapply(seq_len(nrow(coding_bins)), function(gene_idx) {
    
    gene_strand <- coding_bins$strand[gene_idx]
    visible_start <- max(coding_bins$start[gene_idx], track_start)
    visible_end <- min(coding_bins$end[gene_idx], track_end)
    
    head_in_view <- (gene_strand == "+" && coding_bins$end[gene_idx] <= track_end) ||
      (gene_strand == "-" && coding_bins$start[gene_idx] >= track_start)
    
    if (!head_in_view) {
      return(data.frame(
        gene_idx = gene_idx,
        x = c(visible_start, visible_end, visible_end, visible_start),
        y = c(-body_half_height, -body_half_height, body_half_height, body_half_height)
      ))
    }
    
    # A gene shorter than the head is drawn as a bare triangle, as in gggenes.
    head_length <- min(head_length_bp, visible_end - visible_start)
    
    if (gene_strand == "+") {
      tail_x <- visible_start
      tip_x <- visible_end
      head_base_x <- visible_end - head_length
    } else {
      tail_x <- visible_end
      tip_x <- visible_start
      head_base_x <- visible_start + head_length
    }
    
    # Traced from the tail: lower body edge, lower barb, tip, upper barb,
    # upper body edge. The same order works for both strands.
    data.frame(
      gene_idx = gene_idx,
      x = c(tail_x, head_base_x, head_base_x, tip_x, head_base_x, head_base_x, tail_x),
      y = c(-body_half_height, -body_half_height, -head_half_height, 0,
            head_half_height, body_half_height, body_half_height)
    )
  })
  
  dplyr::bind_rows(polygon_per_gene)
}


score_and_plot_candidate <- function(row_idx,
                                     candidate_df,
                                     chr_bins,
                                     lookups,
                                     output_dir) {
  
  site_index_range    <- candidate_df$first_site_idx[row_idx]:candidate_df$last_site_idx[row_idx]
  site_positions      <- lookups$informative_positions[site_index_range]
  site_entropy_values <- lookups$informative_entropy_values[site_index_range]
  total_entropy       <- sum(site_entropy_values)
  
  window_start     <- candidate_df$window_start[row_idx]
  window_end       <- candidate_df$window_end[row_idx]
  left_flank_start <- candidate_df$left_flank_start[row_idx]
  right_flank_end  <- candidate_df$right_flank_end[row_idx]
  
  kde_bandwidth <- tryCatch(
    as.numeric(stats::bw.SJ(site_positions)),
    error = function(e) as.numeric(stats::bw.nrd0(site_positions))
  )
  
  kde_curve <- stats::density(
    x       = site_positions,
    weights = site_entropy_values / total_entropy,
    bw      = kde_bandwidth,
    from    = window_start,
    to      = window_end,
    n       = 512
  )
  
  peak_kde_density <- as.numeric(max(kde_curve$y))
  mean_kde_density <- as.numeric(mean(kde_curve$y))
  peakedness_ratio <- as.numeric(peak_kde_density / mean_kde_density)
  composite_score  <- as.numeric(total_entropy * peakedness_ratio *
                                   candidate_df$site_density[row_idx])
  
  kde_df <- data.frame(
    position     = kde_curve$x,
    density      = kde_curve$y,
    region_label = "Locus"
  )
  
  entropy_stem_df <- data.frame(
    position     = site_positions,
    entropy      = site_entropy_values,
    region_label = "Locus"
  )
  
  region_shading_df <- data.frame(
    xmin         = c(left_flank_start, window_start, window_end + 1),
    xmax         = c(window_start - 1, window_end,   right_flank_end),
    region_label = c("Left Flank", "Locus", "Right Flank")
  )
  
  region_fill_colours <- c("Left Flank" = "#4DBBD5", "Locus" = "#E64B35", "Right Flank" = "#4DBBD5")
  site_colours        <- c("Locus" = "#E64B35", "Flank" = "#4DBBD5")
  
  # Panel 1: entropy stems
  entropy_panel <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = region_shading_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region_label),
      alpha = 0.12
    ) +
    ggplot2::scale_fill_manual(values = region_fill_colours, name = "Region") +
    ggplot2::geom_segment(
      data = entropy_stem_df,
      ggplot2::aes(x = position, xend = position, y = 0, yend = entropy, colour = region_label),
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = entropy_stem_df,
      ggplot2::aes(x = position, y = entropy, colour = region_label),
      size = 1.2
    ) +
    ggplot2::scale_colour_manual(values = site_colours, name = "Site") +
    ggplot2::labs(x = NULL, y = "Site Entropy") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position    = "none",
      axis.text.x        = ggplot2::element_blank(),
      axis.title.x       = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(0, 10, 0, 10)
    )
  
  # Panel 2: weighted KDE
  kde_panel <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = region_shading_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region_label),
      alpha = 0.12
    ) +
    ggplot2::scale_fill_manual(values = region_fill_colours, name = "Region") +
    ggplot2::geom_line(data = kde_df, ggplot2::aes(x = position, y = density),
                       linewidth = 0.8, colour = "grey20") +
    ggplot2::geom_area(
      data = kde_df[kde_df$region_label == "Locus", ],
      ggplot2::aes(x = position, y = density), fill = "#E64B35", alpha = 0.3
    ) +
    #ggplot2::annotate(
    #  "text",
    #  x = (window_start + window_end) / 2,
    #  y = max(kde_df$density) * 0.93,
    #  label = paste0("BW = ", round(kde_bandwidth, 1), " bp"),
    #  size = 3.5, fontface = "italic"
    #) +
    ggplot2::labs(x = NULL, y = "Weighted KDE Density") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position    = "none",
      axis.text.x        = ggplot2::element_blank(),
      axis.title.x       = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(0, 10, 0, 10)
    )
  
  # Panel 3: gene annotation track
  visible_bins <- chr_bins[chr_bins$end >= left_flank_start & chr_bins$start <= right_flank_end, ]
  coding_bins  <- visible_bins[visible_bins$gene != "non-cds", ]
  
  if (nrow(coding_bins) > 0) {
    
    coding_bins$plot_start <- pmax(coding_bins$start, left_flank_start)
    coding_bins$plot_end <- pmin(coding_bins$end, right_flank_end)
    coding_bins$gene_mid <- (coding_bins$plot_start + coding_bins$plot_end) / 2
    
    gene_arrow_df <- build_gene_arrow_polygons(coding_bins, left_flank_start, right_flank_end)
    
    gene_panel <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = region_shading_df,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region_label),
        alpha = 0.12
      ) +
      ggplot2::scale_fill_manual(values = region_fill_colours, guide = "none") +
      # Backbone goes down before the arrows so it runs between genes, not through them.
      ggplot2::annotate(
        "segment",
        x = left_flank_start, xend = right_flank_end, y = 0, yend = 0,
        color = "grey60", linewidth = 0.3
      ) +
      ggplot2::geom_polygon(
        data = gene_arrow_df,
        ggplot2::aes(x = x, y = y, group = gene_idx),
        fill = "#7CAE00", color = "grey30", linewidth = 0.3, linejoin = "mitre"
      ) +
      ggplot2::geom_text(
        data = coding_bins,
        ggplot2::aes(x = gene_mid, y = 0.55, label = gene),
        size = 2.3, fontface = "italic", color = "grey20", check_overlap = TRUE
      ) +
      ggplot2::coord_cartesian(xlim = c(left_flank_start, right_flank_end),
                               ylim = c(-0.6, 0.9)) +
      ggplot2::labs(x = "Contig Position (bp)") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 10, 5, 10)
      )
    
  } else {
    
    gene_panel <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = region_shading_df,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = region_label),
        alpha = 0.12
      ) +
      ggplot2::scale_fill_manual(values = region_fill_colours, guide = "none") +
      # Same backbone as the gene branch, so tracks with and without genes match.
      ggplot2::annotate(
        "segment",
        x = left_flank_start, xend = right_flank_end, y = 0, yend = 0,
        color = "grey60", linewidth = 0.3
      ) +
      ggplot2::annotate(
        "text",
        x = (left_flank_start + right_flank_end) / 2, y = 0.55,
        label = "No coding genes in this region",
        size = 3.5, fontface = "italic", color = "grey50"
      ) +
      ggplot2::coord_cartesian(xlim = c(left_flank_start, right_flank_end),
                               ylim = c(-0.6, 0.9)) +
      ggplot2::labs(x = "Contig Position (bp)") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 10, 5, 10)
      )
  }
  
  plot_title <- paste0(" Size: ", candidate_df$window_size[row_idx], " bp",
                       " | Informative Sites: ", candidate_df$n_informative_sites[row_idx],
                       " | Bandwidth: ", round(kde_bandwidth, 1), " bp",
                       " | Site Density: ", round(candidate_df$site_density[row_idx], 4),
                       " | Peakedness Ratio: ", round(peakedness_ratio, 2),
                       " | Composite Score: ", round(composite_score, 2))
  
  combined_plot <- entropy_panel / kde_panel / gene_panel +
    patchwork::plot_layout(heights = c(1, 1, 0.5)) +
    patchwork::plot_annotation(title = plot_title,
                                theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11.5)))
  
  # ggsave() rather than pdf()/dev.off(): opening a device inside a forked
  # child is the usual source of truncated PDFs.
  ggplot2::ggsave(
    filename = file.path(output_dir, sprintf("locus_%02d.pdf", row_idx)),
    plot     = combined_plot,
    width    = 10,
    height   = 6,
    device   = "pdf"
  )
  
  list(peakedness_ratio = peakedness_ratio, composite_score = composite_score)
}


plot_candidate_site_distribution <- function(candidate_df, layer_label = "clean", output_file = NULL) {
  
  if (is.null(candidate_df) || nrow(candidate_df) == 0) {
    message("No ", layer_label, " candidates to plot")
    return(invisible(NULL))
  }
  
  # One row per informative-site count: how many candidates have that count,
  # and the median window size among those candidates
  site_count_summary <- candidate_df %>%
    dplyr::group_by(n_informative_sites) %>%
    dplyr::summarise(
      n_candidates = dplyr::n(),
      median_window_size = median(window_size),
      .groups = "drop"
    )
  
  # Both axes are counts, so tick marks are restricted to whole numbers
  integer_breaks <- function(axis_limits) unique(round(pretty(axis_limits)))
  
  site_distribution_plot <- ggplot2::ggplot(
    site_count_summary,
    ggplot2::aes(x = n_informative_sites, y = n_candidates, fill = median_window_size)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_x_continuous(breaks = integer_breaks, expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_y_continuous(breaks = integer_breaks, expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_fill_viridis_c(option = "viridis", name = "Median window size (bp)") +
    ggplot2::labs(
      x = "Informative sites",
      y = "Number of candidates"
    ) +
    ggplot2::theme(axis.text.x = element_text(size = 20, lineheight = 0.9),
                   axis.text.y = element_text(size = 20),
                   axis.title.x = element_text(size = 25),
                   axis.title.y = element_text(size = 25),
                   plot.background = element_blank(),
                   panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(linewidth = 1),
                   legend.text = element_text(size = 15),
                   legend.title = element_text(size = 17.5, hjust = 0.5),
                   axis.line = element_blank(),
                   legend.position = "top",
                   legend.title.position = "top",
                   legend.key.width = unit(1, "null"),
                   legend.margin = margin(0, 0, 0, 0))
  
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = site_distribution_plot, width = 7.5, height = 5)
  }
  
  invisible(site_distribution_plot)
}


## Resolve top n
# Accepts a positive number, or the string "all" to keep every candidate. ---
resolve_n_top <- function(n_top, n_candidates) {
  
  if (is.character(n_top)) {
    if (!identical(tolower(trimws(n_top)), "all")) {
      stop('n_top must be a positive number or the string "all", not "', n_top, '".')
    }
    return(n_candidates)
  }
  
  if (!is.numeric(n_top) || length(n_top) != 1L || is.na(n_top) || n_top < 1) {
    stop('n_top must be a positive number or the string "all".')
  }
  
  min(as.integer(n_top), n_candidates)
}

# Main Function ----

entropy_barcoding <- function(site_entropy,
                              new_pos,
                              mges_df,
                              concat_genome,
                              chr_bins,
                              min_locus_size = 100,
                              max_locus_size = 1000,
                              flank_size = 100,
                              sliding_step = 50,
                              min_informative_sites = 5,
                              n_top = 10,
                              n_cores = detect_available_cores(),
                              output_dir = ".") {
  # Libraries
  library(dplyr)
  library(parallel)
  library(Biostrings)
  library(IRanges)
  library(ggplot2)
  library(patchwork)
  library(openxlsx2)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  genome_length <- Biostrings::nchar(concat_genome)
  
  mge_ranges <- IRanges::IRanges(start = mges_df$start, end = mges_df$end)
  
  lookups <- build_entropy_lookups(site_entropy, genome_length, mge_ranges)
  n_informative_total <- length(lookups$informative_positions)
  
  contig_boundaries <- new_pos[["ctg_pos_df"]]
  contig_boundaries <- contig_boundaries[contig_boundaries$length > min_locus_size, ]
  
  scan_tasks <- build_scan_tasks(contig_boundaries,
                                 mge_ranges,
                                 lookups,
                                 min_locus_size,
                                 min_informative_sites)
  
  message("Scanning ", length(scan_tasks), " non-MGE interval(s) across ",
          format(n_informative_total, big.mark = ","), " informative sites")
  
  scan_worker <- function(task) {
    scan_interval_windows(
      task = task,
      lookups = lookups,
      min_locus_size = min_locus_size,
      max_locus_size = max_locus_size,
      flank_size = flank_size,
      sliding_step = sliding_step,
      min_informative_sites = min_informative_sites,
      n_informative_total = n_informative_total
    )
  }
  
  scan_results <- run_tasks_in_parallel(scan_tasks, scan_worker, n_cores)
  
  # Workers return in task order, so enumeration order - and therefore tie
  # breaking in the global collapse below - is preserved.
  clean_candidates <- dplyr::bind_rows(scan_results)
  clean_candidates <- collapse_duplicate_site_sets(clean_candidates, n_informative_total)
  
  ## PLACEHOLDER (dirty layer): when scan_interval_windows() returns both
  ## layers again, replace the two lines above with the pair below.
  # clean_candidates     <- dplyr::bind_rows(lapply(scan_results, `[[`, "clean"))
  # dirty_candidates_raw <- dplyr::bind_rows(lapply(scan_results, `[[`, "dirty"))
  # clean_candidates     <- collapse_duplicate_site_sets(clean_candidates, n_informative_total)
  # dirty_candidates_raw <- collapse_duplicate_site_sets(dirty_candidates_raw, n_informative_total)
  
  if (nrow(clean_candidates) == 0) {
    
    message("No clean candidates found.")
    
    ## PLACEHOLDER (dirty layer): second search layer, run only when the clean
    ## layer comes up empty. contrast_score is computed here.
    # if (nrow(dirty_candidates_raw) > 0) {
    #   dirty_candidates_raw <- finalise_candidate_table(dirty_candidates_raw, lookups, "dirty")
    #   dirty_candidates_raw <- dirty_candidates_raw %>%
    #     dplyr::slice_max(site_density, n = n_top, with_ties = TRUE)
    #   # ... contrast_score calculation ...
    #   openxlsx2::write_xlsx(dirty_candidates_raw,
    #                         file = file.path(output_dir, "dirty_candidates.xlsx"))
    #   return(invisible(dirty_candidates_raw))
    # }
    
    return(invisible(clean_candidates))
  }
  
  n_top_resolved <- resolve_n_top(n_top, nrow(clean_candidates))
  
  clean_candidates <- clean_candidates %>%
    dplyr::slice_max(n_informative_sites, n = n_top_resolved, with_ties = TRUE)
  
  message("")
  message("In total ",nrow(clean_candidates), " clean candidates are found")
  clean_candidates <- finalise_candidate_table(clean_candidates, lookups, "clean")
  
  plot_results <- run_tasks_in_parallel(
    as.list(seq_len(nrow(clean_candidates))),
    function(row_idx) {
      score_and_plot_candidate(row_idx, clean_candidates, chr_bins, lookups, output_dir)
    },
    n_cores
  )
  
  clean_candidates$peakedness_ratio <- vapply(plot_results, `[[`, numeric(1), "peakedness_ratio")
  clean_candidates$composite_score  <- vapply(plot_results, `[[`, numeric(1), "composite_score")
  
  clean_candidates <- clean_candidates[
    , setdiff(names(clean_candidates), c("first_site_idx", "last_site_idx")), drop = FALSE]
  
  # Summarize the distribution of clean candidates
  plot_candidate_site_distribution(
    clean_candidates,
    layer_label = "clean",
    output_file = file.path(output_dir, "locus_candidate_site_distribution.pdf")
  )

  openxlsx2::write_xlsx(clean_candidates,
                        file = file.path(output_dir, "locus_candidates.xlsx"))
  
  invisible(clean_candidates)
}