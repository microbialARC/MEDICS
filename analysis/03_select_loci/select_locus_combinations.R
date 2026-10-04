# Stage 03 - Smallest set of loci that tells the transmission clusters apart.
#
# Scores every locus, and every combination of up to three loci, by how many
# pairs of epidemiological units (transmission cluster vs cluster, and cluster vs
# non-cluster strain) their haplotypes separate. Combinations are compared
# exactly by branch and bound rather than greedily, because the best pair of loci
# need not contain the best single locus.
#
# Needs the external datasets (see docs/data-availability.md):
#   MEDICS_LOCUS_HAPLOTYPES    stage 02.5 haplotype calls
#   MEDICS_THRESHER_SNP_SUMMARY  WGS strain and cluster assignment
#   MEDICS_LOCUS_SELECTION_DIR   working directory for intermediates
#
# Intermediates are cached in the working directory, so a re-run skips the
# expensive pairwise steps. Delete them to recompute.
#
# Run from the repository root:
#   Rscript analysis/03_select_loci/select_locus_combinations.R
#
# Output: best_combo_df.RDS in the working directory, read by stage 04.
#
# Note on what this stage selects: by pair separation the best single locus is
# cand_1. The benchmark locus reported in the grant is cand_13, which separates
# slightly fewer pairs but resolves far more haplotypes (677 vs 166 across the
# collection), i.e. it is chosen for resolution below the cluster level rather
# than for this objective. See docs/results.md.

library(dplyr)
library(data.table)
library(parallel)

source(Sys.getenv("MEDICS_CONFIG", unset = "analysis/config.R"))

## Inputs and settings ----
working_dir <- external_path("locus_selection_dir")
setwd(working_dir)

cand_haplotype <- readRDS(external_path("locus_haplotypes"))
snp_summary <- readRDS(external_path("thresher_snp_summary"))
strain_compo <- snp_summary[["genomes"]]
cluster_compo <- snp_summary[["clusters"]]
cluster2strain_df <- do.call(rbind,lapply(cluster_compo,function(cluster_entry){
  
  data.frame(
    cluster = cluster_entry$cluster,
    strain = cluster_entry$strain
  )
  
}))

# group of each genome: its transmission cluster, or NC_<strain> when its strain is in no cluster
strain_compo$cluster <- sapply(strain_compo$strain,function(strain_entry){
  
  cluster_id <- unique(cluster2strain_df$cluster[cluster2strain_df$strain == strain_entry])
  if(length(cluster_id)){
    cluster_id
  }else{
    paste0("NC_", strain_entry)
  }
  
})

cand_names <- sapply(cand_haplotype, function(cand_entry) cand_entry$seq_haplotype_list[[1]]$candidate)
names(cand_haplotype) <- cand_names

total_cluster <- sort(unique(strain_compo$cluster))
genome_cluster <- setNames(strain_compo$cluster, strain_compo$genome)
# TRUE for THRESHER clusters, FALSE for NC_ strains; thresher_cluster lists the clusters only
total_is_cluster <- total_cluster %in% cluster2strain_df$cluster
thresher_cluster <- total_cluster[total_is_cluster]
message(length(thresher_cluster), " clusters and ", sum(!total_is_cluster), " non-cluster strains (NC_ groups)")

# genome names must match between the two inputs; genomes without a strain are ignored
all_cand_genomes <- unique(unlist(lapply(cand_haplotype, function(cand_entry){
  lapply(cand_entry$seq_haplotype_list, function(haplotype_entry) haplotype_entry$haplotype_genome)
})))

message(
  sum(strain_compo$genome %in% all_cand_genomes), " of ", nrow(strain_compo), " strain_compo genomes appear in cand_haplotype; ",
  sum(!all_cand_genomes %in% strain_compo$genome), " cand_haplotype genomes are ignored"
)

## Haplotype stats and sensitivity and specificity by cluster ----
# For each haplotype and each cluster carrying it, counting only genomes typed at the candidate: 
# tp = cluster genomes with the haplotype
# fn = cluster genomes without it
# fp = other clusters' genomes with it
# tn = other clusters' genomes without it.
# Rows cover NC_ strains too (is_cluster = FALSE); the *_clusters percentages count all groups

# Sensitivity equals the cluster's haplotype frequency used below
if(file.exists("haplotype_cluster_stats.RDS")){
  haplotype_cluster_stats <- readRDS("haplotype_cluster_stats.RDS")
}else{
  
  haplotype_cluster_stats <- do.call(rbind,lapply(cand_names, function(cand_name_entry){
    
    cand_entry <- cand_haplotype[[cand_name_entry]]
    cand_name <- cand_name_entry
    # strain_compo genomes typed at this candidate (carrying any of its haplotypes)
    cand_genomes <- unlist(lapply(cand_entry$seq_haplotype_list, function(haplotype_entry) haplotype_entry$haplotype_genome))
    cand_strain_compo <- strain_compo[strain_compo$genome %in% cand_genomes, ]
    n_cand_genome <- nrow(cand_strain_compo)
    n_cand_clusters <- length(unique(cand_strain_compo$cluster))
    # Cluster specificity and sensitivity of each haplotype if applicable 
    cand_cluster_stats <- do.call(rbind, mclapply(cand_entry$seq_haplotype_list, function(haplotype_entry){
      
      haplotype_genomes <- haplotype_entry$haplotype_genome
      haplotype_genomes <- haplotype_genomes[haplotype_genomes %in% strain_compo$genome]
      
      haplotype_clusters <- unique(strain_compo$cluster[
        strain_compo$genome %in% haplotype_genomes
      ])
      
      # no genome carrying this haplotype is assigned to a strain; rbind skips NULL
      if(length(haplotype_clusters) == 0) return(NULL)
      
      # one value per typed genome: TRUE if it carries this haplotype
      genome_has_haplotype <- cand_strain_compo$genome %in% haplotype_genomes
      
      
      pct_cand_genomes <- 100 * length(haplotype_genomes) / n_cand_genome
      pct_total_genomes <- 100 * length(haplotype_genomes) / nrow(strain_compo)
      pct_cand_clusters <- 100 * length(haplotype_clusters) / n_cand_clusters
      pct_total_clusters <- 100 * length(haplotype_clusters) / length(total_cluster)
      
      haplotype_clusters_sensit_specif <- do.call(rbind, lapply(haplotype_clusters, function(cluster_entry){
        
        # one value per typed genome: TRUE if it belongs to this cluster
        genome_in_cluster <- cand_strain_compo$cluster == cluster_entry
        
        tp <- sum(genome_in_cluster & genome_has_haplotype)
        fn <- sum(genome_in_cluster & !genome_has_haplotype)
        fp <- sum(!genome_in_cluster & genome_has_haplotype)
        tn <- sum(!genome_in_cluster & !genome_has_haplotype)
        
        data.frame(
          cluster = cluster_entry,
          is_cluster = cluster_entry %in% thresher_cluster,
          tp = tp,
          fn = fn,
          fp = fp,
          tn = tn,
          sensitivity = tp / (tp + fn),
          specificity = tn / (tn + fp)
        )
      }))
      
      cbind(
        candidate = cand_name,
        haplotype = haplotype_entry$haplotype_idx,
        haplotype_clusters_sensit_specif,
        pct_cand_genomes = pct_cand_genomes,
        pct_total_genomes = pct_total_genomes,
        pct_cand_clusters = pct_cand_clusters,
        pct_total_clusters = pct_total_clusters)
      
    }, mc.cores = detectCores() - 2))
  }))
  
  saveRDS(haplotype_cluster_stats,"haplotype_cluster_stats.RDS")
}

## Cluster x haplotype frequency by candidate ----
if(file.exists("cluster_haplotype_freq_list.RDS")){
  cluster_haplotype_freq_list <- readRDS("cluster_haplotype_freq_list.RDS")
}else{
  
  cluster_haplotype_freq_list <- mclapply(cand_haplotype, function(cand_entry){
    
    # data frame showing clusters that have this candidate by genome count
    cand_genomes <- unique(unlist(lapply(cand_entry$seq_haplotype_list, function(haplotype_entry) haplotype_entry$haplotype_genome)))
    cand_cluster_freq_dt <- as.data.table(table(factor(genome_cluster[cand_genomes], levels = total_cluster)))
    colnames(cand_cluster_freq_dt) <- c("cluster","cand_genome_count")
    
    # genomes of each cluster carrying each haplotype
    
    cluster_haplotype_count <- do.call(rbind,lapply(cand_entry$seq_haplotype_list, function(haplotype_entry){
      
      haplotype_cluster_freq_dt <- as.data.frame(table(factor(genome_cluster[haplotype_entry$haplotype_genome], levels = total_cluster)))
      haplotype_cluster_freq_dt <- cbind(
        haplotype = haplotype_entry$haplotype_idx,
        haplotype_cluster_freq_dt
      )
      
      # haplotype_genome_count means 
      # the genomes that have this haplotype within this cluster / strain
      colnames(haplotype_cluster_freq_dt)[2:3] <- c("cluster","haplotype_genome_count")
      
      return(haplotype_cluster_freq_dt)
    }))
    
    # cand_genome_count means the total genomes within this cluster that have this candidate
    cluster_haplotype_count$cand_genome_count <- cand_cluster_freq_dt$cand_genome_count[
      match(cluster_haplotype_count$cluster, cand_cluster_freq_dt$cluster)
    ]
    
    # the percentage of the specific haplotype in the candidate 
    cluster_haplotype_count$cand_haplotype_pct <- sapply(seq_len(nrow(cluster_haplotype_count)),function(row_idx){
      if(cluster_haplotype_count$cand_genome_count[row_idx] == 0){
        0
      }else{
        cluster_haplotype_count$haplotype_genome_count[row_idx] /  cluster_haplotype_count$cand_genome_count[row_idx]
      }
    })
    
    # empty cluster x haplotype matrix: rows in total_cluster order, columns named by haplotype_idx
    haplotype_idx <- sapply(cand_entry$seq_haplotype_list, function(haplotype_entry) haplotype_entry$haplotype_idx)
    cluster_haplotype_freq <- matrix(
      0,
      nrow = length(total_cluster),
      ncol = length(haplotype_idx),
      dimnames = list(total_cluster, haplotype_idx)
    )
    
    # place each frequency at its [cluster, haplotype] cell, matched by name
    cluster_haplotype_freq[cbind(as.character(cluster_haplotype_count$cluster), as.character(cluster_haplotype_count$haplotype))] <-
      cluster_haplotype_count$cand_haplotype_pct
    
    cluster_haplotype_freq
    
  },mc.cores = detectCores()-2)
  
  saveRDS(cluster_haplotype_freq_list,"cluster_haplotype_freq_list.RDS")
}


## Cluster x haplotype incidence by candidate ----
# TRUE when the cluster carries the haplotype and we call it valid
cluster_haplotype_valid_list <- lapply(cluster_haplotype_freq_list, function(cluster_haplotype_freq){
  cluster_haplotype_freq > 0
})

## Cluster-pair separation by candidate ----
# only pairs involving at least one transmission cluster need separating;
# two NC_ strains are never compared
cluster_pair_idx <- which(upper.tri(diag(length(total_cluster))), arr.ind = TRUE)
cluster_pair_idx <- cluster_pair_idx[total_is_cluster[cluster_pair_idx[, "row"]] | total_is_cluster[cluster_pair_idx[, "col"]], , drop = FALSE]

if(file.exists("cluster_pair_df.RDS")){
  cluster_pair_df <- readRDS("cluster_pair_df.RDS")
  cluster_pair_sep <- readRDS("cluster_pair_sep.RDS")
}else{
  
  cluster_pair_df <- data.frame(
    cluster_subject = total_cluster[cluster_pair_idx[, "row"]],
    cluster_query = total_cluster[cluster_pair_idx[, "col"]],
    pair_type = ifelse(
      total_is_cluster[cluster_pair_idx[, "row"]] & total_is_cluster[cluster_pair_idx[, "col"]],
      "cluster-cluster",
      "cluster-NC"
    )
  )
  
  # Cluster pair x candidate logical matrix: TRUE if the candidate separates the pair
  cluster_pair_sep <- sapply(cluster_haplotype_valid_list, function(cluster_haplotype_mat){
    
    # haplotype rows of the two clusters in each pair, in cluster_pair_df row order
    pair_haplotype_subject <- cluster_haplotype_mat[cluster_pair_df$cluster_subject, , drop = FALSE]
    pair_haplotype_query <- cluster_haplotype_mat[cluster_pair_df$cluster_query, , drop = FALSE]
    
    # does each cluster carry at least one haplotype of this candidate?
    pair_typed_subject <- rowSums(pair_haplotype_subject) > 0
    pair_typed_query <- rowSums(pair_haplotype_query) > 0
    
    # number of haplotypes both clusters carry: TRUE where both rows are TRUE, counted per pair
    pair_n_shared <- rowSums(pair_haplotype_subject & pair_haplotype_query)
    
    # separated: both clusters typed and no haplotype in common;
    # unname drops the cluster names that rowSums carried over from the row lookups
    unname(pair_typed_subject & pair_typed_query & pair_n_shared == 0)
  })
  
  cluster_pair_df <- cbind(cluster_pair_df, cluster_pair_sep)
  cluster_pair_df$n_cand_sep <- rowSums(cluster_pair_sep)
  saveRDS(cluster_pair_sep,"cluster_pair_sep.RDS")
  saveRDS(cluster_pair_df,"cluster_pair_df.RDS")
  
}

## Unresolvable pairs ----
# no single candidate separates these pairs, so no candidate set can
# they are reported here and left out of the optimization
unresolved_pair_df <- cluster_pair_df[cluster_pair_df$n_cand_sep == 0, ]
message(
  nrow(unresolved_pair_df), " of ", nrow(cluster_pair_df), " pairs cannot be separated by any candidate (",
  sum(unresolved_pair_df$pair_type == "cluster-cluster"), " cluster-cluster, ",
  sum(unresolved_pair_df$pair_type == "cluster-NC"), " cluster-NC)"
)


## Greedy candidate set ----
cover_mat <- unique(cluster_pair_sep[cluster_pair_df$n_cand_sep > 0, , drop = FALSE])

## Best n-candidate combination ----
# The best combination of n candidates need not contain the individually best candidates,
# so combinations are compared exactly (branch and bound): a branch is dropped only when even
# its most optimistic completion separates fewer pairs than the best combination found so far.
# A combination only joins candidates typed in exactly the same genomes, so no genome carries
# a haplotype at some of its candidates but not at others (e.g. at cand_1 but not at cand_13)

# unique() collapses pairs that share a separation pattern, so one row of cover_mat can stand
# for many pairs; pair_weight counts them so coverage is measured in pairs, not patterns
resolvable_sep <- cluster_pair_sep[cluster_pair_df$n_cand_sep > 0, , drop = FALSE]
pattern_key <- apply(resolvable_sep, 1, function(row_entry) paste(which(row_entry), collapse = ","))
pattern_idx <- match(pattern_key, unique(pattern_key))
cover_mat <- resolvable_sep[!duplicated(pattern_idx), , drop = FALSE]
pair_weight <- tabulate(pattern_idx)

# genome x candidate: TRUE when the genome carries a haplotype of the candidate
genome_typed_mat <- sapply(colnames(cover_mat), function(cand_name_entry){
  cand_genomes <- unlist(lapply(cand_haplotype[[cand_name_entry]]$seq_haplotype_list, function(haplotype_entry) haplotype_entry$haplotype_genome))
  strain_compo$genome %in% cand_genomes
})
rownames(genome_typed_mat) <- strain_compo$genome

# candidates that separate exactly the same pairs and are typed in the same genomes are
# interchangeable: search one per group and report the group as "cand_a|cand_b"
cand_sep_key <- apply(cover_mat, 2, function(col_entry) paste(which(col_entry), collapse = ","))
cand_typed_key <- apply(genome_typed_mat, 2, function(col_entry) paste(which(col_entry), collapse = ","))
cand_key <- paste(cand_sep_key, cand_typed_key, sep = ";")
cand_group_idx <- match(cand_key, unique(cand_key))
cand_members <- split(colnames(cover_mat), cand_group_idx)
cand_rep <- !duplicated(cand_group_idx)
cover_cand <- cover_mat[, cand_rep, drop = FALSE]
colnames(cover_cand) <- sapply(cand_members, paste, collapse = "|")

# typing class: candidates typed in exactly the same genomes; a combination stays within one class
cand_class <- match(cand_typed_key[cand_rep], unique(cand_typed_key[cand_rep]))

# leave out candidates that separate no resolvable pair
cand_resolves_any <- colSums(cover_cand) > 0
cover_cand <- cover_cand[, cand_resolves_any, drop = FALSE]
cand_members <- cand_members[cand_resolves_any]
cand_class <- cand_class[cand_resolves_any]
message(
  sum(pair_weight), " resolvable pairs in ", nrow(cover_mat), " separation patterns; ",
  ncol(cover_mat), " candidates collapse to ", ncol(cover_cand), " distinct candidates in ",
  length(unique(cand_class)), " typing classes (largest: ", max(table(cand_class)), " candidates)"
)

# cover: pattern x candidate logical matrix; weight: pairs behind each pattern row
# returns the pairs separated by the best n_cand combination and all tied best combinations (up to max_ties)
find_best_combination <- function(cover, weight, n_cand, max_ties = 100, n_cores = max(1, detectCores() - 2)){
  
  best_n_pair <- -Inf
  best_set_list <- list()
  total_n_pair <- sum(weight)
  
  # a combination separating every pair is listed only if each member is needed,
  # i.e. each member is the only one separating at least one pair
  is_minimal_full_set <- function(set_cand){
    set_cover <- cover[, set_cand, drop = FALSE]
    only_member_cover <- set_cover & rowSums(set_cover) == 1
    all(colSums(only_member_cover) > 0)
  }
  
  record_set <- function(set_cand, set_n_pair){
    keep_set <- set_n_pair < total_n_pair || length(set_cand) == 1 || is_minimal_full_set(set_cand)
    if(set_n_pair > best_n_pair){
      best_n_pair <<- set_n_pair
      best_set_list <<- if(keep_set) list(set_cand) else list()
    }else if(set_n_pair == best_n_pair && keep_set && length(best_set_list) < max_ties){
      best_set_list[[length(best_set_list) + 1]] <<- set_cand
    }
  }
  
  # pairs each allowed candidate would newly separate, largest first; candidates adding nothing are dropped
  rank_gain <- function(allowed, uncovered){
    gain <- unname(colSums(cover[uncovered, allowed, drop = FALSE] * weight[uncovered]))
    gain_order <- order(gain, decreasing = TRUE)
    gain_order <- gain_order[gain[gain_order] > 0]
    list(allowed = allowed[gain_order], gain = gain[gain_order])
  }
  
  # the last two slots are solved exactly: a pair adds gain_a + gain_b minus the pairs both separate
  add_best_two <- function(chosen, ranked, uncovered, n_pair){
    
    n_pair_left <- sum(weight[uncovered])
    
    # a candidate that separates every remaining pair completes the combination alone
    for(cand_rank in which(ranked$gain == n_pair_left)){
      record_set(c(chosen, ranked$allowed[cand_rank]), n_pair + n_pair_left)
    }
    
    pair_cand <- ranked$allowed[ranked$gain < n_pair_left]
    pair_gain <- ranked$gain[ranked$gain < n_pair_left]
    if(length(pair_cand) < 2) return(invisible(NULL))
    
    # gains are sorted, so each candidate's best possible partner is the first one (the second, for the first itself);
    # drop candidates that fall short of the best even with that partner
    best_partner_gain <- c(pair_gain[2], rep(pair_gain[1], length(pair_gain) - 1))
    pair_keep <- n_pair + pair_gain + best_partner_gain >= best_n_pair
    pair_cand <- pair_cand[pair_keep]
    pair_gain <- pair_gain[pair_keep]
    if(length(pair_cand) < 2) return(invisible(NULL))
    
    # rows no remaining candidate separates add nothing to any pair
    pair_cover <- cover[uncovered, pair_cand, drop = FALSE]
    row_keep <- rowSums(pair_cover) > 0
    
    # both_n_pair[a, b]: remaining pairs separated by both a and b.
    # crossprod(x) uses the faster symmetric product; sqrt(weight) weights it and round() removes sqrt's rounding error
    weighted_cover <- pair_cover[row_keep, , drop = FALSE] * sqrt(weight[uncovered][row_keep])
    both_n_pair <- round(crossprod(weighted_cover))
    
    pair_n_pair <- n_pair + outer(pair_gain, pair_gain, "+") - both_n_pair
    pair_n_pair[lower.tri(pair_n_pair, diag = TRUE)] <- -Inf
    
    top_n_pair <- max(pair_n_pair)
    if(top_n_pair < best_n_pair) return(invisible(NULL))
    
    top_pair_idx <- which(pair_n_pair == top_n_pair, arr.ind = TRUE)
    for(top_idx in seq_len(nrow(top_pair_idx))){
      record_set(c(chosen, pair_cand[top_pair_idx[top_idx, ]]), top_n_pair)
    }
  }
  
  # chosen: candidates already in the combination; allowed: candidates this branch may still add;
  # uncovered: pattern rows not yet separated; n_pair: pairs separated by chosen
  search_branch <- function(chosen, allowed, uncovered, n_pair){
    
    slots_left <- n_cand - length(chosen)
    
    # full combination, or every resolvable pair already separated (more candidates would add nothing)
    if(slots_left == 0 || length(uncovered) == 0) return(record_set(chosen, n_pair))
    
    ranked <- rank_gain(allowed, uncovered)
    
    if(slots_left == 2){
      add_best_two(chosen, ranked, uncovered, n_pair)
    }else{
      search_children(chosen, ranked, uncovered, n_pair, seq_along(ranked$allowed))
    }
  }
  
  # add each candidate at child_rank in turn and search deeper; a branch only adds lower-ranked
  # candidates afterwards, so every combination is visited once
  search_children <- function(chosen, ranked, uncovered, n_pair, child_rank){
    
    slots_left <- n_cand - length(chosen)
    
    for(cand_rank in child_rank){
      
      # a candidate's gain can only shrink as others join, so this candidate plus the next
      # (slots_left - 1) gains is the most this branch can reach; gains are sorted, so the bound
      # only drops for later candidates
      bound_rank <- cand_rank:min(length(ranked$gain), cand_rank + slots_left - 1)
      if(n_pair + sum(ranked$gain[bound_rank]) < best_n_pair) break
      
      cand_idx <- ranked$allowed[cand_rank]
      search_branch(
        chosen = c(chosen, cand_idx),
        allowed = ranked$allowed[-seq_len(cand_rank)],
        uncovered = uncovered[!cover[uncovered, cand_idx]],
        n_pair = n_pair + ranked$gain[cand_rank]
      )
    }
  }
  
  all_cand <- seq_len(ncol(cover))
  all_row <- seq_len(nrow(cover))
  
  if(n_cand <= 2 || n_cores <= 1){
    search_branch(chosen = integer(0), allowed = all_cand, uncovered = all_row, n_pair = 0)
  }else{
    # branches starting from different first candidates are independent: deal them out across cores
    # (interleaved so strong and weak first candidates are spread evenly) and merge the results
    root_ranked <- rank_gain(all_cand, all_row)
    root_rank_chunk <- split(seq_along(root_ranked$allowed), rep_len(seq_len(n_cores), length(root_ranked$allowed)))
    
    chunk_result_list <- mclapply(root_rank_chunk, function(rank_chunk){
      search_children(chosen = integer(0), ranked = root_ranked, uncovered = all_row, n_pair = 0, child_rank = rank_chunk)
      list(n_pair = best_n_pair, set_list = best_set_list)
    }, mc.cores = n_cores)
    
    chunk_failed <- sapply(chunk_result_list, inherits, "try-error")
    if(any(chunk_failed)) stop("search failed on a core: ", chunk_result_list[[which(chunk_failed)[1]]])
    
    best_n_pair <- max(sapply(chunk_result_list, function(chunk_entry) chunk_entry$n_pair))
    best_set_list <- do.call(c, lapply(chunk_result_list, function(chunk_entry){
      if(chunk_entry$n_pair == best_n_pair) chunk_entry$set_list else list()
    }))
    best_set_list <- head(best_set_list, max_ties)
  }
  
  list(n_cand = n_cand, n_pair = best_n_pair, set_list = best_set_list)
}

# best combination of n_cand candidates from a single typing class: each class is searched on the
# pairs its own candidates can separate, and the best combinations across classes are kept
find_best_typed_combination <- function(cover, weight, n_cand, cand_class, max_ties = 100, n_cores = max(1, detectCores() - 2)){
  
  class_result_list <- lapply(unique(cand_class), function(class_entry){
    
    class_cand <- which(cand_class == class_entry)
    # a class with fewer candidates than n_cand cannot form the combination
    if(length(class_cand) < n_cand) return(NULL)
    
    class_row <- which(rowSums(cover[, class_cand, drop = FALSE]) > 0)
    class_result <- find_best_combination(
      cover[class_row, class_cand, drop = FALSE], weight[class_row], n_cand,
      max_ties = max_ties,
      # forking only pays off for large classes
      n_cores = if(length(class_cand) >= 50) n_cores else 1
    )
    
    # candidate indices back to columns of cover
    list(
      n_pair = class_result$n_pair,
      set_list = lapply(class_result$set_list, function(set_entry) class_cand[set_entry])
    )
  })
  class_result_list <- class_result_list[!sapply(class_result_list, is.null)]
  
  if(length(class_result_list) == 0){
    message("No typing class has ", n_cand, " candidates, so no combination of ", n_cand, " is possible")
    return(list(n_cand = n_cand, n_pair = NA, set_list = list()))
  }
  
  best_n_pair <- max(sapply(class_result_list, function(class_entry) class_entry$n_pair))
  best_set_list <- do.call(c, lapply(class_result_list, function(class_entry){
    if(class_entry$n_pair == best_n_pair) class_entry$set_list else list()
  }))
  
  list(n_cand = n_cand, n_pair = best_n_pair, set_list = head(best_set_list, max_ties))
}

n_cand_range <- 1:3
best_combo_list <- lapply(n_cand_range, function(n_cand_entry){
  find_best_typed_combination(cover_cand, pair_weight, n_cand_entry, cand_class)
})

## Best combination summary ----
# one row per tied best combination; fewer than n candidates are listed when they already separate
# every pair their typing class can separate (any other candidate of the class fills the rest)
best_combo_df <- do.call(rbind, lapply(best_combo_list, function(combo_entry){
  do.call(rbind, lapply(seq_along(combo_entry$set_list), function(tie_idx){
    
    set_entry <- combo_entry$set_list[[tie_idx]]
    # any member of an interchangeable group separates the same pairs
    set_cand <- sapply(cand_members[set_entry], function(member_entry) member_entry[1])
    pair_resolved <- rowSums(cluster_pair_sep[, set_cand, drop = FALSE]) > 0
    
    data.frame(
      n_cand = combo_entry$n_cand,
      tie_idx = tie_idx,
      candidates = paste(colnames(cover_cand)[set_entry], collapse = " + "),
      # genomes with a haplotype at every candidate of the combination
      n_genome_typed = sum(rowSums(genome_typed_mat[, set_cand, drop = FALSE]) == length(set_cand)),
      n_pair_resolved = sum(pair_resolved),
      n_cluster_cluster_resolved = sum(pair_resolved & cluster_pair_df$pair_type == "cluster-cluster"),
      n_cluster_nc_resolved = sum(pair_resolved & cluster_pair_df$pair_type == "cluster-NC"),
      pct_resolvable_pairs = 100 * sum(pair_resolved) / sum(pair_weight),
      pct_all_pairs = 100 * sum(pair_resolved) / nrow(cluster_pair_df)
    )
  }))
}))
print(best_combo_df)

# Stage 04 reads this to know which loci to report on
saveRDS(best_combo_df, "best_combo_df.RDS")

# pairs still unresolved by the first best combination of the largest n
top_combo <- best_combo_list[[length(best_combo_list)]]
if(length(top_combo$set_list) > 0){
  top_set_cand <- sapply(cand_members[top_combo$set_list[[1]]], function(member_entry) member_entry[1])
  top_set_unresolved_df <- cluster_pair_df[
    rowSums(cluster_pair_sep[, top_set_cand, drop = FALSE]) == 0,
    c("cluster_subject", "cluster_query", "pair_type", "n_cand_sep")
  ]
  print(table(top_set_unresolved_df$pair_type))
}