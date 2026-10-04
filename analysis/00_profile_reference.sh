#!/usr/bin/env bash
# Stage 00 - Build the per-site entropy profile of the reference genome.
#
# SILOSim's profiler maps every genome in the collection against the reference
# assembly, calls SNPs, annotates mobile genetic elements (MGEs), and writes the
# per-site Shannon entropy that stage 01 scans. Its outputs are the files in
# data/reference/, which are distributed with this repository, so stage 01 can be
# re-run without repeating this step.
#
# Needs the genome assemblies (not distributed) and the Bakta / WhatsGNU
# databases. See docs/data-availability.md.
#
# SILOSim: https://github.com/microbialARC/SILOSim
set -euo pipefail

conda activate silosim

silosim profiler \
  --input_genome GCF_000013425.1_ASM1342v1_genomic.fna \
  --output profiler_all_nicu/ \
  --local_query_dir all_nicu/ \
  --species sau \
  --bakta_db_path bakta_db/db/ \
  --whatsgnu_db_path Sau_WhatsGNU_Ortholog_db/Sau_WhatsGNU_Ortholog_db.pickle \
  -t 6 \
  --conda_prefix conda_envs/

# Of the profiler's outputs, stage 01 reads five:
#   *_entropy.RDS     per-site Shannon entropy
#   *_new_pos.RDS     contig boundaries in concatenated coordinates
#   *_mges.RDS        MGE intervals, excluded from the locus search
#   *_concat.fasta    the concatenated reference
#   *_chr_bins.RDS    gene annotation, for the locus figures
# These are in data/reference/. The remaining outputs (per-genome SNP tables,
# position coverage, plot frames) are large and regenerable, so they are left out
# of the repository; re-run this command to recreate them.
