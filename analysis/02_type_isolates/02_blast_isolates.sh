#!/usr/bin/env bash
# Stage 02.2 - Query every isolate assembly against the candidate-locus database.
#
# One SLURM array task per isolate. Each task BLASTs one assembly against the
# database built in stage 02.1 and writes a tabular hit file; stage 02.3 reads
# those and decides which loci are present in which isolate.
#
# The alignment parameters are tuned for short, near-identical nucleotide
# matches within a species: a discontiguous-megablast template with a small word
# size, so a locus is still found when it carries several SNPs.
#
# Submit with:
#   sbatch --array=1-$(wc -l < "$MANIFEST") analysis/02_type_isolates/02_blast_isolates.sh
#
#SBATCH -c 1
#SBATCH --mem=500M
#SBATCH --time=0:30:00
#SBATCH --job-name=medics_locus_blast
set -euo pipefail

module load BLAST+/2.14.0-gompi-2022b.lua

# Set these for your site. MANIFEST lists one assembly path per line.
: "${LOCUS_DB:?Set LOCUS_DB to the makeblastdb-indexed FASTA from stage 02.1}"
: "${MANIFEST:?Set MANIFEST to a file listing one assembly path per line}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR to the directory for BLAST hit tables}"

ASSEMBLY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
if [[ -z "${ASSEMBLY:-}" || ! -f "$ASSEMBLY" ]]; then
  echo "No assembly for task ${SLURM_ARRAY_TASK_ID} (got: '${ASSEMBLY}')" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
GENOME=$(basename "${ASSEMBLY%.*}")

echo "[$(date -Is)] Task ${SLURM_ARRAY_TASK_ID}: $GENOME"

blastn \
  -word_size 11 \
  -gapopen 5 -gapextend 2 \
  -reward 2 -penalty -3 \
  -template_type coding -template_length 18 \
  -window_size 40 \
  -num_threads 1 \
  -outfmt 6 \
  -db "$LOCUS_DB" \
  -query "$ASSEMBLY" \
  -out "${OUTPUT_DIR}/${GENOME}_blastn.tsv"

echo "[$(date -Is)] Done."
