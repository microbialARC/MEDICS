#!/usr/bin/env bash
# Stage 02.4 - Align each locus across the isolates carrying it.
#
# One MAFFT alignment per locus, over the per-locus FASTA from stage 02.3. The
# alignment is what makes haplotypes comparable: stage 02.5 calls a haplotype by
# collapsing identical aligned sequences, so indels have to be represented as
# gaps in a common frame rather than as different sequence lengths.
#
# One SLURM array task per locus.
#   sbatch --array=1-78 analysis/02_type_isolates/04_align_loci.sh
#
#SBATCH -c 10
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --job-name=medics_locus_mafft
set -euo pipefail

module load MAFFT

: "${INPUT_DIR:?Set INPUT_DIR to the per-locus FASTAs from stage 02.3}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR to the directory for alignments}"

mkdir -p "$OUTPUT_DIR"

LOCUS="cand_${SLURM_ARRAY_TASK_ID}"
INPUT_FASTA="${INPUT_DIR}/${LOCUS}_extract_seq.fasta"

if [[ ! -f "$INPUT_FASTA" ]]; then
  echo "No extracted sequences for ${LOCUS} at ${INPUT_FASTA}" >&2
  exit 1
fi

echo "[$(date -Is)] Aligning ${LOCUS}"
mafft --auto --thread 10 "$INPUT_FASTA" > "${OUTPUT_DIR}/${LOCUS}_aln.fasta"
echo "[$(date -Is)] Done."
