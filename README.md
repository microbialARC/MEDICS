# MEDICS

**M**aximum **E**ntropy **D**iagnostic **I**nfection **C**ontrol **S**ystem

## Aim

Whole-genome sequencing (WGS) is the gold standard for resolving
healthcare-associated transmission, but its long turnaround time and high cost
limit its scalability for routine surveillance. The objective of MEDICS is to
determine whether a **single maximum-entropy locus** can serve as a first-line
alternative to WGS for surveillance. Our central hypothesis is that this
locus is inexpensive enough to sequence from every isolate while retaining
sufficient discriminatory power to prioritize isolates for confirmatory WGS.

**Approach.** Our pipeline, [SILOSim](https://github.com/microbialARC/SILOSim),
identifies the target locus by scanning per-site Shannon entropy along the
reference genome for the most variable window flanked by zero-entropy regions.
The **conserved flanks** provide primer-binding sites for amplification, and the
variable interior carries the typing signal. Isolates are assigned *haplotypes*
based on exact sequence identity at the locus, eliminating the need for distance
thresholds or clustering.

## Preliminary findings

Benchmarked against WGS on the CHOP NICU *S. aureus* collection, including 1,594 isolates
over 2.5 years (Sep 2022 – Mar 2025), with 101 transmission clusters identified by
[THRESHER](https://github.com/microbialARC/THRESHER).


| Finding | Value |
| --- | --- |
| Isolates typed by the single locus | 1,060 of 1,594 |
| Distinct haplotype types | 519 |
| Transmission clusters reached | 93 of 101 (**92.08%**) |
| Sequence types (MLST) among typed isolates | 58 |
| STs subdivided by ST-specific haplotypes | 46 |
| ST-specific haplotypes per ST | median 2.5 (IQR 1–8.75), mean 9.83 (SD 20.47), max 124 (ST398) |

See details at [docs/results.md](docs/results.md).

## Repository map

```
analysis/                     the pipeline, in the order it runs
  config.R                    all paths and thresholds; every script sources this
  00_profile_reference.sh     SILOSim profiler to per-site entropy  (stage 00)
  01_scan_loci/               entropy scan to candidate loci        (stage 01)
  02_type_isolates/           BLAST, extract, align, call haplotypes (stage 02)
  03_select_loci/             which loci separate which clusters     (stage 03)
  04_benchmark/               benchmark vs WGS and MLST, figures     (stage 04)
data/reference/               SILOSim profile of S. aureus NCTC 8325 (GCF_000013425.1)
docs/                         results, data availability
results/figures/              candidate locus plots, Sankey, monthly capacity
results/tables/               candidate table and aggregate benchmark tables
```

Only **stage 01 runs on what is in this repository**. Stages 00 and 02–04 need
the CHOP NICU isolate collection, which is not shared here. See
[docs/data-availability.md](docs/data-availability.md).

## Reproducing the locus scan

Needs R with `dplyr`, `Biostrings`, `IRanges`, `ggplot2`, `patchwork` and
`openxlsx2`. From the repository root:

```bash
Rscript analysis/01_scan_loci/run_locus_scan.R
```

This re-derives the 78 candidate loci in `results/tables/locus_candidates.xlsx`
from the entropy profile in `data/reference/`, and redraws each candidate's
figure. 

The benchmark locus is row 13 of that table: a 500-bp region at
contig_1:549,440–549,939 with 31 informative sites and 100 bp of invariant
flank on each side. It lies within the coding sequence of the Ser-Asp rich
fibrinogen-binding bone sialoprotein-binding protein.
