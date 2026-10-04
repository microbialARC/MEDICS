# Reference entropy profile

SILOSim `profiler` output for the reference genome used throughout this analysis:

**_Staphylococcus aureus_ subsp. _aureus_ NCTC 8325** — assembly
[GCF_000013425.1][assembly] (ASM1342v1).

Produced by `analysis/00_profile_reference.sh`. Filenames keep the prefix the
profiler writes (`<assembly>_<output>`), so they can be matched against a fresh
profiler run without renaming.

| File | Contents | Read by |
| --- | --- | --- |
| `*_entropy.RDS` | Shannon entropy at every reference position | stage 01 |
| `*_new_pos.RDS` | Contig boundaries in concatenated coordinates | stage 01 |
| `*_mges.RDS` | Mobile genetic element intervals, excluded from the scan | stage 01 |
| `*_concat.fasta` | Concatenated reference sequence | stages 01, 02.1 |
| `*_chr_bins.RDS` | Gene annotation, for the candidate figures | stage 01 |
| `*_mges.csv` | MGE intervals, same content in readable form | — |
| `*_mges_seq.fasta` | MGE sequences | — |
| `*_bin_summary.csv` | Per-bin annotation summary | — |
| `mge_entropy/` | Per-MGE entropy, one file per element (149) | — |

The entropy profile reflects the diversity of the **collection** aligned against
this reference, not of the reference itself. It is therefore specific to the CHOP
NICU *S. aureus* collection: re-profiling a different collection against the same
reference gives different entropy and, in general, different candidate loci.

Large profiler outputs that no later stage reads — per-genome SNP tables,
position coverage, plot frames — are not committed. See
[../../docs/data-availability.md](../../docs/data-availability.md).

[assembly]: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013425.1/
