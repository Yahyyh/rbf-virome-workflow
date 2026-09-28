# Curated R analysis modules

The public scripts replace the original interactive R notebooks. They do not use `setwd()`, personal paths, spreadsheet-only intermediates, manual P-value
labels, or hard-coded sample ordering. Every module writes tab-separated source tables, SVG/PDF figures where applicable, and `sessionInfo.txt`.

## Input conventions

- Abundance matrices are **feature by sample**: the first column contains a
  stable feature ID and every remaining column is a sample.
- Metadata joins use exact `sample_id` values. Group membership is never parsed
  from filename substrings.
- DNA and RNA matrices are analysed separately. The synthetic
  `metadata/samples.example.tsv` demonstrates the schema only; its dates and
  paths are not study data.
- `metadata/samples.ncbi-submission.tsv` maps the 108 submitted sample aliases
  to their library aliases and paired FASTQ filenames. Before using it as
  `metadata/samples.tsv`, filter to one assay, populate the exact
  `sampling_date`, and verify `pair_id_candidate`; the NCBI submission
  workbooks record collection month only and do not explicitly encode pairing.
- `pair_id` identifies river and groundwater samples collected on the same
  date. Source contrasts within a season use this pairing; season contrasts
  compare different sampling dates and are unpaired.
- For the AMG long-table module, pass metadata already filtered to the intended
  nucleic-acid dataset and include zero-TPM rows so every intended sample is
  represented. This prevents a completely negative sample from disappearing.
- Count-based richness estimators require non-negative integer counts. TPM or
  other abundance is converted to sample-wise relative abundance for
  composition analyses. TPM is not an absolute concentration.

## Modules

| Script | Purpose | Main safeguards |
|---|---|---|
| `01_abundance_diversity.R` | Count-based richness, Shannon/Simpson/Pielou, Bray-Curtis PCoA, planned PERMANOVA, PERMDISP | count-only richness estimators; natural-log Pielou; date-pair permutation restrictions; Lingoes-corrected PCoA; BH correction |
| `02_virus_host_procrustes.R` | Community-level host/virus Procrustes and PROTEST | sample-wise Hellinger transform; symmetric fit; fixed seed; within-group sensitivity analysis; no claim of individual-host validation |
| `03_amg_lifestyle.R` | AMG load, abundance, composition, and predicted-lifestyle association | genome-level de-duplication; zero-completed sample summaries; virus-label permutation test; standardized residuals |
| `04_coverm_counts_to_tpm.R` | Merge per-sample counts and calculate TPM | fails on unknown/duplicate IDs; writes length separately; validates each TPM library sum |
| `05_differential_abundance.R` | Planned feature-wise source and season comparisons | paired source tests; prevalence filter; BH within contrast; no invalid confidence intervals |
| `06_virus_host_covariation.R` | Spearman covariation for **a priori** predicted links | does not search all possible pairs; prevalence/constant-vector checks; BH by stratum; covariation is not host validation |
| `07_lifestyle_composition.R` | Predicted-lifestyle abundance, taxonomy, and within-class diversity | explicit unclassified category; no duplicated genome counts; paired source tests; correct Pielou base |
| `08_group_intersections.R` | UpSet-style source-by-season presence intersections | metadata-driven groups; documented detection threshold; one membership row per feature |
| `09_temporal_turnover.R` | Bray-Curtis time-decay by source | real `sampling_date`; unique upper-triangle pairs; Mantel permutations instead of treating pairwise distances as independent replicates |

All scripts source `scripts/r/lib/common.R`, which fixes the random seed at
`20260925`, validates schemas, aligns samples, implements restricted
permutations, and saves vector figures.

## Commands

```bash
# CoverM counts to TPM
Rscript scripts/r/04_coverm_counts_to_tpm.R \
  input/lengths.tsv input/counts results/dna

# Diversity and community structure
Rscript scripts/r/01_abundance_diversity.R \
  results/dna.counts.tsv results/dna.tpm.tsv metadata/dna.tsv \
  results/dna_diversity dna

# Community-level virus-host concordance
Rscript scripts/r/02_virus_host_procrustes.R \
  input/virus_abundance.tsv input/host_abundance.tsv metadata/dna.tsv \
  results/procrustes dna_host

# AMG/lifestyle analysis
Rscript scripts/r/03_amg_lifestyle.R \
  input/amg_abundance_long.tsv metadata/dna.tsv results/amg dna

# Feature-wise planned comparisons
Rscript scripts/r/05_differential_abundance.R \
  results/dna.tpm.tsv metadata/dna.tsv results/differential dna 0.20

# Covariation only among precomputed putative links
Rscript scripts/r/06_virus_host_covariation.R \
  input/virus_abundance.tsv input/host_abundance.tsv input/predicted_links.tsv \
  metadata/dna.tsv results/covariation dna_host 0.20

# Predicted-lifestyle composition
Rscript scripts/r/07_lifestyle_composition.R \
  results/dna.tpm.tsv input/virus_annotations.tsv metadata/dna.tsv \
  results/lifestyle dna

# Source-by-season intersections
Rscript scripts/r/08_group_intersections.R \
  results/dna.tpm.tsv metadata/dna.tsv results/intersections dna 0

# Temporal turnover
Rscript scripts/r/09_temporal_turnover.R \
  results/dna.tpm.tsv metadata/dna.tsv results/temporal dna
```

