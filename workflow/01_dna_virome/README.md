# Module 01: DNA virome

## Purpose

Generate non-redundant DNA viral operational taxonomic units (vOTUs), their
relative abundance profiles, taxonomic assignments, predicted lifestyles, and
auxiliary metabolic gene annotations.

## Inputs

- paired DNA metavirome FASTQ files;
- DNA sample metadata; and
- required tool and reference-database versions.

## analysis

1. Quality-control paired reads with fastp (`Q20`, minimum length 20,
   poly-G/poly-X trimming).
2. Co-assemble reads with metaSPAdes.
3. Retain contigs at least 3,000 bp long with SeqKit.
4. Identify candidate viral contigs using VirSorter2 and DeepVirFinder.
5. Apply the custom candidate-reconciliation filter and assess quality using
   CheckV.
6. Pool retained viral contigs and perform all-versus-all BLASTn.
7. Calculate pairwise ANI/alignment fractions using the CheckV-distributed
   `anicalc.py` utility.
8. Cluster DNA vOTUs using the CheckV-distributed `aniclust.py` utility with
   `--min_ani 95 --min_tcov 85 --min_qcov 0`.
9. Extract the first-column representative IDs from the cluster table using
   SeqKit.
10. Calculate per-sample counts with CoverM and convert them to TPM.
11. Assign taxonomy and predict lifestyle/AMG potential using the finalized
    taxonomy workflow, DeePhage, and DRAM-v.

## Outputs

- pooled screened DNA viral contigs;
- DNA vOTU cluster table and representative FASTA;
- CheckV quality table;
- count and TPM abundance matrices;
- taxonomy and predicted lifestyle tables; and
- DRAM-v annotation and AMG tables.

