# Module 02: RNA virome

## Purpose

Generate non-redundant RNA vOTUs and their abundance and taxonomic profiles.

## Inputs

- paired RNA metavirome FASTQ files;
- RNA sample metadata;
- SortMeRNA, VirBot, CheckV, Stampede-ClusterGenomes, CoverM, and EsViritu
  resources.

## Analysis

1. Quality-control reads with fastp.
2. Remove rRNA reads using SortMeRNA.
3. Co-assemble reads and retain contigs at least 1,000 bp long.
4. Identify RNA viral contigs using VirBot.
5. Run CheckV and extract the retained viral contigs.
6. Pool RNA viral contigs across samples.
7. Activate the `stampede-clustergenomes` environment and run:

   ```bash
   Cluster_genomes.pl -f all_rnavir.fasta -c 85 -i 95
   ```

8. Reclassify representative sequences using VirBot in sensitive mode.
9. Calculate CoverM counts and TPM relative abundance.
10. Identify eukaryotic RNA viral signals using EsViritu.

## Outputs

- screened RNA viral contigs;
- RNA vOTU clusters and representative sequences;
- count and TPM abundance matrices; and
- VirBot and EsViritu classification tables.

