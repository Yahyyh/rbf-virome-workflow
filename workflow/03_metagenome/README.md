# Module 03: metagenome

## Purpose

Assemble the metagenome, predict genes, characterize community
composition and functional potential, and provide contigs for host inference.

## Inputs

- paired whole-metagenome FASTQ files;
- metagenome sample metadata; and
- UniVec and functional/taxonomic reference databases.

## Analysis

1. Quality-control reads with fastp. The metagenome section uses Q30 and a
   minimum read length of 35.
2. Remove UniVec-matching reads using BWA.
3. Assemble reads with metaSPAdes using k-mers
   21, 33, 55, 77, 95, 107, and 121.
4. Evaluate assembly statistics with QUAST.
5. Retain contigs at least 300 bp long and rename them using anvi'o.
6.Predict ORFs from the filtered metagenomic contigs using Prodigal in metagenomic 
   mode (-p meta), generating protein FASTA, nucleotide CDS FASTA, and GFF files.
7. Construct a non-redundant nucleotide gene catalogue using CD-HIT at 95%
   identity.
8. Quantify metagenome features using CoverM.
9. Perform functional annotation with eggNOG-mapper and contig taxonomy with
   CAT where required.

## Outputs

- quality-controlled metagenomic reads;
- assembled and filtered prokaryotic contigs;
- predicted proteins (`.faa`), gene nucleotide sequences (`.fna`), gene
  coordinates (`.gff`), and a non-redundant gene catalogue;
- abundance matrices; and
- functional and taxonomic annotation tables.
