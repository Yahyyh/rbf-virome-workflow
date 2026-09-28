# Module 04: Putative virus-host associations

## Purpose

Infer putative hosts for DNA vOTUs from CRISPR spacer matches, classify the
matched host contigs, and characterize their abundance and functional
potential. CRISPR matches record potential historical virus-host interactions;
they do not demonstrate active infection or experimentally confirmed host
range.

## Inputs

- DNA vOTU representative sequences from Module 01;
- metagenomic contigs retained at a minimum length of 300 bp from Module 03;
- the contig-renaming tables linking each contig to its source sample; and
- the CRISPRone command-line distribution and its reference files.

## Analysis

1. Merge the renamed metagenomic contigs and cluster them with CD-HIT-EST at
   95% sequence identity and at least 90% coverage of the shorter sequence
   (`-c 0.95 -aS 0.90 -n 10`). Retain the `.clstr` membership file and the
   original contig-to-sample mapping.
2. Predict CRISPR arrays, spacers, and Cas loci from the clustered host contigs
   with the CRISPRone command-line driver, `crisprone-local.py`. The local copy
   disabled the final cleanup block so that the CRISPRone-generated
   `spacer/all-spacer.fa` file was retained for downstream matching; the
   equivalent public configuration is recorded as `cleanup: false`.
3. Use CRISPRone's `spacer/all-spacer.fa` directly. No additional CD-HIT
   clustering of spacer sequences was performed.
4. Build a nucleotide BLAST database from the DNA vOTU representatives using
   `makeblastdb -dbtype nucl`.
5. Align the CRISPRone spacers against the DNA vOTU database using
   `blastn -task blastn-short`. Report query and subject IDs, percentage
   identity, alignment length, mismatches, gap openings, coordinates, E-value,
   bitscore, `qcovs`, and `qcovhsp`.
6. Apply the original result-producing hit rule: percentage identity at least
   97%, no more than one mismatch, and bitscore at least 50. The original
   filter did not impose a query-coverage or gap-opening threshold; these
   fields are retained so stricter sensitivity analyses can be performed
   without changing the documented legacy result.
7. Recover the source host contig for every retained spacer hit. The original
   R snippet removed the terminal underscore-delimited component of each
   spacer ID; a public rerun must first construct an explicit CRISPRone
   spacer-to-contig mapping and validate every recovered ID against the
   CD-HIT-EST representative FASTA.
8. Extract the matched host contigs with SeqKit.
9. Classify the matched host contigs with `CAT_pack contigs`, attach official
   taxon names with `CAT_pack add_names`, and generate a classification summary
   with `CAT_pack summarise`. The recorded CAT database directory was
   `20240422_CAT_nr`.




## Downstream host characterization

These analyses describe the predicted host contigs but do not provide
additional evidence for an individual virus-host link.

1. Quantify matched host contigs with CoverM in count mode using 95% minimum
   read identity and 90% minimum aligned-read coverage.
2. Obtain contig lengths with SeqKit and convert count tables to TPM through
   [`../../scripts/r/04_coverm_counts_to_tpm.R`](../../scripts/r/04_coverm_counts_to_tpm.R).
3. Predict host-contig genes with Prodigal in metagenomic mode (`-p meta`) and
   annotate the proteins with eggNOG-mapper (`-m diamond`, seed-ortholog
   E-value at most 1e-5).
4. Compare matched virus and host community structures using Procrustes and
   PROTEST through
   [`../../scripts/r/02_virus_host_procrustes.R`](../../scripts/r/02_virus_host_procrustes.R).

## Outputs

- CD-HIT-EST host-contig representatives and cluster-membership table;
- CRISPRone CRISPR/Cas annotations and `spacer/all-spacer.fa`;
- spacer-to-vOTU BLAST table and filtered putative association table;
- validated spacer-to-host-contig mapping;
- CAT_pack host-taxonomy tables;
- host-contig count and TPM abundance matrices;
- host functional-annotation tables; and
- optional Procrustes/PROTEST summaries.
