# Test plan

The repository does not yet contain publishable test data. A minimal dataset
should be added only after confirming that it contains no sensitive sample or
partner information.

The first automated checks should cover:

1. paired-end sample-sheet validation;
2. feature/count/length identifier consistency;
3. TPM columns summing to approximately one million for non-empty samples;
4. sample alignment across abundance and metadata tables;
5. deterministic vOTU representative selection;
6. exact CRISPRone spacer-to-host-contig ID recovery and validation against the
   CD-HIT-EST representative FASTA;
7. the legacy CRISPR hit filter (identity, mismatch, and bitscore) plus any
   separately labelled coverage/gap sensitivity filter;
8. expected row counts at major workflow checkpoints; and
9. R script parsing in the declared Conda environment;
10. `estimateR()` output dimensions and Pielou calculation on a hand-checkable
    sample-by-feature count matrix;
11. source tests restricted within `pair_id` and season tests operating on
    independent dates;
12. BH correction applied within each declared feature/correlation family;
13. unique upper-triangle temporal pairs built from `sampling_date`; and
14. AMG genomic counts de-duplicated across repeated sample-abundance rows.

Once the Conda environment is available, run the current syntax gate from the
repository root:

```bash
Rscript tests/parse_r_scripts.R
```
