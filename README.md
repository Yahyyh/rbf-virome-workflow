# Riverbank filtration virome workflow

This repository is contains the computational analyses supporting the manuscript **“Riverbank filtration couples viral attenuation with
ecological selection.”** The analyses are organized into separate DNA virome, RNA virome, metagenome, and virus–host association modules.

## Four workflow modules

| Module | Start here |
|---|---|
| 01 DNA virome | [`workflow/01_dna_virome/`](workflow/01_dna_virome/) |
| 02 RNA virome | [`workflow/02_rna_virome/`](workflow/02_rna_virome/) |
| 03 Prokaryotic metagenome | [`workflow/03_metagenome/`](workflow/03_metagenome/) |
| 04 Virus-host associations | [`workflow/04_virus_host/`](workflow/04_virus_host/) |

## Repository layout

```text
.
├── config/
│   ├── config.example.yaml
│   ├── dna-virome.example.yaml
│   ├── rna-virome.example.yaml
│   ├── metagenome.example.yaml
│   └── virus-host.example.yaml
├── workflow/
│   ├── 01_dna_virome/
│   ├── 02_rna_virome/
│   ├── 03_metagenome/
│   └── 04_virus_host/
├── scripts/r/
├── envs/
├── metadata/                 # Example sample sheet
├── docs/
│   ├── workflow.md
│   └── r-analysis.md
├── tests/                    # Validation plan
```

## Shared downstream R modules

Create the analysis environment:

```bash
conda env create -f envs/r-analysis.yml
conda activate rbf-r-analysis
```

The curated modules cover CoverM count-to-TPM conversion, alpha/beta diversity, restricted-permutation PERMANOVA and PERMDISP, virus-host
Procrustes/PROTEST, AMG/lifestyle analyses, planned differential abundance, putative-link covariation, lifestyle composition, group intersections, and
temporal turnover. See the [`R analysis guide`](docs/r-analysis.md) for exact input schemas, commands, outputs, and the legacy analyses deliberately excluded
on statistical grounds. TPM is relative library abundance, not absolute concentration.

## Reproducibility boundary

- `metadata/samples.example.tsv` is a synthetic schema example; its year-2000 dates and file paths do not represent actual study observations.
- `metadata/samples.ncbi-submission.tsv` contains 108 records linking sample aliases one-to-one with SRA library aliases, together with assay classes,
  collection months, and paired-end FASTQ filenames, as recorded in the BioSample and SRA submission workbooks.

## Data availability

The manuscript lists NCBI BioProject **PRJNA1493781**. 
