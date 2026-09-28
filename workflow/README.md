# Workflow modules

The computational workflow is divided into four modules because they use
different inputs, software, thresholds, and outputs.

| Module | Main purpose | Documentation |
|---|---|---|
| 01 | DNA virome assembly, detection, vOTU clustering, abundance and annotation | [`01_dna_virome/`](01_dna_virome/) |
| 02 | RNA virome assembly, VirBot detection, Stampede clustering and annotation | [`02_rna_virome/`](02_rna_virome/) |
| 03 | Prokaryotic metagenome assembly, gene catalogue and annotation | [`03_metagenome/`](03_metagenome/) |
| 04 | CRISPRone spacer-supported putative host associations, CAT_pack taxonomy, and downstream host characterization | [`04_virus_host/`](04_virus_host/) |

The dependency order is:

```text
01 DNA virome ----------------┐
                             ├──> 04 Virus-host associations
03 Prokaryotic metagenome ----┘

02 RNA virome -> independent RNA-vOTU analysis
```

Shared downstream R scripts remain in `scripts/r/` to avoid duplicating the
same TPM, diversity, and statistical functions in multiple modules.
