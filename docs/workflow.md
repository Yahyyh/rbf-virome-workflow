# Workflow architecture

The workflow is documented in four independent modules:

1. [`DNA virome`](../workflow/01_dna_virome/)
2. [`RNA virome`](../workflow/02_rna_virome/)
3. [`Prokaryotic metagenome`](../workflow/03_metagenome/)
4. [`Putative virus-host associations`](../workflow/04_virus_host/)

```text
DNA metavirome reads                         RNA metavirome reads
        |                                             |
        v                                             v
01 DNA virome                                  02 RNA virome
        |
        | DNA vOTU representatives
        v
04 Putative virus-host associations <--- 03 Prokaryotic metagenome
                                            |
                                            └── host contigs and sample mapping
```
