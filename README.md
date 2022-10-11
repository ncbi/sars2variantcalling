# The NCBI SRA-B (Sequence Read Archive Bioinformatics)

### Contact:
email: sra@ncbi.nlm.nih.gov

### Change Log
Please check the [CHANGES.md](CHANGES.md) file for change history.

----

## The SRA-B SARS COV2 virus variations pipeline 
The NCBI SARS COV2 virus variations pipeline is a collection of [snakemake](https://snakemake.readthedocs.io/en/stable/) workflows.
See their overview [here](https://www.ncbi.nlm.nih.gov/sra/docs/sars-cov-2-variant-calling/).

Results are available via following NCBI resources:
1) [NCBI SARS2-COV-2](https://www.ncbi.nlm.nih.gov/sars-cov-2/)
2) [Weekly ACTIV-TRACE reports](https://ftp.ncbi.nlm.nih.gov/pub/ACTIV-TRACE/) 

#### Content:
 * Wrapper script to run workflow: 
   * toolbox/workflow/run.sh
 * Reference sequence and snpEff databases:
   * toolbox/static/reference/
   * toolbox/static/snpEff/
   
#### Usage:
```bash
  run.sh --platform --accession|--list [--instrument] [--conf] [--workdir] [--help]
    --platform:   platform, choices = [illumina, ont, pacbio, genbank], default = illumina
    --instrument: applicable to ONT platform only, use with single accession option, default = PromethION
    --accession:  single accession, optional, either --accession or --list must be specified
    --list:       file with accession list, optional, either --accession or --list must be specified
                  NOTE: ONT file list is expected to have at least two columns: <acc> <instrument>
    --conf:       optional custom configfile to override default
    --workdir:    optional working directory, default = ./workdir
    --help:       to display this help
``` 

#### Usage examples:
1) ILLUMINA 
   ```bash
   toolbox/workflow/run.sh --accession SRR21830388
   # in case of using docker image
   /pipelines/toolbox/workflow/run.sh --accession SRR21830388 --conf /pipelines/toolbox/workflow/extra.config.yaml 
   ```
2) OXFORD-NANOPORE
   ```bash
   toolbox/workflow/run.sh SRR15965069 --platform ont --instrument GridION
   ```
3) PacBio
   ```bash
   toolbox/workflow/run.sh SRR14895419 --platform pacbio
   ```
