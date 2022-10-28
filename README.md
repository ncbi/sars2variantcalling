# The NCBI SRA-B (Sequence Read Archive Bioinformatics)

### Contact:
email: sra@ncbi.nlm.nih.gov

----

## The SRA-B SARS COV2 virus variations pipeline 
The NCBI SARS COV2 virus variations pipeline is a collection of [snakemake](https://snakemake.readthedocs.io/en/stable/) workflows.
See their overview [here](https://www.ncbi.nlm.nih.gov/sra/docs/sars-cov-2-variant-calling/).

Results are available via following NCBI resources:
1) [NCBI SARS2-COV-2](https://www.ncbi.nlm.nih.gov/sars-cov-2/)
2) [Weekly ACTIV-TRACE reports](https://ftp.ncbi.nlm.nih.gov/pub/ACTIV-TRACE/) 

#### Content:
 * Wrapper script to run workflow: 
   * [run.sh](toolbox/workflow/run.sh)
 * Reference sequence and snpEff databases:
   * toolbox/static/reference/
   * toolbox/static/snpEff/
   
#### Dependencies:
Following third party tools assumed to be installed in user environment:
* [samtools](https://github.com/samtools):
  * htslib
  * samtools
  * bcftools 
* [mummer](https://github.com/mummer4/mummer)
* [snpEff](https://sourceforge.net/projects/snpeff/files/snpEff_v4_2_core.zip/download)
* [trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip)
* [minimap2](https://github.com/lh3/minimap2)
* [SRA-toolkit](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/)
* [Picard](https://github.com/broadinstitute/picard)
* [GATK](https://github.com/broadinstitute/gatk)
* [medaka](https://github.com/nanoporetech/medaka)
* hisat2
* [vcftools](https://sourceforge.net/projects/vcftools/files/vcftools_0.1.13.tar.gz/download)

For the comprehensive list of tools used and their versions, please refer to [Dockerfile](toolbox/Dockerfile) 

#### Reference sequence indexing 
Pipeline expect reference sequence to be indexed, to create these indexes please run
```bash
gatk CreateSequenceDictionary -R toolbox/static/reference/NC_045512.2.fa
samtools faidx toolbox/static/reference/NC_045512.2.fa
gatk IndexFeatureFile --input toolbox/static/reference/NC_045512.2.known_sites.vcf
hisat2-build toolbox/static/reference/NC_045512.2.fa toolbox/static/reference/NC_045512.2.fa
```
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
   toolbox/workflow/run.sh --accession SRR15965069 --platform ont --instrument GridION
   ```
3) PacBio
   ```bash
   toolbox/workflow/run.sh --accession SRR14895419 --platform pacbio
   ```
