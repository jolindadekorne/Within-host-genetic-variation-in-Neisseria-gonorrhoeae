# Within-host-genetic-variation-in-Neisseria-gonorrhoeae
This repository contains the Snakemake pipeline used to assess within-host genetic variation in _Neisseria gonorrhoeae_. Isolate were obtained during a clinical trial on novel antibiotic treatment options for gonorrhoea. Isolate pairs were formed that consisted of isolates obtained from a single individual from multiple time points (reflecting antibiotic treatment failure) and from multiple anatomical locations. 

## Input
This pipeline uses forward and reverse raw Illumina sequencing reads which are located in the folder `raw_data`. The raw data files should be named `{id}_R1.fastq.gz` and `{id}_R2.fastq.gz`

## Pipeline 
The pipeline includes the following steps and tools:

| Step     | Tool     |   
| ---------|----------|
| Filter low quality raw reads + trim adapters | fastp |
| Assembly | skesa | 
| Assembly quality check | QUAST |
| Read mapping against reference genome FA1090 | BWA-MEM2 |
| Calculation of percentage of bases covered and coverage depth | samtools |
| Call variants using reference FA1090 + create core genome alignment | snippy |
| Remove recombination from variant alignment | Gubbins |
| Mask the recombination sites in the core genome alignment | maskrc-svg script |
| Calculate recombination filtered- and unfiltered SNP distances | snp-dists |
!Add script for visualization of SNPs??

## Other files needed
- Reference genome FA1090 is used for calculating coverage and calling variants: NC_002946.2. The reference genome should be located in the same directory as the Snakefile.


