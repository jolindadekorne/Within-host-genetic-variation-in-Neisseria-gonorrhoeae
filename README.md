# Within-host-genetic-variation-in-Neisseria-gonorrhoeae
This repository contains the Snakemake pipeline used to assess within-host genetic variation in _Neisseria gonorrhoeae_. Isolate were obtained during a clinical trial on novel antibiotic treatment options for gonorrhoea. Isolate pairs were formed that consisted of isolates obtained from a single individual from multiple time points (reflecting antibiotic treatment failure) and from multiple anatomical locations. 

## Dependencies
This pipeline uses the following dependencies:
- Conda
- Snakemake
- Python3

## Input
This pipeline uses forward and reverse raw Illumina sequencing reads which are located in the folder `raw_data`. The raw data files should be named `{id}_R1.fastq.gz` and `{id}_R2.fastq.gz`. 

## Other files needed
- [maskrc-svg script](https://github.com/kwongj/maskrc-svg); expected path: 'scripts/maskrc-svg.py'
- Reference genome FA1090 is used for calculating coverage and calling variants: NC_002946.2. The reference genome should be located in the same directory as the Snakefile.

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
| Mask the recombination sites in the core genome alignment | [maskrc-svg script](https://github.com/kwongj/maskrc-svg) |
| Calculate recombination filtered- and unfiltered SNP distances | snp-dists |

## Visualization of SNPs using Artemis
The script that was used to create input files for visualization of SNPs in Artemis is also added to this repository, together with the files needed to run this script:
 - List containing all isolate pairs to compare
 - Annotation files of reference FA1090
 - List of genes annotated in FA1090

The script uses the core genome alignment, created by Snippy in the pipeline: 'snippy_refFA1090_out/clean.full.aln'




