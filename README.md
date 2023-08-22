# Estimation of germline mutation rate in the Epaulette shark

This repository contains scripts necessary to re-conduct the de novo mutation rate estimation pipeline presented in our manuscript:
Slow evolution of the epaulette shark illustrated by its low rate of de novo mutations.

## Overview
The pipeline was run on Linux computing resources at Uppsala Multidisciplinary Center for Advanced Computational Science (UPPMAX) with jobs submitted using a slurm job scheduling system. Scripts may therefore require some editing to run on a different computing resource. Nonetheless, scripts should be run in the following order:

1. Read processing (XXXXX): Filter fastq files to remove ....
2. Read mapping (): Map each lane to the reference genome (bwa2 mem) producing multiple bam files per offspring.
3. Individual genotyping (block_individual_genotyping.sh): Call genotypes per individual with GATK HaplotypeCaller in BP-RESOLUTION. Genotypes calling is split into ten genomic blocks for speed. Genomic block ranges are defined in file XXXXXX.
4. Joint genotype calling (block_joint_genotyping.sh): Combine individual sample gVCFs with GATK CombineGVCFs and conduct joint genotyping with GATK GenotypeGVCF. Also split into ten genomic blocks for speed.
5. Merge genomic block VCFs (merge_blocks.sh): Merge the ten genomic block VCFs into a single VCF file.
6. Apply hard filters (apply_hard_filters.sh): Extract biallelic and monomorphic sites, remove genotypes with GQ < 20, remove sites within repeat reions or regions of low Mappability, remove sites with missing parental genotypes.
7. Define in-house filters (identify_high_conf_heterozygotes.sh): Identify high confidence heterozygous sites and use these to define quality filter cut-offs. 
8. Apply in-house filters (apply_in-house_filters.sh): Apply in-house filters and identify putative de novo mutations. 
9. Apply secondary genotype calling (secondary_genotyping.sh): Perform secondary genotype calling at putative mutations using bcftools mpileup and call. Only putative mutations where GATK and bcftools genotypes match are considered true candidate de novo mutations. 

## Requirements 

Default software used to run the pipeline are:
- bwa2 vX.X
- samtools
- bcftools
- picard tools
