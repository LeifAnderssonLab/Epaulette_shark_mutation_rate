# Estimation of germline mutation rate in the Epaulette shark

This repository contains scripts necessary to re-conduct the de novo mutation rate estimation pipeline presented in our manuscript:
Pedigree-based estimate of de novo mutation rate in the epaulette shark defines the lowest known vertebrate mutation rate.

## Overview

The pipeline was design for a slurm job scheduling system on Linux. Using another job scheduling system would require some changes in the scripts. Scripts should be run in the following order:

1. Read processing (01_XXXXX): Filter fastq files to remove ....
2. Read mapping (): Map each lane to the reference genome (bwa2 mem) producing multiple bam files per offspring.
3. Individual genotyping (XX_block_individual_genotyping.sh): Call genotypes per individual with GATK HaplotypeCaller in BP-RESOLUTION. Genotypes calling is split into ten genomic blocks for speed. Genomic block ranges are defined in file XXXXXX.
4. Joint genotype calling (XX_block_joint_genotyping.sh): Combine individual sample gVCFs with GATK CombineGVCFs and conduct joint genotyping with GATK GenotypeGVCF. Also split into ten genomic blocks for speed.
5. Merge genomic block VCFs (XX_merge_blocks.sh): Merge the ten genomic block VCFs into a single VCF file.
6. Apply hard filters (04_apply_hard_filters.sh): Extract biallelic and monomorphic sites, remove genotypes with GQ < 20, remove sites within repeat reions or regions of low Mappability, remove sites with missing parental genotypes.
7. Define in-house filters (XX_identify_high_conf_heterozygotes.sh): Identify high confidence heterozygous sites and use these to define quality filter cut-offs. 
8. Apply in-house filters (XX_apply_in-house_filters.sh): Apply in-house filters and identify putative de novo mutations. 
9. Apply secondary genotype calling (XX_secondary_genotyping.sh): Perform secondary genotype calling at putative mutations using bcftools mpileup and call. Only putative mutations where GATK and bcftools genotypes match are considered true candidate de novo mutations. 

## Requirements 

Default software used to run the pipeline are:
- bwa2 vX.X
- samtools
- bcftools
- picard tools
