# Estimation of germline mutation rate in the Epaulette shark

This repository contains scripts necessary to perform the in-house genotype filtering pipeline detailed in our manuscript:
Extremely low de novo mutation rate in the epaulette shark. This pipeline was used to identify candidate de novo mutations within a shark pedigree.

## Overview
The pipeline was run on Linux computing resources at Uppsala Multidisciplinary Center for Advanced Computational Science (UPPMAX) with jobs submitted using a slurm job scheduling system. This script requires a single GATK BP-RESOLUTION VCF file as input. Scripts may require some editing to run on a different computing resource. Nonetheless, scripts should be run in the following order:

1. Apply hard filters (apply_hard_filters.sh): Extract biallelic and monomorphic sites, remove genotypes with GQ < 20, remove sites within repeat reions or regions of low Mappability, remove sites with missing parental genotypes.
2. Identify high confidence heterozygous sites (identify_high_conf_heterozygotes.sh): Identify high confidence heterozygous sites to be used to define quality filter cut-offs.
3. Apply in-house filters (apply_in-house_filters.sh): Define and apply in-house filters and identify putative de novo mutations. 
4. Apply secondary genotype calling (secondary_genotyping.sh): Perform secondary genotype calling at putative mutations using bcftools mpileup and call. Only putative mutations where GATK and bcftools genotypes match are considered true candidate de novo mutations. 

## Requirements 

Default software used to run the pipeline are:
- samtools v1.5
- vcftools v0.1.16
- bcftools v1.17
- bedtools v2.29.2
- picard tools v3.1.0

## Additional script
This directory "accessory_scripts" also contains the sript "calc_mean_pi.R" which was used to calculate mean nucleotide diversity presented in Fig. 7. 
