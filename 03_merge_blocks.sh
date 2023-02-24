#!/bin/bash
# ------------------------------------------------------------------
# Author:       A.Sendell-Price
# Date:         February 2023
# Title:        xx_merge_blocks.sh
# Description:  Merges individual "genomic block" VCFs into a single
#               file.
# slurm:        #SBATCH -p core -n 1
#               #SBATCH -t 2-00:00:00
# ------------------------------------------------------------------

# Load required modules
module load bioinfo-tools GATK/4.2.0.0 bcftools

# Merge genomic block VCFs into a single VCF
cd joint_genotyping
ls block_*[0-9].vcf.gz > block.vcf.list
bcftools concat -f block.vcf.list -o raw_merged.vcf.gz -O z
gatk IndexFeatureFile -I raw_merged.vcf.gz
rm block_*[0-9].vcf.gz*
