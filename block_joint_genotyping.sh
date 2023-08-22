#!/bin/bash
# ------------------------------------------------------------------
# Author:       A.Sendell-Price
# Date:         February 2023
# Title:        block_joint_genotyping.sh
# Description:  Performs joint sample genotyping of parental and
#               offspring samples for a given "genomic block".
# slurm:        #SBATCH -p core -n 1
#               #SBATCH -t 2-00:00:00
#               #SBATCH --array=1-10:1
# ------------------------------------------------------------------

# Load required modules
module load bioinfo-tools GATK/4.2.0.0 bcftools

# Create / move into directory "joint_genotyping"
if [[ ! -d joint_genotyping ]]
then
    mkdir joint_genotyping
    cd joint_genotyping
else
    cd joint_genotyping
fi

# Set block number using slurm array task ID
BLOCK_NO=$SLURM_ARRAY_TASK_ID

# Combine GVCFs
gatk --java-options "-Xmx120g" CombineGVCFs \
-R /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.fasta \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/ind1722_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/ind1835_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/ind1895_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/ind1923_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/ind1925_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/ind1983_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/ind2023_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/ind2024_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/ind2046_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/sHemOce2_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/sHemOce3_block_${BLOCK_NO}.g.vcf.gz \
--convert-to-base-pair-resolution \
-O block_${BLOCK_NO}.g.vcf.gz

# Conduct joint genotyping of gVCF
gatk --java-options "-Xmx120g" GenotypeGVCFs \
--include-non-variant-sites \
-R /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.fasta \
-V block_${BLOCK_NO}.g.vcf.gz \
-O block_${BLOCK_NO}.vcf.gz

# Remove gVCF and its index as no longer needed
rm block_${BLOCK_NO}.g.vcf.gz*

# END!
