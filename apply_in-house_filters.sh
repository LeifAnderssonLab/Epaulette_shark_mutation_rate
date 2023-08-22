#!/bin/bash
# ------------------------------------------------------------------
# Author:       A.Sendell-Price
# Date:         February 2023
# Title:        06_apply_in-house_filters.sh
# Description:  Applies in house genotype filtering criteria to dataset
# slurm:        #SBATCH -p core -n 1
#               #SBATCH -t 2-00:00:00
# Requirements: requires installation of Simon Martin's
#               genomics general scripts:
#               https://github.com/simonhmartin/genomics_general
# ------------------------------------------------------------------

# Load required modules
module load bioinfo-tools GATK/4.2.0.0 bcftools BEDTools/2.29.2 samtools/1.5 vcftools/0.1.16

# Make directory and move into it
mkdir putative_mutations_GATK
cd putative_mutations_GATK

# Based on the distribution of genotype quality annotations at
# high confidence heterozygous sites, define in-house quality
# filtering criteria using 5th and 95th percentiles as cut-offs. 
# mapping quality: MQ (lower bound only)
# mapping quality rank sum: MQRankSum
# base quality rank sum: BaseQRankSum
# read position rank sum: ReadPosRankSum
# quality by depth: QD (lower bound only)
# Filtering will then be applied to the following individual annotations:
# genotype quality: GQ (lower bound only)
# sample read depth: DP
Rscript ../accessory_scripts/generate.custom.filters.R \
  ../high_conf_heterozygotes/HighConf_heterozygous_sites.vcf.gz \
  ../hard_filtering/called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
  called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.PASSED_inhouse_filters.vcf.gz \
  11

# Apply filtering criteria by running the script created by "generate.custom.filters.R"
source apply.custom.filters.sh
# This file looks like this:
# zcat ../hard_filtering/called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
#   | bcftools view -e 'INFO/BaseQRankSum < -0.908 || INFO/BaseQRankSum > 1.07' \
#   | bcftools view -e 'INFO/ReadPosRankSum < -0.728 || INFO/ReadPosRankSum > 0.861' \
#   | bcftools view -e 'INFO/MQ < 58.19' \
#   | bcftools view -e 'INFO/QD < 20.75' \
#   | bcftools filter -S . -e 'FMT/DP[0] < 28 || FMT/DP[0] > 76' \
#   | bcftools filter -S . -e 'FMT/DP[1] < 24 || FMT/DP[1] > 67' \
#   | bcftools filter -S . -e 'FMT/DP[2] < 14 || FMT/DP[2] > 45' \
#   | bcftools filter -S . -e 'FMT/DP[3] < 27 || FMT/DP[3] > 75' \
#   | bcftools filter -S . -e 'FMT/DP[4] < 14 || FMT/DP[4] > 48' \
#   | bcftools filter -S . -e 'FMT/DP[5] < 17 || FMT/DP[5] > 52' \
#   | bcftools filter -S . -e 'FMT/DP[6] < 19 || FMT/DP[6] > 56' \
#   | bcftools filter -S . -e 'FMT/DP[7] < 22 || FMT/DP[7] > 62' \
#   | bcftools filter -S . -e 'FMT/DP[8] < 17 || FMT/DP[8] > 50' \
#   | bcftools filter -S . -e 'FMT/DP[9] < 45 || FMT/DP[9] > 99' \
#   | bcftools filter -S . -e 'FMT/DP[10] < 101 || FMT/DP[10] > 159' \
#   | bgzip > called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.PASSED_inhouse_filters.vcf.gz

# From the custom filtered dataset, extract positons where both parents are
# homozygous ref and at least one offspring is hereozygous
zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.PASSED_inhouse_filters.vcf.gz \
  | bcftools view -i 'GT[9]="ref" & GT[9]="hom"' \
  | bcftools view -i 'GT[10]="ref" & GT[10]="hom"' \
  | bcftools view -i 'GT[0]="het" || GT[1]="het" || GT[2]="het" || GT[3]="het" || GT[4]="het" || GT[5]="het" || GT[6]="het" || GT[7]="het" || GT[8]="het"' \
  | bcftools sort | bgzip > temp1.vcf.gz

#From the custom filtered dataset, extract positons where both parents are homozygous
# homozygous alt and at least one offspring is hereozygous
  zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.PASSED_inhouse_filters.vcf.gz \
  | bcftools view -i 'GT[9]="alt" & GT[9]="hom"' \
  | bcftools view -i 'GT[10]="alt" & GT[10]="hom"' \
  | bcftools view -i 'GT[0]="het" || GT[1]="het" || GT[2]="het" || GT[3]="het" || GT[4]="het" || GT[5]="het" || GT[6]="het" || GT[7]="het" || GT[8]="het"' \
  | bcftools sort | bgzip > temp2.vcf.gz

#Combine the two VCFs into a single file
java -jar /proj/snic2020-2-19/private/shark/users/ash/BIN/picard/build/libs/picard.jar MergeVcfs \
  I=temp1.vcf.gz \
  I=temp2.vcf.gz \
  O=putative_mutations.vcf.gz
rm temp1.vcf.gz* temp2.vcf.gz*

#For putative mutations, set sample genotype to missing if allele balance is:
#< 0.2 for heterozygote calls
#< 0.9 for homozygote calls 
#This is implemented using the jvarkit tool "vcffilterjdk"
source ../accessory_scripts/filter_het_homo_by_AB.sh \
  putative_mutations.vcf.gz \
  putative_mutations.AB_filtered.vcf.gz

#Convert AB-filterd putative mutations to "genotype" format i.e. A,T,C,G
python ../genomics_general/VCF_processing/parseVCF.py \
  -i putative_mutations.AB_filtered.vcf.gz > putative_mutations.AB_filtered.GATK_genotypes.txt
  
# END!
