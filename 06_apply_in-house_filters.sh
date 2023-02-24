#!/bin/bash
# ------------------------------------------------------------------
# Author:       A.Sendell-Price
# Date:         February 2023
# Title:        06_apply_in-house_filters.sh
# Description:  ADD
# slurm:        #SBATCH -p core -n 1
#               #SBATCH -t 2-00:00:00
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
  ../high_conf_heterozygotes/called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.vcf.gz \
  called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.PASSED_inhouse_filters.vcf.gz \
  11

# Apply filtering criteria by running the script created by "generate.custom.filters.R"
source apply.custom.filters.sh

# From the custom filtered dataset, extract positons where both parents are
# homozygous ref and at least one offspring is hereozygous
zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.PASSED_inhouse_filters.vcf.gz \
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

#Combine the two VCFs into a single file, these represent putative mutations
java -jar /proj/snic2020-2-19/private/shark/users/ash/BIN/picard/build/libs/picard.jar MergeVcfs \
  I=temp1.vcf.gz \
  I=temp2.vcf.gz \
  O=putative_mutations.vcf.gz

#Remove temp files
rm temp1.vcf.gz* temp2.vcf.gz*

# END!
