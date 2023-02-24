#!/bin/bash
# ------------------------------------------------------------------
# Author:       A.Sendell-Price
# Date:         February 2023
# Title:        04_apply_hard_filters.sh
# Description:  ADD
# slurm:        #SBATCH -p core -n 1
#               #SBATCH -t 2-00:00:00
# ------------------------------------------------------------------

# Load required modules
module load bioinfo-tools GATK/4.2.0.0 bcftools BEDTools/2.29.2 samtools/1.5 vcftools/0.1.16

# Make directory and move into it
mkdir hard_filtering
cd hard_filtering

# From raw VCF file extract invariant sites
# plus remove any site that doesnt pass filtering (i.e. contains LowQual in filter column)
gatk SelectVariants \
  -R ../assembly/sHemOce1.mat.decon.20210528.fasta \
  -V ../joint_genotyping/raw_merged.vcf.gz \
  --select-type-to-include NO_VARIATION \
  -O called_by_GATK_invariant.vcf.gz

# From raw VCF file extract biallelic sites
# plus remove any site that doesnt pass filtering (i.e. contains LowQual in filter column)
gatk SelectVariants \
  -R ../assembly/sHemOce1.mat.decon.20210528.fasta \
  -V ../joint_genotyping/raw_merged.vcf.gz \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  -O called_by_GATK_biallelic.vcf.gz

# Combine invariant and biallelic VCFs using 
# picard tools (https://broadinstitute.github.io/picard/)
java -jar ../BIN/picard.jar MergeVcfs \
  I=called_by_GATK_invariant.vcf.gz \
  I=called_by_GATK_biallelic.vcf.gz \
  O=called_by_GATK_invariant_plus_biallelic.vcf.gz

# Remove sites where parental genotypes are missing as these aren't informative.
zcat called_by_GATK_invariant_plus_biallelic.vcf.gz \
  | bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
  | bgzip > called_by_GATK_invariant_plus_biallelic.informative_sites.vcf.gz

# Remove sites from repeat regions / regions of low mappability
bedtools intersect -a called_by_GATK_invariant_plus_biallelic.informative_sites.vcf.gz \
  -b ../resources/HiglyMappable_Unmasked_ranges.txt -header \
  | bgzip > called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.vcf.gz

# Set sample genotypes with GQ < 20 to missing "./." and again 
# remove any sites where parental genotypes are missing
zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.vcf.gz \
  | bcftools filter -S . -e 'FMT/GQ[@all.samples] < 20' \
  | bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
  | bgzip > called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.vcf.gz

# END!
