#!/bin/bash
# ------------------------------------------------------------------
# Author:       A.Sendell-Price
# Date:         February 2023
# Title:        ientify_high_conf_heterozygotes.sh
# Description:  Identifies sites where both parents are homozygous for the same
#               allele and offspring are all heteozygous.
# slurm:        #SBATCH -p core -n 1
#               #SBATCH -t 2-00:00:00
# ------------------------------------------------------------------

# Load required modules
module load bioinfo-tools GATK/4.2.0.0 bcftools BEDTools/2.29.2 samtools/1.5 vcftools/0.1.16

# Make directory and move into it
mkdir high_conf_heterozygotes
cd high_conf_heterozygotes

# Extract sample names from VCF
bcftools query -l \
  ../hard_filtering/called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
  > all.samples

# Set sample genotype to missing if allele balance is:
# < 0.2 for heterozygote calls
# < 0.9 for homozygote calls 
# This is implemented using the jvarkit tool "vcffilterjdk"
# Remove any sites where parental genotypes are now missing
source ../accessory_scripts/filter_het_homo_by_AB.sh \
  ../hard_filtering/called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
  temp.vcf.gz
zcat temp.vcf.gz | bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
  | bgzip > called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.vcf.gz
rm temp.vcf.gz

# Extract positions where parents are homozygous for different alleles
# and all offspring are heterozygous and all sample genotypes are non-missing.
zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.vcf.gz \
  | bcftools filter -e 'GT[@all.samples]="mis"' \
  | bcftools filter -i 'GT[0]="het"' \
  | bcftools filter -i 'GT[1]="het"' \
  | bcftools filter -i 'GT[2]="het"' \
  | bcftools filter -i 'GT[3]="het"' \
  | bcftools filter -i 'GT[4]="het"' \
  | bcftools filter -i 'GT[5]="het"' \
  | bcftools filter -i 'GT[6]="het"' \
  | bcftools filter -i 'GT[7]="het"' \
  | bcftools filter -i 'GT[8]="het"' \
  | bcftools view -i 'GT[9]="ref" & GT[9]="hom"' \
  | bcftools view -i 'GT[10]="alt" & GT[10]="hom"' \
  | bcftools sort | bgzip > temp1.vcf.gz

zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.vcf.gz \
  | bcftools filter -e 'GT[@all.samples]="mis"' \
  | bcftools filter -i 'GT[0]="het"' \
  | bcftools filter -i 'GT[1]="het"' \
  | bcftools filter -i 'GT[2]="het"' \
  | bcftools filter -i 'GT[3]="het"' \
  | bcftools filter -i 'GT[4]="het"' \
  | bcftools filter -i 'GT[5]="het"' \
  | bcftools filter -i 'GT[6]="het"' \
  | bcftools filter -i 'GT[7]="het"' \
  | bcftools filter -i 'GT[8]="het"' \
  | bcftools view -i 'GT[9]="alt" & GT[9]="hom"' \
  | bcftools view -i 'GT[10]="ref" & GT[10]="hom"' \
  | bcftools sort | bgzip > temp2.vcf.gz

bcftools index temp1.vcf.gz
bcftools index temp2.vcf.gz

java -jar ../BIN/picard.jar MergeVcfs \
    I=temp1.vcf.gz \
    I=temp2.vcf.gz \
    O=HighConf_heterozygous_sites.vcf.gz

#Remove temp files
rm temp1.vcf.gz* temp2.vcf.gz*

# END!
