#!/bin/bash
# ------------------------------------------------------------------
# Author:       A.Sendell-Price
# Date:         February 2023
# Title:        07_conduct_secondary_genotyping.sh
# Description:  ADD
# slurm:        #SBATCH -p core -n 1
#               #SBATCH -t 2-00:00:00
# Requirements: requires installation of Simon Martin's
#               genomics general scripts:
#               https://github.com/simonhmartin/genomics_general
# ------------------------------------------------------------------

# Load required modules
module load bioinfo-tools bcftools samtools/1.5

# Make directory and move into it
mkdir secondary_genotyping
cd secondary_genotyping

# Make subdirectories
mkdir bams
mkdir mpileup_VCFs

# Extract mutation positions from VCF
bcftools query -f '%CHROM\t%POS\n' \
  ../putative_mutations_GATK/putative_mutations.vcf.gz \
  > mutation.positions

#Extract sample IDs
bcftools query -l ../putative_mutations_GATK/putative_mutations.vcf.gz \
  > sample.list

# For each mutation do the following:
for MUTATION_NO in $(seq 1 $(wc -l < mutation.positions))
do
    SCAFF=$(head -n $MUTATION_NO mutation.positions | tail -n 1 | cut -f 1)
    POS=$(head -n $MUTATION_NO mutation.positions | tail -n 1 | cut -f 2)
    UPSTREAM=$(expr $POS - 5000)
    DOWNSTREAM=$(expr $POS + 5000)

    for SAMPLE in $(cat sample.list)
    do
        #Path to sample bam file list
        BAMS=/proj/snic2020-2-19/private/shark/users/mats/genotyping/${SAMPLE}_bam.list

        #Count number of BAM files in list
        BAM_COUNT=$(wc -l < $BAMS)

        #For each bam extract reads within region of interest
        for BAM_NO in $(seq 1 $BAM_COUNT)
        do
            samtools view -h $(head -n $BAM_NO $BAMS | tail -n 1 | cut -d " " -f 2) \
            "${SCAFF}:${UPSTREAM}-${DOWNSTREAM}" > temp_${BAM_NO}.bam
        done

        #Merge temp bam files
        ls temp_*.bam > temp.bam.list
        samtools merge -b temp.bam.list merged.bam -p

        #Fix sample names in header
        samtools view -H merged.bam | sed "s/SM:[^\t]*/SM:$SAMPLE/g" \
        | samtools reheader - merged.bam \
        > bams/${SCAFF}_${POS}_${SAMPLE}.bam

        #Tidy up
        rm merged.bam temp*
    done
    
    #Create list of sample bams for region of interest
    ls bams/${SCAFF}_${POS}_*.bam > bam.list

    #Re-call snps within region of interest using bcftools mpileup
    bcftools mpileup -Ou -f $REF -b bam.list -q 58 -Q 20 | bcftools call -m -Ov -o mpileup_VCFs/${SCAFF}_${POS}.vcf

    #Tidy up
    rm bam.list
done


##########################################################################
##########################################################################

#Pull out mutations from VCFs and comibine into a single new vcf file
#Extract header lines
cat $(ls mpileup_VCFs/scaffold_* | head -n 1) | grep "#" > putative_mutations.mpileup.vcf

#Extract mutation mpileup genotypes if present
for MUTATION_NO in $(seq 1 $(wc -l < mutation.positions))
do  
    SCAFF=$(head -n $MUTATION_NO mutation.positions | tail -n 1 | cut -f 1)
    POS=$(head -n $MUTATION_NO mutation.positions | tail -n 1 | cut -f 2)
    cat mpileup_VCFs/${SCAFF}_${POS}.vcf | grep -v "#" | grep ${POS} >> putative_mutations.mpileup.vcf
    echo "Processing complete for mutation no." $MUTATION_NO
done

#Convert to geno format
python ../genomics_general/VCF_processing/parseVCF.py \
-i putative_mutations.mpileup.vcf > putative_mutations.mpileup.genotypes.txt

# END!
