#!/usr/bin/env Rscript
options(warn=-1)

args = commandArgs(trailingOnly=TRUE)

#Extract variables from bash submission
VCF1=args[1]
VCF2=args[2]
OUT=args[3]
NO_SAMPLES=args[4]

#Initiate filter command expression
filter.cmd <- "#!/bin/bash"
tmp.cmd <- paste0("zcat ", VCF2, " \\")
filter.cmd <- append(filter.cmd, tmp.cmd, after = length(filter.cmd))

###############################################################################################
# SITE FILTERING CRITERIA
###############################################################################################

#Estimate filtering bounds for various "all sample" annotations
message(paste0("Estimating site filtering ranges for ..."))
#filter by upper and lower bounds
for (STATISTIC in c("BaseQRankSum","ReadPosRankSum")){
  message(paste0("  - ", STATISTIC))
  #Extract statistic of interest from VCF file
  system(paste0("bcftools query -f '%",STATISTIC,"\n' ",VCF1," > stats.tmp"))
  #Load statistics
  results <- na.omit(as.numeric(readLines("stats.tmp")))
  #calculate 5th and 95th quartile
  results.cutoffs <- quantile(results, c(0.05, 0.95))
  #Plot histogram
  png(file = paste0(STATISTIC,".png"))
  hist(results, col="lightblue")
  abline(v = results.cutoffs[[1]], col="red", lwd=3, lty=2)
  abline(v = results.cutoffs[[2]], col="red", lwd=3, lty=2)
  dev.off()
  #Append expression to command
  tmp.cmd <- paste0("| bcftools view -e 'INFO/",STATISTIC," < ",results.cutoffs[[1]]," || INFO/",STATISTIC," > ",results.cutoffs[[2]],"' \\")
  filter.cmd <- append(filter.cmd, tmp.cmd, after = length(filter.cmd))
  write.table(results, paste0(STATISTIC,".txt"))
}

#Filter by lower bound only
for (STATISTIC in c("MQ","QD")){
  message(paste0("  - ", STATISTIC))
  #Extract statistic of interest from VCF file
  system(paste0("bcftools query -f '%",STATISTIC,"\n' ",VCF1," > stats.tmp"))
  #Load statistics
  results <- na.omit(as.numeric(readLines("stats.tmp")))
  #calculate 5th and 95th quartile
  results.cutoffs <- quantile(results, c(0.05, 0.95))
  #Plot histogram
  png(file = paste0(STATISTIC,".png"))
  hist(results, col="lightblue")
  abline(v = results.cutoffs[[1]], col="red", lwd=3, lty=2)
  dev.off()
  #Append expression to command
  tmp.cmd <- paste0("| bcftools view -e 'INFO/",STATISTIC," < ",results.cutoffs[[1]],"' \\")
  filter.cmd <- append(filter.cmd, tmp.cmd, after = length(filter.cmd))
  write.table(results, paste0(STATISTIC,".txt"))
}

###############################################################################################
# INDIVIDUAL GENOTYPE FILTERING CRITERIA
###############################################################################################

message(paste0("Estimating individual DP filtering ranges for ..."))
for (SAMPLE in 1:NO_SAMPLES){
  message(paste0("  -  individual ", SAMPLE))
  #####################################################################
  # DP filtering
  #####################################################################
  #Extract field from VCF file
  system(paste0("bcftools query -f '[%DP,]\n' ",VCF1," | cut -d ',' -f ", SAMPLE," > stats.tmp"))
  #Load statistics
  results <- na.omit(as.numeric(readLines("stats.tmp")))
  #calculate 5th and 95th quartile
  results.cutoffs <- quantile(results, c(0.05, 0.95))
  #Append expression to command
  tmp.cmd <- paste0("| bcftools filter -e 'FMT/DP[",SAMPLE-1,"] < ",results.cutoffs[[1]]," || FMT/DP[",SAMPLE-1,"] > ",results.cutoffs[[2]],"' \\")
  filter.cmd <- append(filter.cmd, tmp.cmd, after = length(filter.cmd))
  #Plot histogram
  png(file = paste0("DP_",SAMPLE,".png"))
  hist(results, col="lightblue")
  abline(v = results.cutoffs[[1]], col="red", lwd=3, lty=2)
  abline(v = results.cutoffs[[2]], col="red", lwd=3, lty=2)
  dev.off()
  write.table(results, paste0("DP_",SAMPLE,".txt"))
}

#Add output line to command
tmp.cmd <- paste0("| bgzip > ", OUT)
filter.cmd <- append(filter.cmd, tmp.cmd, after = length(filter.cmd))

#Write filtering script to file
write(filter.cmd, "apply.custom.filters.sh")
message("Completed! Apply filtering using command 'source apply.custom.filters.sh'")
