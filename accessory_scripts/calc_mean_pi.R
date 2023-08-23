library(dplyr)

#Load windowed pi estimates for each parental sample
pixy_pi_sHemOce2 <- read.delim("pixy_pi.txt") %>% filter(pop == 1) %>% filter(no_sites >= 5000) 
pixy_pi_sHemOce3 <- read.delim("pixy_pi.txt") %>% filter(pop == 2) %>% filter(no_sites >= 5000) 

#Calculate mean and SE across windows
se <- function(x) sqrt(var(x) / length(x))
mean(na.omit(pixy_pi_sHemOce2$avg_pi)) #0.001802193
mean(na.omit(pixy_pi_sHemOce3$avg_pi)) #0.002132612
se(na.omit(pixy_pi_sHemOce2$avg_pi)) #1.177419e-05
se(na.omit(pixy_pi_sHemOce3$avg_pi)) #1.528566e-05

#Do all the above but combined
pixy_pi <- read.delim("pixy_pi.txt") %>% filter(no_sites >= 5000) 
mean(na.omit(pixy_pi$avg_pi)) #0.001967192
se(na.omit(pixy_pi$avg_pi)) #9.66323e-06
