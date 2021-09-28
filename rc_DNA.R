# quantify assembly processes for RNA dataset
# Xiu Jia at University of Groningen built on 28-01-2019
# updated on 19-11-2019

# load packages --------------------------------------------------------------------------------------------
rm(list=ls())

library(vegan)

source("raup_crick_abundance.R")

SUBSET = 'DNA'

# read in com table and convert to a data.frame
print("read table")
com <- read.csv(paste("schier_", SUBSET, ".csv", sep=""), sep=",", header=T, row.names=1, check.names = FALSE)
com <- t(com)
dim(com)
com[1:5, 1:2]

iteration = 999
cat("calculate RC-bray for", SUBSET)
df <- raup_crick_abundance(com, reps = iteration)
print("write csv")
write.csv(as.data.frame(as.matrix(df)), paste("RC-bray", SUBSET, iteration, ".csv", sep = ""))
str(df)




