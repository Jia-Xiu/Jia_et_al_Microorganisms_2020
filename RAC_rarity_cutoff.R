# using the rank abundance curve to calculate rarity cutoff per sample
# Built on 13-05-2019 by Jia, Xiu at the University of Groningen
# x.jia@rug.nl

library(vegan)

dataset.name = "KBr"

# load the rarefied otu table -----------------------------------------------------------------------
com <- read.csv("feature_table_rarefied.csv", sep=",", header=1, row.names=1, check.names = FALSE)
com <- t(com)
Chao <- as.data.frame(t(estimateR(com)))
Chao$slope <- Chao$S.obs/Chao$S.chao1

# using the rank abundance curve to calculate rarity cutoff per sample ------------------------------

# built empty matrix
cutoffs <- matrix(NA, nrow(com), 2)
row.names(cutoffs) <- row.names(com)
cutoffs[,2] <- Chao$slope
colnames(cutoffs) <- c("Rarity.cutoffs", "Slopes")

# find rarity cutoffs 
for (j in 1:nrow(com)) {
  com_j = sort(as.numeric(com[j,]), decreasing = TRUE)
  com_j <- com_j[com_j!=0]
  slope <- cutoffs[j,2]
  for (i in 1:length(com_j)){
    if (com_j[i]>=i*slope){
      H=i
    }
  }
  cutoffs[j, 1] <- H
}
cutoffs <- as.data.frame(cutoffs)
head(cutoffs)


# the rare biosphere
df <- com
for (j in 1:nrow(df)) {
  for (i in 1:ncol(df)) {
    if (df[j, i] > cutoffs[j,1]) {
      df[j, i] <- NA
    }
  }
}

write.csv(t(df), "rare.csv")

# the common biosphere
df <- com
for (j in 1:nrow(df)) {
  for (i in 1:ncol(df))
    if (df[j, i] <= cutoffs[j,1]) {
      df[j,i] <- NA
    }
}

write.csv(t(df), "common.csv")
