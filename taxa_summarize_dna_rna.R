# taxa summarize for the DNA/RNA datasets
# Author: Jia Xiu
# Date: 15-08-2019

rm(list=ls())

# Load the directory
setwd()

library(reshape2)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(doBy)
library(VennDiagram)
library(cowplot)
display.brewer.all()


# read whole dataset
com <- read.csv("feature_table_rarefied_taxon.csv", sep=",", header=1, row.names=1, check.names = FALSE)

# split community info and taxa info
taxa <- com[, c(122:127)]
com <- t(com[, 1:120])
cat("rarefaction depth", rarefactiondepth <- mean(rowSums(com[1:120,])))

com_rna <- subset(com, grepl('^cDNA_', row.names(com)))
row.names(com_rna) <- gsub("^cDNA_", "", row.names(com_rna))
com_dna <- subset(com, grepl('^DNA_', row.names(com)))
row.names(com_dna) <- gsub("^DNA_", "", row.names(com_dna))
str(com_dna)
com_dna[1:15, 1:2]



# phyla changes along the successional gradients for the RNA/DNA dataset
com_taxa <- t(com)
com_taxa[is.na(com_taxa)] <- 0
com_taxa <- com_taxa[rowSums(com_taxa) > 0, ]

# add taxa info
head(taxa)

com_taxa <- transform(merge(taxa, com_taxa, by="row.names", all=FALSE), row.names=Row.names, Row.names=NULL)
com_taxa <- com_taxa[, -c(2:6)]
com_taxa[1:5, 1:5]

phyla <- aggregate(com_taxa[,-1], list(com_taxa$Phylum), sum)
row.names(phyla) <- phyla$Group.1
phyla <- phyla[, -1]
phyla <- 100*phyla/rarefactiondepth
phyla$mean <- rowMeans(phyla)
phyla <- subset(phyla, mean > 0.1) 
#phyla <- phyla[with(phyla, order(phyla$cDNA_110_11_A)), ]
phyla <- phyla[with(phyla, order(-mean)), ]
phyla$mean <- NULL
phyla <- t(phyla)
phyla <- as.data.frame(phyla)
colnames(phyla)
phyla[1:5, 1:5]

# split treatment info
group_info<-data.frame(row.names=rownames(phyla), t(as.data.frame(strsplit(rownames(phyla),"_"))))
head(group_info)

# combine treatment info with phyla relative abundance
phyla2 <- data.frame(Datasets = as.factor(group_info[,1]),
                     Year = as.factor(group_info[,2]),
                     Month = as.factor(group_info[,3]),
                     replicates = as.factor(group_info[,4]),
                     phyla)
row.names(phyla2) <- row.names(phyla)

phyla2 <- melt(phyla2, id=c("Datasets", "Year", "Month", "replicates"), variable = 'Phyla')

phyla2$Datasets <- factor(phyla2$Datasets, levels = c("cDNA", "DNA"), labels = c("RNA", "DNA"))

phyla2$Year <- factor(phyla2$Year, levels=c("0", "10", "40", "70", "110"))

phyla2$Month <- factor(phyla2$Month, levels=c("5", "7", "9", "11"), 
                       labels=c("M", "J", "S", "N"))#labels=c("May", "Jul", "Sep", "Nov"))

names(table(phyla2$Phyla))
str(phyla2)


# Line plot
pd <- position_dodge(0.2)  
(my_palette = c(brewer.pal(9, "Set1")[c(1,2)]))


dstats<-function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data <- summaryBy(value ~ Datasets + Year + Phyla, data=phyla2, FUN=dstats)
head(data)

mytheme <- theme_bw()+
  theme(text = element_text(size = 12),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold", size = 11),
        strip.background = element_blank(),
        strip.text=element_text(face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

(p1 <- ggplot(data, aes(x = Year, y = value.mean, group = Datasets, colour = Datasets)) + 
    geom_errorbar(aes(ymin=value.mean-value.se, ymax=value.mean+value.se), colour="gray", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3) +
    scale_color_brewer(palette="Set1") +
    facet_wrap(~Phyla, nrow=6, scales="free_y") +
    labs(title="", x = "Stage of succession (Year)", y = "Relative abundance (%)")+
    mytheme + theme(legend.position = c(0.88, 0.06)))

#ggsave("Phyla_RNA_DNA.jpg", width = 17, height = 24, units = "cm", p1, scale = 1, device = "jpeg",dpi = 400)
