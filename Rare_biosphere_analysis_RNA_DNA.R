# Split the rare biosphere and common biosphere based on RNA/DNA data 
# Author: Jia Xiu
# Date: 2020-03-04

rm(list=ls())

# load packages
library(ggplot2)
library(VennDiagram)
library(ggpubr) 
library(reshape2)

# change directory
getwd()

# plot theme
mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face = "bold"),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# Set the cutoff for rarity
cutoff = 0.1/100

# load the rarefied otu table
com <- read.csv("feature_table_rarefied_taxon.csv", sep=",", header=1, row.names=1, check.names = FALSE)
com <- t(com[, 1:120])

# source the trucate function
source("TruncateTable.r") # see https://github.com/Jia-Xiu/rare_biosphere_assembly

# The truncated datasets can be stored as follows: 
truncated_ds_dominant <-TruncateTable(com, cutoff, typem="dominant") 
str(truncated_ds_dominant)
# write.csv(t(truncated_ds_dominant), paste("schier_dna_cdna_common", cutoff, "cutoff.csv", sep="_"))

truncated_ds_rare_without_dominant <-TruncateTable(com, cutoff, typem="rare") 
str(truncated_ds_rare_without_dominant)
# write.csv(t(truncated_ds_rare_without_dominant),  paste("schier_dna_cdna_rare", cutoff, "cutoff.csv", sep="_"))


# load the common and rare biospheres dataset
common <- read.csv(paste("schier_dna_cdna_common", cutoff, "cutoff.csv", sep="_"), sep=",", header=1, row.names=1)
common <- t(common)
common[is.na(common)] <- 0
common <- common[, colSums(common != 0) > 0]
common_dna <- subset(common, grepl('^DNA_', row.names(common)))
common_dna <- common_dna[, colSums(common_dna != 0) > 0]
common_cdna <- subset(common, grepl('^cDNA_', row.names(common)))
common_cdna <- common_cdna[, colSums(common_cdna != 0) > 0]
str(common_cdna)
cat("\nthe number of samples is:", nrow(common), "\nthe number of species/ASVs is:", ncol(common),
    "\nthe range of sequence number among samples is:", range(rowSums(common)))

rare <- read.csv(paste("schier_dna_cdna_rare", cutoff, "cutoff.csv", sep="_"), sep=",", header=1, row.names=1)
rare <- t(rare)
rare[is.na(rare)] <- 0
rare <- rare[, colSums(rare != 0) > 0]
rare_dna <- subset(rare, grepl('^DNA_', row.names(rare)))
rare_dna <- rare_dna[, colSums(rare_dna != 0) > 0]
rare_cdna <- subset(rare, grepl('^cDNA_', row.names(rare)))
rare_cdna <- rare_cdna[, colSums(rare_cdna != 0) > 0]
str(rare)
cat("\nthe number of samples is:", nrow(rare), "\nthe number of species/ASVs is:", ncol(rare),
    "\nthe range of sequence number among samples is:", range(rowSums(rare)))

# Venn digram for the overlap between rare and common 
venn.plot <- venn.diagram(
  x = list(
    common_RNA = colnames(common_cdna),
    common_DNA = colnames(common_dna),
    rare_RNA = colnames(rare_cdna),
    rare_DNA = colnames(rare_dna)
  ),
  filename = "Venn_4set_pretty.tiff",
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", 
                "white", "white", "white", "white", "darkblue", "white", 
                "white", "white", "white", "darkgreen", "white"),
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.07,
  cat.fontfamily = "serif",
  rotation.degree = 270,
  margin = 0.2
)




# calculate the rRNA:DNA ratio for the two biospheres 
# common biosphere
common <- read.csv(paste("schier_dna_cdna_common", cutoff, "cutoff.csv", sep="_"), sep=",", header=1, row.names=1)
common <- t(common)
common <- subset(common, grepl('^DNA_', row.names(common)))
row.names(common) <- gsub("^DNA_", "", row.names(common))
common[is.na(common)] <- 0
common <- common[, colSums(common != 0) > 0]

com_common_whole <- com_cdna[, colnames(com_cdna) %in% colnames(common)]
str(com_common_whole)
common[common == 0] <- NA
common <- com_common_whole/common
common[1:6, 1:2]
common[is.na(common)] <- 0
common <- common[, colSums(common != 0) > 0]
range(colSums(common)); dim(common); str(common)

group_info <- data.frame(row.names=rownames(common), t(as.data.frame(strsplit(rownames(common),"_"))))
head(group_info)

df <- data.frame(Year = as.factor(group_info[,1]),
                 Month = as.factor(group_info[,2]),
                 replicates = as.factor(group_info[,3]),
                 common)
df$Year <- factor(df$Year, levels = c("0", "10", "40", "70", "110")) 
df$Month <- factor(df$Month, levels = c("5", "7", "9", "11"), 
                   labels=c("May", "Jul", "Sep", "Nov"))
df <- melt(df, id.vars = c("Year", "Month", "replicates"))
df <-df[!(df$value == 0),]
str(df)

(p1 <- ggplot(df, aes(x = Year, y = log10(value), fill = Month))+ 
    geom_violin(trim=FALSE, position=position_dodge(0.8))+ 
    geom_boxplot(width=0.1, position=position_dodge(0.8))+
    geom_hline(yintercept = c(0), linetype = "dashed", colour = "#CC3300") + 
    scale_fill_brewer(palette="Pastel1")+
    labs(x = "Successional stages (year)", y = "rRNA/rDNA (log10)", title = "The common biosphere")+
    mytheme)

(f2 <- ggplot(data = df, aes(x=value)) + 
    geom_density(aes(fill=Month), alpha = 0.5)+ 
    facet_wrap( ~ Year, nrow = 1)+
    geom_vline(xintercept = c(1), linetype = "dashed", colour = "#CC3300") + 
    scale_x_continuous(limits = c(0, 30))+
    scale_fill_brewer(palette="Pastel1")+
    labs(x = "rRNA/rDNA", y = "Density", title = "The common biosphere")+
    mytheme)


# rare biosphere
rare <- read.csv(paste("schier_dna_cdna_rare", cutoff, "cutoff.csv", sep="_"), sep=",", header=1, row.names=1)
rare <- t(rare)
rare <- subset(rare, grepl('^DNA_', row.names(rare)))
row.names(rare) <- gsub("^DNA_", "", row.names(rare))
rare[is.na(rare)] <- 0
rare <- rare[, colSums(rare != 0) > 0]

com_rare_whole <- com_cdna[, colnames(com_cdna) %in% colnames(rare)]
str(com_rare_whole)
rare[rare == 0] <- NA # OR 1
rare <- com_rare_whole/rare
rare[1:6, 1:2]
rare[is.na(rare)] <- 0
rare <- rare[, colSums(rare != 0) > 0]
range(colSums(rare)); dim(rare); str(rare)

group_info <- data.frame(row.names=rownames(rare), t(as.data.frame(strsplit(rownames(rare),"_"))))
head(group_info)

df <- data.frame(Year = as.factor(group_info[,1]),
                 Month = as.factor(group_info[,2]),
                 replicates = as.factor(group_info[,3]),
                 rare)
df$Year <- factor(df$Year, levels = c("0", "10", "40", "70", "110")) #labels = c("0yr", "10yr", "40yr", "70yr", "110yr")
df$Month <- factor(df$Month, levels = c("5", "7", "9", "11"), labels=c("May", "Jul", "Sep", "Nov"))
df <- melt(df, id.vars = c("Year", "Month", "replicates"))
df <- df[!(df$value == 0),]
head(df)

(p1 <- ggplot(df, aes(x = Year, y = log10(value), fill = Month))+ 
    geom_violin(trim=FALSE, position=position_dodge(0.8))+ 
    geom_boxplot(width=0.1, position=position_dodge(0.8))+
    #geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+
    geom_hline(yintercept = c(0), linetype = "dashed", colour = "#CC3300") + 
    scale_fill_brewer(palette="Pastel1")+
    labs(x = "Successional stages (year)", y = "rRNA/rDNA (log10)", title = "The rare biosphere")+
    mytheme)

(f1 <- ggplot(data = df, aes(x=value)) + geom_density(aes(fill=Month), alpha = 0.3)+ 
    facet_wrap( ~ Year, nrow = 1)+
    geom_vline(xintercept = c(1), linetype = "dashed", colour = "#CC3300") + 
    scale_x_continuous(limits = c(0, 30))+
    scale_fill_brewer(palette="Pastel1")+
    labs(x = "rRNA/rDNA", y = "Density", title = "The rare biosphere")+
    mytheme)


p <- ggarrange(p1, p2, labels = c("A", "B"), 
               align = "v", font.label = list(size = 14, face = "bold"), 
               common.legend = TRUE, legend = "right", ncol = 1, nrow = 2)
p

ggsave("Relative_activity_common_rare_log.png", width = 17, height = 10, units = "cm", p, scale = 2, dpi = 300)
ggsave("Relative_activity_common_rare_log.pdf", width = 17, height = 10, units = "cm", p, scale = 2)

f <- ggarrange(f1, f2, labels = c("A", "B"), 
               align = "v", font.label = list(size = 14, face = "bold"), 
               common.legend = TRUE, legend = "right", ncol = 1, nrow = 2)
f

ggsave("Relative_activity_common_rare_density.png", width = 17, height = 10, units = "cm", f, scale = 2, dpi = 300)
ggsave("Relative_activity_common_rare_density.pdf", width = 17, height = 10, units = "cm", f, scale = 2)


