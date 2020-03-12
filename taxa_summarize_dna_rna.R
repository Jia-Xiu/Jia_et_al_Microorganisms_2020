# taxa summarize for the DNA/RNA datasets
# Author: Jia Xiu
# Date: 15-08-2019

rm(list=ls())

# Load the directory
directory = 'C:/Users/P278113/Dropbox'
# directory = '~/Dropbox/' 
subfolder = 'Schier/DNA_RNA'

setwd(paste(directory, subfolder, sep="/"))
getwd()

library(reshape2)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(doBy)
library(VennDiagram)
library(cowplot)
display.brewer.all()

marker = list(color = colorRampPalette(brewer.pal(11,"Spectral"))(100)) # balck back ground

mytheme <- theme_bw()+
  theme(text = element_text(size=12),
        legend.text=element_text(size=11),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1, 0.5, 1, 0.5),"cm")) 

# Check the overlap between the total and active subcommunities ------------------------------------------
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



# phyla changes along the successional gradients for the RNA/DNA dataset-----------------------------------------------------------------------------------
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


# loess fit plot 
phyla2$Year <- as.numeric(as.character(phyla2$Year))

(p2 <- ggplot(phyla2, aes(x = Year, y =value, color = Datasets, shape = Datasets)) + 
    geom_point(size = 2, alpha = 0.7)+
    geom_smooth(aes(fill = Datasets), method = loess, linetype="dashed") +
    facet_wrap(~Phyla, nrow=6, scales="free_y") +
    scale_color_manual(values = my_palette) + 
    scale_size_manual(values = c(16, 15)) +
    scale_x_continuous(breaks = c(0, 10, 40, 70, 110))+
    labs(x="Stage of succession (Years)", y="Richness")+
    mytheme + theme(legend.position = c(0.88, 0.06)))

ggsave("Phyla_loess_RNA_DNA.jpg", width = 17, height = 24, units = "cm", p2, scale = 1.1, device = "jpeg", dpi = 400)




# stacked-bar plot ----
colourCount <- length(levels(data$Phyla))
my_palette = c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"))
set.seed(1112)
my_palette <- sample(my_palette, colourCount)
pie(rep(1,colourCount), my_palette)
my_palette 

(f1 <- ggplot(data[data$Datasets == "cDNA", ], aes(x=Month, y=value.mean, fill=Phyla)) + 
    geom_bar(stat="identity", width=0.8) + #colour="black", 
    scale_y_continuous(expand = c(0, 0), limits = c(0,100))+
    facet_grid(~Year, switch = "x", scales = "free_x") +
    scale_fill_manual(values = my_palette)+
    guides(fill=guide_legend(reverse = TRUE, title = NULL)) +
    labs(x=" ",y="Relative abundance (%)", title="rRNA") +
    mytheme)

(f2 <- ggplot(data[data$Datasets == "DNA", ], aes(x=Month, y=value.mean, fill=Phyla)) + 
    geom_bar(stat="identity", width=0.8) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,100))+
    facet_grid(~Year, switch = "x", scales = "free_x") +
    scale_fill_manual(values = my_palette)+
    guides(fill=guide_legend(reverse = TRUE, title = NULL)) +
    labs(x=" ",y="Relative abundance (%)", title="rDNA") +
    mytheme)

f <- ggarrange(f1, f2, labels = c("A", "B"), common.legend = TRUE, legend = "bottom", ncol = 2)
f

ggsave("Phyla_RNA_DNA.png", width = 17, height = 11, units = "cm", f, scale = 1.5, dpi = 300)
ggsave("Phyla_RNA_DNA.pdf", width = 17, height = 11, units = "cm", f, scale = 1.5)






# calculate the rRNA:DNA ratio for each phylum ===========================================
com_dna[com_dna == 0] <- NA
df <- com_rna/com_dna
df[is.na(df)] <- 0
df <- df[, colSums(df != 0) > 0]
df <- t(df)
range(colSums(df)); dim(df)
str(df)
df[1:4, 1:2]

# inmport taxa info 
taxa$Phylum <- factor(gsub("Candidate division.", "", taxa$Phylum))
taxa$Phylum <- factor(gsub("Candidate division ", "", taxa$Phylum))
taxa$Order <- factor(gsub("uncultured bacterium", "", taxa$Order))
taxa$Family <- factor(gsub("uncultured", "", taxa$Family))
taxa$Genus <- factor(gsub("uncultured", "", taxa$Genus))
taxa$Species <- factor(gsub("uncultured bacterium", "", taxa$Species))
#taxa$info <- paste("p_", taxa$Phylum, "; c_", taxa$Class, "; o_", taxa$Order, "; f_", taxa$Family, "; g_", taxa$Genus, "; s_", taxa$Species)
head(taxa)
str(taxa)

taxa_Phylum <- as.data.frame(taxa$Phylum, row.names = row.names(taxa))
head(taxa_Phylum)

# merge the RNA/DNA ratio dataset with phylum info
df1 <- transform(merge(taxa_Phylum, df, by="row.names", all=TRUE), row.names=Row.names, Row.names=NULL) 
df1 <- na.omit(df1)
df1[1:4, 1:4]
levels(df1$taxa.Phylum)

df2 <- melt(df1, id=c("taxa.Phylum"))
df2 <- df2[!(df2$value == 0),]
df2$variable <- gsub("^X", "", df2$variable)
df2 <- df2[!(df2$taxa.Phylum == ""), ]
df2$taxa.Phylum <- factor(df2$taxa.Phylum, 
                          levels = c("Cyanobacteria", "JL-ETNP-Z39", "BHI80-139", "Fusobacteria", "KB1", 
                                     "Spirochaetae", "NPL-UPA2", "SM2F11", "Proteobacteria", "Dictyoglomi", 
                                     "Elusimicrobia", "Armatimonadetes", "Deinococcus-Thermus", "SHA-109", 
                                     "Tenericutes", "TA06", "BRC1", "Chlamydiae", "Chloroflexi", 
                                     "Deferribacteres", "Bacteroidetes", "Lentisphaerae", "Nitrospirae", 
                                     "Planctomycetes", "Firmicutes", 
                                     "Acidobacteria", "Actinobacteria", 
                                     "Fibrobacteres", "GAL08", "Gemmatimonadetes", "GOUTA4", 
                                     "OP11", "OP8", "RsaHF231", 
                                     "SR1", "TM7",  "WS3", "Verrucomicrobia", 
                                     "Chlorobi", "BD1-5", "TM6", "OD1", "OP3", "CKC4","WCHB1-60"))
str(df2)
head(df2)

(p <- ggplot(df2, aes(x = taxa.Phylum, y = log10(value)))+ 
    geom_jitter(shape=16, position=position_jitter(0.2), color = "darkgray", alpha = 0.6)+
    geom_violin(trim=FALSE, position=position_dodge(0.8), fill = "white", alpha = 0.4)+ 
    geom_boxplot(width=0.1, position=position_dodge(0.8))+
    geom_hline(yintercept = c(0), linetype = "dashed", colour = "#CC3300") + 
    #scale_fill_brewer(palette="Pastel1")+
    labs(x = "Phyla", y = "rRNA/rDNA (log10)", title = "")+
    theme_bw()+
    theme(text = element_text(size=14),
          panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text=element_text(size=11),
          legend.box.background = element_rect(),
          legend.box.margin = margin(1, 1, 1, 1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

ggsave("Relativeactivity_Phyla.jpg", width = 12, height = 5, units = "cm", p, scale = 2, device = "jpeg", dpi = 300)




