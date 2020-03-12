# beta diversity analysis for the DNA/RNA dataset
# Update: 11-03-2020
# Author: Jia, Xiu

rm(list=ls())

# Load libraries
library(vegan) # for multivariable analysis
library(ape) # for pcoa 
library(ggplot2) # for graphing
library(ggpubr) # combine figures
library(RColorBrewer) # for color bar
library(plyr) # for rename
display.brewer.all()

# Load the data into the R workspace and run the scripts: community structure
directory = 'C:/Users/P278113/Dropbox'
# directory = '~/Dropbox/' # mac
subfolder = 'Schier/DNA_RNA'

setwd(paste(directory, subfolder, sep="/"))
getwd()

# ggplot theme
mytheme <- theme_bw()+
  theme(text = element_text(size=12),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  

# function for do poca analysis
pcoa.unfrac <- function(dist) {
  outlist <- list()
  # change row and column names
  row.names(dist) <- gsub("\\-", "\\_", row.names(dist))
  colnames(dist) <- gsub("\\-", "\\_", colnames(dist))
  dist[1:5, 1:5]
  # change as distance matirx
  dist <- as.dist(as.matrix(dist))
  str(dist)
  outlist[["dist"]] <- dist
  
  # PCoA 
  re <- pcoa(dist, correction="none", rn=NULL)
  str(re)
  
  # data.frame for ploting
  group_info <- data.frame(row.names=rownames(re$vectors), t(as.data.frame(strsplit(rownames(re$vectors),"_"))))
  head(group_info)
  
  df <- data.frame(pc1 = re$vectors[,1], 
                   pc2 = re$vectors[,2],
                   pc3 = re$vectors[,3],
                   Datasets=as.factor(group_info[,1]),
                   Year=as.factor(group_info[,2]),
                   Month=as.factor(group_info[,3]),
                   replicates=as.factor(group_info[,4]))
  
  df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))
  df$Datasets <- factor(df$Datasets, levels = c("cDNA", "DNA"), labels = c("RNA", "DNA"))
  
  outlist[["df"]] <- df
  
  re_eig <- c(round(re$values$Relative_eig[1] * 100, 2), 
              round(re$values$Relative_eig[2] * 100, 2), 
              round(re$values$Relative_eig[3] * 100, 2))
  
  outlist[["re.eig"]] <- re_eig
  
  return(outlist)
}

# read Unifrac distances
dist <- read.table('weighted_unifrac_distance_matrix.tsv', sep = '\t', header=1, row.names=1, check.names = FALSE)

data <- pcoa.unfrac(dist)

(p1 <- ggplot(data$df, aes(pc1, pc2, shape=Datasets, fill=Year))+
    geom_point(size=3, alpha=0.7)+ 
    labs(x=paste("PCoA1 (", data$re.eig[1], "%)", sep=""), 
         y=paste("PCoA2 (", data$re.eig[2], "%)", sep=""), 
         title="Weighted Unifrac")+
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme)

(p2 <- ggplot(data$df, aes(pc1, pc3, shape=Datasets, fill=Year))+
    geom_point(size=3, alpha=0.7)+ 
    labs(x=paste("PCoA1 (", data$re.eig[1], "%)", sep=""), 
         y=paste("PCoA3 (", data$re.eig[3], "%)", sep=""), 
         title="Weighted Unifrac")+
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme)

(p3 <- ggplot(data$df, aes(pc2, pc3, shape=Datasets, fill=Year))+
    geom_point(size=3, alpha=0.7)+ 
    labs(x=paste("PCoA2 (", data$re.eig[2], "%)", sep=""), 
         y=paste("PCoA3 (", data$re.eig[3], "%)", sep=""), 
         title="Weighted Unifrac")+
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme)

# Unweighted Unifrac
dist <- read.table('unweighted_unifrac_distance_matrix.tsv', sep = '\t', header=1, row.names=1, check.names = FALSE)

data <- pcoa.unfrac(dist)
str(data)

(p4 <- ggplot(data$df, aes(pc1, pc2, shape=Datasets, fill=Year))+
    geom_point(size=3, alpha=0.7)+ 
    labs(x=paste("PCoA1 (", data$re.eig[1], "%)", sep=""), 
         y=paste("PCoA2 (", data$re.eig[2], "%)", sep=""), 
         title="Unweighted Unifrac")+
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme)

(p5 <- ggplot(data$df, aes(pc1, pc3, shape=Datasets, fill=Year))+
    geom_point(size=3, alpha=0.7)+ 
    labs(x=paste("PCoA1 (", data$re.eig[1], "%)", sep=""), 
         y=paste("PCoA3 (", data$re.eig[3], "%)", sep=""), 
         title="Unweighted Unifrac")+
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme)

(p6 <- ggplot(data$df, aes(pc2, pc3, shape=Datasets, fill=Year))+
    geom_point(size=3, alpha=0.7)+ 
    labs(x=paste("PCoA2 (", data$re.eig[2], "%)", sep=""), 
         y=paste("PCoA3 (", data$re.eig[3], "%)", sep=""), 
         title="Unweighted Unifrac")+
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme)

(p <- annotate_figure((ggarrange(p1, p4, p2, p5, p3, p6, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
  common.legend = TRUE, legend = "right", ncol = 2, nrow = 3, align = "hv")), right = " "))

#ggsave("PCoA_Unifrac_Datasets_Year.pdf", width = 12, height = 15, units = "cm", p, scale = 1.5)
#ggsave("PCoA_Unifrac_Datasets_Year.jpg", width = 12, height = 15, units = "cm", device = "jpeg", p, scale = 1.5, dpi = 300)


# Bray-Curtis
com <- read.csv("feature_table_rarefied_taxon.csv", sep=",", header=1, row.names=1, check.names = FALSE)
com <- t(com[, 1:120])

dist <- vegdist(com, method="bray", binary=FALSE, diag=1) 

re <- pcoa(dist, correction="none", rn=NULL)

group_info <- data.frame(row.names=rownames(re$vectors), t(as.data.frame(strsplit(rownames(re$vectors),"_"))))

df <- data.frame(pc1 = re$vectors[,1], 
                 pc2 = re$vectors[,2],
                 pc3 = re$vectors[,3],
                 Datasets=as.factor(group_info[,1]),
                 Year=as.factor(group_info[,2]),
                 Month=as.factor(group_info[,3]),
                 replicates=as.factor(group_info[,4]))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))
df$Datasets <- factor(df$Datasets, levels = c("cDNA", "DNA"), labels = c("RNA", "DNA"))


(p1 <- ggplot(df, aes(pc1, pc2, shape=Datasets, fill=Year))+
    geom_point(size=3, alpha=0.7)+ 
    labs(x=paste("PCoA1 (", round(re$values$Relative_eig[1] * 100, 2), "%)", sep=""), 
         y=paste("PCoA2 (", round(re$values$Relative_eig[2] * 100, 2), "%)", sep=""), 
         title="")+
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme)

(p2 <- ggplot(df, aes(pc1, pc3, shape=Datasets, fill=Year))+
    geom_point(size=3, alpha=0.7)+ 
    labs(x=paste("PCoA1 (", round(re$values$Relative_eig[1] * 100, 2), "%)", sep=""), 
         y=paste("PCoA3 (", round(re$values$Relative_eig[3] * 100, 2), "%)", sep=""), 
         title="")+
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme) 

(p3 <- ggplot(df, aes(pc2, pc3, shape=Datasets, fill=Year))+
    geom_point(size=3, alpha=0.7)+ 
    labs(x=paste("PCoA2 (", round(re$values$Relative_eig[2] * 100, 2), "%)", sep=""), 
         y=paste("PCoA3 (", round(re$values$Relative_eig[3] * 100, 2), "%)", sep=""), 
         title="")+
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme) 

(p <- annotate_figure((ggarrange(p1, p2, p3, labels = c("(a)", "(b)", "(c)"), 
                common.legend = TRUE, legend = "right", ncol = 3, nrow = 1)), right = " "))
  
#ggsave("PCoA_bray_Datasets_Year.jpg", width = 17, height = 5, units = "cm", device = "jpeg", p, scale = 1.5, dpi = 300)

######################################################################################################
# Turnover of DNA- and RNA-inferred community across successional stages and sampling month
######################################################################################################


# pairwise similarity
dist <- data$dist
df <- as.matrix(dist)
df <- data.frame(as.table(df))[lower.tri(df, diag = FALSE), ]
cat("For no. of df pairs: Observed = Predict", (120*120-120)/2 == length(df$Freq), "\n")
row.names(df) <- paste(df$Var1, df$Var2, sep = "_")

# split treatment info
group <- data.frame(row.names=rownames(df), t(as.data.frame(strsplit(as.character(row.names(df)), "_"))))
head(group)

df <- transform(merge(df, group[, -c(4, 8)], by="row.names"), row.names=Row.names, Row.names=NULL, check.names=FALSE)
df <- df[, -c(1:2)]
df <- df[which(df$X1==df$X5),]
head(df)

# Year
dfy <- df[which(df$X2==df$X6),]
dfy <- dfy[which(dfy$X3!=dfy$X7),]
dfy <- dfy[, c(1:3)]
dfy$Year <- factor(dfy$X2, levels=c("0", "10", "40", "70", "110"))
dfy$Datasets <- factor(dfy$X1, levels=c("cDNA", "DNA"), labels=c("RNA", "DNA"))
str(dfy)
head(dfy)

# Wilcoxin - Year
years <- c("0", "10", "40", "70", "110")
for (yr in years) {
  cat("for Year:", yr)
  res <- wilcox.test(Freq ~ Datasets, data = dfy[which(dfy$Year == yr),], exact = FALSE)
  print(res)
}

pd <- position_dodge(0.2)
(my_palette = c(brewer.pal(9, "Set1")[c(1,2)]))

(p1 <- ggplot(dfy, aes(x = Year, y = Freq, fill = Datasets))+ 
    geom_point(size = 2, shape = 16, alpha = 0.3, position = position_jitterdodge()) +
    geom_boxplot(alpha = 0.8, outlier.size=-1) +
    scale_fill_manual(values = my_palette) + 
    scale_y_continuous(limits=c(0.04, 0.5)) + 
  labs(x="Stage of succession (Years)", y="Temporal turnover (weighted UniFrac)", title="") + 
    stat_compare_means(aes(group = supp), label = "p.signif",  label.y = c(16, 25, 29)) + ???
    mytheme)

#ggsave("Tunrover_unifrac_Datasets_x_year_raw.pdf", width = 6, height = 5, p1, scale = 1)


# Month
dfm <- df[which(df$X2!=df$X6),]
dfm <- dfm[which(dfm$X3==dfm$X7),]
dfm <- dfm[, c(1:2, 4)]
dfm$Month <- factor(dfm$X3, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "Jul", "Sep", "Nov"))
dfm$Datasets <- factor(dfm$X1, levels=c("cDNA", "DNA"), labels=c("RNA", "DNA"))
str(dfm)
head(dfm)

# Wilcoxin - month
months <- c(c("May", "Jul", "Sep", "Nov"))
for (mon in months) {
  cat("for Month:", mon)
  res <- wilcox.test(Freq ~ Datasets, data = dfm[which(dfm$Month == mon),], exact = FALSE)
  print(res)
}

(p2 <- ggplot(dfm, aes(x = Month, y = Freq, fill = Datasets))+ 
    geom_point(size = 2, shape = 16, alpha = 0.3, position = position_jitterdodge()) +
    geom_boxplot(alpha = 0.8, outlier.size=-1) +
    scale_fill_manual(values = my_palette) + 
    scale_y_continuous(limits=c(0.04, 0.6)) + 
    labs(x="Sampling month", y="Spatial turnover (weighted UniFrac)", title="") + 
    mytheme)

#ggsave("Tunrover_unifrac_Datasets_x_month_raw.pdf", width = 5.8, height = 5, p2, scale = 1)

(p <- annotate_figure((ggarrange(p1, p2, labels = c("(a)", "(b)", widths = c(1, 0.8)), 
                                 common.legend = TRUE, legend = "right", ncol = 2, nrow = 1)), right = " "))

#ggsave("Tunrover_unifrac_Datasets_raw.pdf", width = 15, height = 6, units = "cm", p, scale = 1.5)



# two-way anova of turnover ------------------------------------------

#kruskal-wallis
kt <- kruskal.test(Freq ~ Datasets, data = dfy[which(dfy$Year == "0"),])
kt
#summary(aov)

# check the homogeneity of variances
bartlett.test(Freq ~ Datasets, data=dfy)
bartlett.test(Freq ~ Year, data=dfy)
bartlett.test(Freq ~ interaction(Datasets, Year), data=dfy)
# p > 0.05, suggesting that the variance across groups is not statistically significantly different.

# check normality
ggqqplot(dfy$Freq)

# Shapiro-Wilk test of normality for univariate
shapiro.test(dfy$Freq)
# p-value > 0.05, implying that the distribution of the data are not significantly different from normal distribution. 



# NMDS ----------------------------------------------------------------------------------------------------
# ref. ?
# First step is to calculate a distance matrix, i.e. bray-curtis distance, which is recommended for abundance data
dist <- vegdist(com,  method = "bray")

# In this part, we define a function NMDS.scree() that automatically 
# performs a NMDS for 1-10 dimensions and plots the nr of dimensions vs the stress
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# *A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions,
# < 0.1 is great, < 0.2 is good/ok, and stress < 0.3 provides a poor representation.

# Use the function that we just defined to choose the optimal nr of dimensions
NMDS.scree(dist)

# Because the final result depends on the initial random placement of the points 
# we`ll set a seed to make the results reproducible
set.seed(2)

# Here, we perform the final analysis and check the result
# If you don`t provide a dissimilarity matrix, metaMDS automatically applies Bray-Curtis. 
NMDS1 <- metaMDS(dist, k = 2, trymax = 100, trace = F, autotransform = F)
# Do you know what the trymax = 100 and trace = F means?
# Let's check the results
NMDS1

# stressplot draws a Shepard diagram which is a plot of ordination distances and monotone or linear fit line against 
# original dissimilarities. In addition, it displays two correlation-like statistics on the goodness of fit in the graph. 
# The nonmetric fit is based on stress S and defined as R2 = 1-S*S. 
# The ???linear fit??? is the squared correlation between fitted values and ordination distances. 
# For monoMDS, the ???linear fit??? and R2 from ???stress type 2??? are equal.
stressplot(NMDS1)

# There is a good non-metric fit between observed dissimilarities (in our distance matrix) and 
# the distances in ordination space. Also the stress of our final result was ok 
# (do you know how much the stress is?). So we can go further and plot the results:
plot(NMDS1, type = "t")

# There are no species scores (same problem as we encountered with PCoA). 
# We can work around this problem, by giving metaMDS the original community matrix as input and specifying the distance measure.

NMDS <- metaMDS(com, k = 2, trymax = 100, trace = F, autotransform = F, distance="bray")
plot(NMDS, display = "sites", type = "n")
points(NMDS, display = "sites", col = "red", cex = 1.25)

str(NMDS)
head(NMDS$points)
group_info <- data.frame(row.names=row.names(NMDS$points), 
                         t(as.data.frame(strsplit(as.character(row.names(NMDS$points)), "_"))))
head(group_info)

#Make a new data frame, with Sites and Groups information there, to be useful for coloring, and shape of points
df <- data.frame(x=NMDS$points[,1],
                 y=NMDS$points[,2],
                 Fraction=as.factor(group_info[,1]),
                 Year=as.factor(group_info[,2]),
                 Month=as.factor(group_info[,3]),
                 replicates=as.factor(group_info[,4]))

df$Dataset <- factor(df$Dataset, levels = c("cDNA", "DNA"), labels = c("rRNA subset", "rDNA subset"))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
str(df); head(df)


(f <- ggplot(df, aes(x, y, shape=Fraction, fill=Year))+
    geom_point(size=6, alpha=0.7)+ 
    labs(x= "NMDS1", y = "NMDS2", title = " ", shape="Biosphere")+
    annotate('text', label = paste("Stress =", round(NMDS3$stress, 2)), x = Inf, y = Inf, size = 4, vjust =1.6, hjust = 1.2) +
    scale_fill_brewer(palette="Accent", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(21, 22))+ 
    mytheme
)

#ggsave("NMDS_bray_Fraction_x_year..pdf", width = 8.5, height = 6, units = "cm", f, scale = 2)
#ggsave("NMDS_bray_Fraction_x_year.png", width = 8.5, height = 6, units = "cm", f, scale = 2, dpi = 300)

# ==================================================================================================================
# PERMANOVA based on Unifrac 
# ==================================================================================================================

# whole dataset
dist <- read.table('unweighted_unifrac_distance_matrix.tsv', sep = '\t', header=1, row.names=1, check.names = FALSE)
row.names(dist) <- gsub("\\-", "\\_", row.names(dist))
colnames(dist) <- gsub("\\-", "\\_", colnames(dist))
dist[1:5, 1:5]

df <- data.frame(row.names=rownames(dist), t(as.data.frame(strsplit(rownames(dist),"_"))))
df$Dataset <- factor(df$X1, levels = c("cDNA", "DNA"), labels = c("RNA", "DNA"))
df$Year <- factor(df$X2, levels=c("0", "10", "40", "70", "110"))
df$Month <- factor(df$X3, levels=c("5", "7", "9", "11"), labels=c("May", "July", "September", "November"))
head(df)


# change as distance matirx
dist <- as.dist(as.matrix(dist))
str(dist)

# two way permanova (successional stages vs. season)
set.seed(123)
(result_whole <- adonis(dist ~ Dataset*Year*Month, data=df, method="bray", permutation=9999))
result_whole
#write.csv(result_whole$aov.tab, "permanova_unweightedUnifrac_dataset_x_year_x_month.csv")


# calculate the impact on each successional year ---------------------------------------------------------
permanova.list <- list()
dispersion.list <- list()
dispersion.sites.list <- list()

years <- c("0", "10", "40", "70", "110")

dist <- read.table('unweighted_unifrac_distance_matrix.tsv', sep = '\t', header=1, row.names=1, check.names = FALSE)
row.names(dist) <- gsub("\\-", "\\_", row.names(dist))
colnames(dist) <- gsub("\\-", "\\_", colnames(dist))
dist[1:5, 1:5]

for (yr in years) {
  dist_sub <- subset(dist, grepl(paste("_", yr, "_", sep = ""), row.names(dist)))
  dist_sub <- t(dist_sub)
  dist_sub <- subset(dist_sub, grepl(paste("_", yr, "_", sep = ""), row.names(dist_sub)))
  try(if(nrow(dist_sub) !=24) stop("Error in subset com!"))
  df <- data.frame(row.names=rownames(dist_sub), t(as.data.frame(strsplit(rownames(dist_sub),"_"))))
  df$Dataset <- factor(df$X1, levels = c("cDNA", "DNA"), labels = c("RNA", "DNA"))
  df$Year <- factor(df$X2, levels=c("0", "10", "40", "70", "110"))
  df$Month <- factor(df$X3, levels=c("5", "7", "9", "11"), labels=c("May", "July", "September", "November"))
  head(df)
  
  cat("\nTwo-way permanova (Dataset vs. Month) for year:", yr)
  set.seed(123)
  result_whole <- adonis(dist_sub ~ Dataset*Month, data=df, permutation=9999) # two way permanova
  result_whole$aov.tab
  permanova.list[[yr]] <- as.data.frame(result_whole$aov.tab)
  
  cat("\nDispersion between DNA & RNA for year:", yr) 
  dist_sub <- as.dist(as.matrix(dist_sub))
  dispersion <- betadisper(dist_sub, group = df$Dataset)
  dispersion.list[[yr]] <- permutest(dispersion, permutation=9999)$tab
  p <- plot(dispersion, hull=FALSE, ellipse=TRUE, main = paste("Dispersion (", yr, " yr)", sep = "")) ##sd ellipse
  dispersion.sites.list[[yr]] <- p$sites
}

str(permanova.list)
df <- do.call(rbind.data.frame, permanova.list)
df$Year <- gsub("\\..*", "", row.names(df))
df$ii <- sub('.*\\.', '', row.names(df))
#write.csv(df, "permanova_weighted_unifrac_dataset_years.csv")

str(dispersion.list)
df <- do.call(rbind.data.frame, dispersion.list )
df$Year <- gsub("\\..*", "", row.names(df))
df$ii <- sub('.*\\.', '', row.names(df))
write.csv(df, "dispersion_permanova_weighted_unifrac_dataset_years.csv")

str(dispersion.sites.list )
dispersion.sites <- do.call(rbind.data.frame, dispersion.sites.list )
dispersion.sites$Year <- gsub("\\..*", "", row.names(dispersion.sites))
dispersion.sites$sample <- sub('.*\\.', '', row.names(dispersion.sites))
#write.csv(df, "dispersion_weighted_unifrac_sites_dataset_years.csv")


# RNA and DNA subset -----------------------------------------------------------------------------------------------
dist <- read.table('unweighted_unifrac_distance_matrix.tsv', sep = '\t', header=1, row.names=1, check.names = FALSE)
row.names(dist) <- gsub("\\-", "\\_", row.names(dist))
colnames(dist) <- gsub("\\-", "\\_", colnames(dist))
dist[1:5, 1:5]

# RNA subset 
dist_RNA <- subset(dist, grepl('^cDNA_', row.names(dist)))
dist_RNA <- dist_RNA[, grepl('^cDNA_', names(dist_RNA))]
dim(dist_RNA)
dist_RNA[1:5, 1:2]

df <- data.frame(row.names=rownames(dist_RNA), t(as.data.frame(strsplit(rownames(dist_RNA),"_"))))
df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))
df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))
df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), labels=c("May", "July", "September", "November"))
head(df)

dist_RNA <- as.dist(as.matrix(dist_RNA))
str(dist_RNA)

# two way permanova (successional stages vs. season)
set.seed(123)
result_RNA <- adonis(dist_RNA ~ Month*Year, data=df, method="bray", permutation=9999) # two way permanova
result_RNA

# DNA subset
dist_DNA <- subset(dist, grepl('^DNA_', row.names(dist)))
dist_DNA <- dist_DNA[, grepl('^DNA_', names(dist_DNA))]
dim(dist_DNA)
dist_DNA[1:5, 1:2]

df <- data.frame(row.names=rownames(dist_DNA), t(as.data.frame(strsplit(rownames(dist_DNA),"_"))))
df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))
df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))
df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), labels=c("May", "July", "September", "November"))
head(df)

dist_DNA <- as.dist(as.matrix(dist_DNA))
str(dist_DNA)

# two way permanova (successional stages vs. season)
set.seed(123)
result_DNA <- adonis(dist_DNA ~ Month*Year, data=df, method="bray", permutation=9999) # two way permanova
result_DNA

(results <- rbind(result_RNA$aov.tab, result_DNA$aov.tab))
#write.csv(results, "permanova_unifrac_DNA_RNA_subsets.csv")

######################################################################################################
# dispersion analysis across successional years
######################################################################################################


# calculate the impact of dataset of each year ----------
permanova.list <- list()
dispersion.list <- list()
dispersion.sites.list <- list()

years <- c("0", "10", "40", "70", "110")
yr <- "70"
for (yr in years) {
  com_sub <- subset(com, grepl(paste("_", yr, "_", sep = ""), row.names(com)))
  try(if(nrow(com_sub) !=24) stop("Error in subset com!"))
  df <- data.frame(row.names=rownames(com_sub), t(as.data.frame(strsplit(rownames(com_sub),"_"))))
  df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))
  df
  
  cat("\nTwo-way permanova (Dataset vs. Month) for year:", yr)
  set.seed(123)
  result_whole <- adonis(com_sub ~ Dataset*Month, data=df, method="bray", permutation=9999) # two way permanova
  result_whole$aov.tab
  permanova.list[[yr]] <- as.data.frame(result_whole$aov.tab)
  
  cat("\nDispersion between DNA & RNA for year:", yr) 
  dist <- vegdist(com_sub, method="bray", binary=FALSE, diag=1) 
  dispersion <- betadisper(dist, group = df$Dataset)
  dispersion.list[[yr]] <- permutest(dispersion, permutation=9999)$tab
  p <- plot(dispersion, hull=FALSE, ellipse=TRUE, main = paste("Dispersion (", yr, " yr)", sep = "")) ##sd ellipse
  dispersion.sites.list[[yr]] <- p$sites
}

str(permanova.list)
df <- do.call(rbind.data.frame, permanova.list)
df$Year <- gsub("\\..*", "", row.names(df))
df$ii <- sub('.*\\.', '', row.names(df))
write.csv(df, "permanova_dataset_years.csv")

str(dispersion.list)
df <- do.call(rbind.data.frame, dispersion.list )
df$Year <- gsub("\\..*", "", row.names(df))
df$ii <- sub('.*\\.', '', row.names(df))
write.csv(df, "dispersion_permanova_dataset_years.csv")

str(dispersion.sites.list )
dispersion.sites <- do.call(rbind.data.frame, dispersion.sites.list )
dispersion.sites$Year <- gsub("\\..*", "", row.names(dispersion.sites))
dispersion.sites$sample <- sub('.*\\.', '', row.names(dispersion.sites))
write.csv(df, "dispersion_sites_dataset_years.csv")

# all communities
df <- data.frame(row.names=dispersion.sites$sample, t(as.data.frame(strsplit(dispersion.sites$sample,"_"))))

df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))

df$Dataset <- factor(df$Dataset, levels = c("cDNA", "DNA"), labels = c("rRNA", "rDNA"))

df <- data.frame(df, dispersion.sites[, 1:2])

head(df)

my_palette <- c(brewer.pal(9, "Set1")[c(1:2)])

(f1 <- ggplot(df, aes(PCoA1, PCoA2, shape=Month, fill=Dataset))+
    facet_grid(~Year, scales = "free_x") +
    geom_point(size=4, alpha=0.7)+ 
    labs(x= "PCoA1", y = "PCoA2", title = " ")+
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(24, 22, 23, 21))+
    mytheme
)

# RNA subset ----------------------------------------------
com_RNA <- subset(com, grepl('^cDNA_', row.names(com)))
com_RNA[is.na(com_RNA)] <- 0
com_RNA <- com_RNA[, colSums(com_RNA != 0) > 0]
dim(com_RNA)
com_RNA[1:5, 1:2]
#write.csv(t(com_RNA), "schier_RNA.csv")

df <- data.frame(row.names=rownames(com_RNA), t(as.data.frame(strsplit(rownames(com_RNA),"_"))))

df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
head(df)

# two way permanova (successional stages vs. season)
set.seed(123)
result_RNA <- adonis(com_RNA ~ Month*Year, data=df, method="bray", permutation=9999) # two way permanova
result_RNA

# DNA subset
com_DNA <- subset(com, grepl('^DNA_', row.names(com)))
com_DNA[is.na(com_DNA)] <- 0
com_DNA <- com_DNA[, colSums(com_DNA != 0) > 0]
dim(com_DNA)
com_DNA[1:5, 1:2]
#write.csv(t(com_DNA), "schier_DNA.csv")

df <- data.frame(row.names=rownames(com_DNA), t(as.data.frame(strsplit(rownames(com_DNA),"_"))))

df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))
df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
head(df)

set.seed(123)
result_DNA <- adonis(com_DNA ~ Month*Year, data=df, method="bray", permutation=9999) # two way permanova
result_DNA

(results <- rbind(result_RNA$aov.tab, result_DNA$aov.tab))
#write.csv(results, "permanova_bray_DNA_RNA_subsets.csv")

######################################################################################################
# dispersion analysis across sampling time (Month)
######################################################################################################

# calculate the impact of dataset of each year ----------
permanova.list <- list()
dispersion.list <- list()
dispersion.sites.list <- list()

months <- c("5", "7", "9", "11")
#mo <- "7"
for (mo in months) {
  com_sub <- subset(com, grepl(paste("_", mo, "_", sep = ""), row.names(com)))
  try(if(nrow(com_sub) !=30) stop("Error in subset com!"))
  df <- data.frame(row.names=rownames(com_sub), t(as.data.frame(strsplit(rownames(com_sub),"_"))))
  df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))
  head(df)
  
  cat("\nTwo-way permanova (Dataset vs. Year) for Month:", mo)
  set.seed(123)
  result_whole <- adonis(com_sub ~ Dataset*Year, data=df, method="bray", permutation=9999) # two way permanova
  result_whole$aov.tab
  permanova.list[[mo]] <- as.data.frame(result_whole$aov.tab)
  
  cat("\nDispersion between DNA & RNA for Month:", mo) 
  dist <- vegdist(com_sub, method="bray", binary=FALSE, diag=1) 
  dispersion <- betadisper(dist, group = df$Dataset)
  dispersion.list[[mo]] <- permutest(dispersion, permutation=9999)$tab
  p <- plot(dispersion, hull=FALSE, ellipse=TRUE, main = paste("Dispersion (", mo, " mo)", sep = "")) ##sd ellipse
  dispersion.sites.list[[mo]] <- p$sites
}

str(permanova.list)
df <- do.call(rbind.data.frame, permanova.list)
df$Year <- gsub("\\..*", "", row.names(df))
df$ii <- sub('.*\\.', '', row.names(df))
write.csv(df, "permanova_dataset_months.csv")

str(dispersion.list)
df <- do.call(rbind.data.frame, dispersion.list )
df$Year <- gsub("\\..*", "", row.names(df))
df$ii <- sub('.*\\.', '', row.names(df))
write.csv(df, "dispersion_permanova_dataset_months.csv")

str(dispersion.sites.list )
dispersion.sites <- do.call(rbind.data.frame, dispersion.sites.list )
dispersion.sites$Year <- gsub("\\..*", "", row.names(dispersion.sites))
dispersion.sites$sample <- sub('.*\\.', '', row.names(dispersion.sites))
write.csv(df, "dispersion_sites_dataset_months.csv")

# all communities
df <- data.frame(row.names=dispersion.sites$sample, t(as.data.frame(strsplit(dispersion.sites$sample,"_"))))

df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))

df$Dataset <- factor(df$Dataset, levels = c("cDNA", "DNA"), labels = c("rRNA", "rDNA"))

df <- data.frame(df, dispersion.sites[, 1:2])
head(df)

(my_palette <- c(brewer.pal(9, "Set1")[c(1:2)]))

(f2 <- ggplot(df, aes(PCoA1, PCoA2, shape=Year, fill=Dataset))+
    facet_grid(~Month, scales = "free_x") +
    geom_point(size=4, alpha=0.7)+ 
    labs(x= "PCoA1", y = "PCoA2", title = " ")+
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(24, 25, 22, 23, 21))+
    mytheme
)

(f <- ggarrange(f1, ggarrange(f2, NULL, ncol = 2,  widths = c(2, 0.485), heights = c(1, 1)), labels = c("A", "B"), common.legend = FALSE, legend = "right", 
             ncol = 1, nrow = 2, align = "v"))

ggsave("PCoA_dispersion_dataset.png", width = 14, height = 7.5, units = "cm", f, scale = 2, dpi = 300)

####################################################################################################
# procrustes analysis
####################################################################################################
plot(procrustes(rda(sptrans), cca(varespec)))

