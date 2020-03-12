# alpha diversity analysis for the DMA/RNA dataset
# Date: 08-08-2019 
# Update: 19-09-2019
# Author: Jia, Xiu

rm(list=ls())

# load the directory
directory = 'C:/Users/P278113/Dropbox'
#directory = '~/Dropbox'
subfolder = 'Schier/DNA_RNA'

setwd(paste(directory, subfolder, sep="/"))
getwd()

# load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
display.brewer.all()
library(reshape2)
library(doBy) 
library(plyr)

# ggplot theme
mytheme <- theme_bw()+
  theme(text = element_text(size = 12),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text=element_text(face="bold", size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# load the rarefied otu table -----------------------------------------------------------------------
com <- read.csv("feature-table-rarified.csv", sep=",", header=1, row.names=1, check.names = FALSE)
com <- t(com[, 1:120])

# calculate alpha-diversity  ------------------------------------------------------------------------
pd <- read.csv("PD_DNA_RNA.csv", sep=",", header=1, row.names=1, check.names = FALSE); head(pd)
shannon <- diversity(com)
simpson <- diversity(com, "simpson")
richness <- specnumber(com)
Chao <- estimateR(com)
Chao1 <- as.data.frame(t(Chao))
# Pielou's J = H'/ln(S) where H' is Shannon Weiner diversity and S is the total number of species in a sample,
pielous.evenness <- shannon/log(richness)

group_info<-data.frame(row.names=rownames(com), t(as.data.frame(strsplit(rownames(com),"_"))))
head(group_info)

df <- data.frame(cbind(richness, Chao1, shannon, simpson, pielous.evenness),
                 PD = pd$PD,
                 Datasets = as.factor(group_info[,1]),
                 Year = as.factor(group_info[,2]),
                 Month = as.factor(group_info[,3]),
                 replicates = as.factor(group_info[,4]))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))
df$Year <- as.numeric(as.character(df$Year))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
df$Datasets <- factor(df$Datasets, levels = c("cDNA", "DNA"), labels = c("RNA", "DNA"))
df$ro1 <- df$richness*100/df$S.chao1
df$ro2 <- df$richness*100/df$S.ACE
str(df)
head(df)


## Coefficient of variation (CV) --------------------------------------------------------------------

# Coefficient of Variation = (Standard Deviation / Mean) * 100.
df.richness <- ddply(df, c("Datasets", "Year"), summarise,
            N    = length(richness),
            mean = mean(richness),
            sd   = sd(richness),
            se   = sd / sqrt(N),
            cv   = sd*100/mean,
            div  = "richness")

df.shannon <- ddply(df, c("Datasets", "Year"), summarise,
            N    = length(shannon),
            mean = mean(shannon),
            sd   = sd(shannon),
            se   = sd / sqrt(N),
            cv   = sd*100/mean,
            div  = "shannon")

df.PD <- ddply(df, c("Datasets", "Year"), summarise,
            N    = length(PD),
            mean = mean(PD),
            sd   = sd(PD),
            se   = sd / sqrt(N),
            cv   = sd*100/mean,
            div  = "PD")

df.pielous.evenness <- ddply(df, c("Datasets", "Year"), summarise,
            N    = length(pielous.evenness),
            mean = mean(pielous.evenness),
            sd   = sd(pielous.evenness),
            se   = sd / sqrt(N),
            cv   = sd*100/mean,
            div  = "pielous.evenness")

data <- rbind(df.richness, df.shannon, df.PD, df.pielous.evenness)
head(data)

#write.csv(data, "coefficient_of_variation_alpha_div.csv")

# Line plot for four alpha-diversity measurements --------------------------------------------------
data$div <- factor(data$div, levels = c("richness", "shannon", "pielous.evenness", "PD"),
                       labels = c("Richness", "Shannon", "Pylogenetic diversity", "Pielou's evenness"))

pd <- position_dodge(0.8)  
(p <- ggplot(data, aes(x = Year, y = mean, group = Datasets, colour = Datasets)) + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=3, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, alpha=0.9) +
    scale_color_brewer(palette="Set1") +
    scale_x_continuous(breaks = c(0, 10, 40, 70, 110))+
    facet_wrap(~div, nrow=2, scales="free_y") +
    labs(y="", x = "Stage of succession (Year)")+
    mytheme)

#ggsave("alpha_div_four_RNA_DNA.jpg", width = 16, height = 12, units = "cm", p, scale = 1, device = "jpeg",dpi = 600)

# line plot for Sampling depth -------------------------------------------------------------------
df.ro1 <- ddply(df, c("Datasets", "Year"), summarise,
                N    = length(ro1),
                mean = mean(ro1),
                sd   = sd(ro1),
                se   = sd / sqrt(N),
                cv   = sd*100/mean,
                div  = "Richness/Chao1")

df.ro2 <- ddply(df, c("Datasets", "Year"), summarise,
                N    = length(ro2),
                mean = mean(ro2),
                sd   = sd(ro2),
                se   = sd / sqrt(N),
                cv   = sd*100/mean,
                div  = "Richness/ACE")

data2 <- rbind(df.ro1, df.ro2)
#data2$Year <- factor(data2$Year, levels=c("0", "10", "40", "70", "110"))

str(data2)
pd <- position_dodge(0)
(p <- ggplot(data2, aes(x = Year, y = mean, group = Datasets, colour = Datasets)) + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.8, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3) +
    scale_color_brewer(palette="Set1") +
    #scale_x_continuous(breaks = c(0, 10, 40, 70, 110))+
    facet_wrap(div~., nrow=1, scales="free_y") +
    labs(y="Sampling depth", x = "Stage of succession (Year)")+
    mytheme)

#ggsave("sampling_depth_RNA_DNA.jpg", width = 16, height = 7, units = "cm", p, scale = 1, device = "jpeg",dpi = 600)

# box plot for all diversity indexes
df1 <- melt(df[, -c(2, 4, 6, 14)], id.vars = c("Datasets", "Year", "Month"))
df1$Year <- factor(df1$Year, levels=c("0", "10", "40", "70", "110"))
df1$variable <- factor(df1$variable, levels = c("richness", "S.chao1","S.ACE", "PD", "shannon", "simpson", "pielous.evenness"),
                       labels = c("Richness", "Chao1","ACE", "Pylogenetic diversity", "Shannon", "Simpson", "Pielou's evenness"))
head(df1)

(my_palette = c(brewer.pal(9, "Set1")[c(1,2)]))
( p <- ggplot(df1, aes(x = Year, y = value, fill = Datasets))+ 
    geom_point(size = 2, shape = 16, alpha = 0.3, position = position_jitterdodge()) +
    geom_boxplot(alpha = 0.8, outlier.size=-1) +
    facet_wrap(.~variable, scales = "free_y") +
    scale_fill_manual(values = my_palette) + 
    labs(x="Stage of succession (Years)", y="", title="") + 
    mytheme)
#ggsave("alpha_div_all_boxplot.jpg", width = 16, height = 12, units = "cm", device = "jpeg", p, scale = 1.5, dpi = 300)

# statistic test -------------------------------------------------------------------------------------
# check the homogeneity of variances
bartlett.test(richness ~ interaction(Year, Month), data=df[which(df$Datasets == "RNA"),])
# p > 0.05, suggesting that the variance across groups is not statistically significantly different.

# Shapiro-Wilk normality test for each dataset
shapiro.test(df[which(df$Datasets == "RNA"),]$simpson) 

# only RNA dataset is normally distributed, p-value=0.01733; DNA pvalue= 0.0001144
ggqqplot(df[which(df$Datasets == "RNA"),]$richness, ylab = "richness")


#method=c("pearson", "kendall", "spearman")
#Kendall's rank correlation tau
#tau is the Kendall correlation coefficient
(res1 <- cor.test(df[which(df$Datasets == "DNA"),]$richness, 
                  df[which(df$Datasets == "DNA"),]$Year,  method="kendall"))
(res2 <- cor.test(df[which(df$Datasets == "RNA"),]$richness, 
                  df[which(df$Datasets == "RNA"),]$Year,  method="pearson"))

(res3 <- cor.test(df[which(df$Datasets == "DNA"),]$shannon, 
                  df[which(df$Datasets == "DNA"),]$Year,  method="kendall"))
(res4 <- cor.test(df[which(df$Datasets == "RNA"),]$shannon, 
                  df[which(df$Datasets == "RNA"),]$Year,  method="kendall"))

# correlation analysis
(my_palette = c(brewer.pal(9, "Set1")[c(1,2)]))

# loess fit plots with density plots----------------------------------------------------------------------------
(p1 <- ggplot(df, aes(x = Year, y = richness, color = Datasets, shape = Datasets)) + 
   geom_point(size = 3, alpha = 0.7)+
   geom_smooth(aes(fill = Datasets), method = loess,  linetype="dashed") +
   scale_color_manual(values = my_palette) + 
   scale_size_manual(values = c(16, 15)) +
   scale_x_continuous(breaks = c(0, 10, 40, 70, 110))+
   annotate('text', label = paste("Pearson's r = ", round(res2$statistic, 3), ", p-value = ", round(res2$p.value,3), sep = ""), 
            x = -Inf, y = -Inf, size = 3, hjust = -.24, vjust = -2.1, color = "#E41A1C") +
   annotate('text', label = paste("Kendall's tau = ", round(res1$statistic, 3), ", p-value = ", round(res1$p.value,3), sep = ""), 
            x = -Inf, y = -Inf, size = 3, hjust = -.23, vjust = -.8, color = "#377EB8") +
   labs(x="Stage of succession (Years)", y="Richness")+
   mytheme)

# Marginal density plot of y 
(f1 <- ggplot(df, aes(richness, fill=Datasets)) + 
    geom_density(alpha=.6) + 
    scale_fill_manual(values=my_palette) + #values = c('#999999','#E69F00')) + # some nice colour
    labs(x = "Richness", y = "Density")+
    mytheme +
    theme(legend.position="none"))

(p2 <- ggplot(df, aes(x = Year, y = simpson, color = Datasets, shape = Datasets)) + 
    geom_point(size = 3, alpha = 0.7)+
    geom_smooth(aes(fill = Datasets), method = loess, linetype="dashed") +
    scale_color_manual(values = my_palette) + 
    scale_size_manual(values = c(16, 15)) +
    scale_x_continuous(breaks = c(0, 10, 40, 70, 110))+
    annotate('text', label = paste("Kendall's tau = ", round(res4$statistic, 3), ", p-value = ", round(res4$p.value,3), sep = ""), 
             x = -Inf, y = -Inf, size = 3, hjust = -.25, vjust = -2.1, color = "#E41A1C") +
    annotate('text', label = paste("Kendall's tau = ", round(res3$statistic, 3), ", p-value = ", round(res3$p.value,4), sep = ""), 
             x = -Inf, y = -Inf, size = 3, hjust = -.25, vjust = -.8, color = "#377EB8") +
    labs(x="Stage of succession (Years)", y="Simpson")+
    mytheme)

# Marginal density plot of y 
(f2 <- ggplot(df, aes(simpson, fill = Datasets)) + 
    geom_density(alpha=.6) + 
    scale_fill_manual(values=my_palette) + #values = c('#999999','#E69F00')) + 
    labs(x = "Simpson", y = "Density")+
    mytheme +
    theme(legend.position="none"))

(f <- annotate_figure((ggarrange(f1, f2, p1, p2, labels = c("", "", "A", "B"), ncol = 2, nrow = 2, 
                                 heights = c(0.7, 2),  align = "v", common.legend = TRUE, legend = "right")), right = " "))

#ggsave("Correlation_density_rich_sim_dataset_year.jpg", width = 14, height = 7, units = "cm", device = "jpeg", f, scale = 1.5, dpi = 300)
#ggsave("Correlation_density_rich_sim_dataset_year.pdf", width = 15, height = 7.5, units = "cm", f, scale = 1.5, dpi = 300)

# Stats for each stage
(res <- wilcox.test(simpson ~ Datasets, data = df[which(df$Year == "110"),], exact = FALSE))


# FORCE TO DO Pearson correlation ------------------------------------------------------------------
(p1 <- ggscatter(df, x = "Year", y = "richness", color = "Datasets",
                 add = "reg.line", conf.int = TRUE) + 
   #stat_cor(aes(color = Datasets), label.x = 3) +
   scale_color_manual(values=my_palette) + 
   scale_x_continuous(breaks = c(0, 10, 40, 70, 110))+
   #stat_cor(aes(color = Datasets, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3) +
   mytheme)

(p2 <- ggscatter(df, x = "Year", y = "simpson", color = "Datasets",
                 add = "reg.line", conf.int = TRUE) + 
    #stat_cor(aes(color = Datasets), label.x = 3) +
    scale_color_manual(values=my_palette) + 
    scale_x_continuous(breaks = c(0, 10, 40, 70, 110))+
    #stat_cor(aes(color = Datasets, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3) +
    mytheme)

f <- ggarrange(p1, p2, labels = c("A", "B"), common.legend = TRUE, legend = "right", ncol = 2)
f

#ggsave("Correlation_alpha_div_year.png", width = 17, height = 7, units = "cm", f, scale = 1.5, dpi = 300)
#ggsave("Correlation_alpha_div_year.pdf", width = 17, height = 7, units = "cm", f, scale = 1.5)



# Two-way ANOVA -----------------------------------------------------------------------------------------

my_palette <- c(brewer.pal(9, "Set1")[c(1:2)])

# Line plots with multiple groups
(f1 <- ggline(df, x = "Year", y = "richness", 
              linetype = "Datasets", shape = "Datasets", color = "Datasets",
              point.size = 3,
              add = c("mean_se", "dotplot"),
              palette = my_palette)+
    labs(x="Stage of succession (Years)", y="Richness", title=" ") +
    mytheme)

(f2 <- ggline(df, x = "Year", y = "simpson", 
              linetype = "Datasets", shape = "Datasets", color = "Datasets",
              point.size = 3,
              add = c("mean_se", "dotplot"),
              palette = my_palette)+
    labs(x="Stage of succession (Years)", y="simpson", title=" ") +
    mytheme)

f <- ggarrange(f1, f2, labels = c("A", "B"), common.legend = TRUE, legend = "right", ncol = 2)
f

#ggsave("Richness_simpson_Datasets_x_year_lineplot.png", width = 17, height = 7.5, units = "cm", f, scale = 1.5, dpi = 300)
#ggsave("Richness_simpson_Datasets_x_year_lineplot.pdf", width = 17, height = 7.5, units = "cm", f, scale = 1.5)

# Two-way interaction plot
interaction.plot(x.factor = df$Year, trace.factor = df$Datasets, 
                 response = df$richness, fun = mean, 
                 type = "b", legend = TRUE, 
                 xlab = "Stage of succession (Years)", ylab="Richness",
                 pch=c(1,19), col = my_palette)


# richness, shannon, simpson and PD
aov <- aov(richness ~ Datasets * Year, data = df) # NO interaction effect
summary(aov)
capture.output(summary(aov), file = "anova_results_richness.csv")

# check the homogeneity of variances
bartlett.test(richness ~ interaction(Datasets, Year), data=df)
# p > 0.05, suggesting that the variance across groups is not statistically significantly different.
plot(aov, 1)
# there is no evident relationships between residuals and fitted values (the mean of each groups), which is good. 
# So, we can assume the homogeneity of variances.

# check normality
plot(aov, 2) # draw the correlation between a given sample and the normal distribution.
ggqqplot(df$richness)

# Shapiro-Wilk test of normality for univariate
shapiro.test(df$richness)
# p-value > 0.05, implying that the distribution of the data are not significantly different from normal distribution. 


