                                    ######### PHYLOSEQ pipeline for visualizing the diversity  ###########
                                    
                                    follow tutorial from https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
## install PHYLOSEQ package.
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")                                   
install.packages("dplyr")
install.packages("readxl") 
install.packages("ggplot2") 
install.packages("SpiecEasi") 
install.packages("ggpubr") # calculate T.test (stat_compare_means)
install.packages("microbiome") 
install.packages("vegan")
install.packages("utils")
## loading libraries
library(phyloseq) # processing the data.
packageVersion("phyloseq")
library(ggplot2) # for visualization (for high quality graphics) 
library(vegan) # for analysis
library(igraph) # for network analysis
library(dplyr) # To manipulate dataframes (filter and reformat data frames)
library(readxl) # To read Excel files into R
library(RColorBrewer)
library(remotes)
library(SpiecEasi)
library(Biostrings)
library("tibble") # Needed for converting column to row names
library(ggpubr)
theme_set(theme_bw())
library(vegan)
library(utils)             

# Taxa information
taxonomyGTDB44 <- read.csv("E:/sequencing data/microbiota/dada2 output/GTDB/Taxonomy.GTDB44.csv")
HowManyKingdom <- length(unique(unlist(taxonomyGTDB44[,c("Kingdom")]))) # 3
HowManyPhyla <- length(unique(unlist(taxonomyGTDB44[,c("Phylum")]))) # 104
HowManyClass <- length(unique(unlist(taxonomyGTDB44[,c("Class")]))) # 227
HowManyOrder <- length(unique(unlist(taxonomyGTDB44[,c("Order")]))) # 573
HowManyFamily <- length(unique(unlist(taxonomyGTDB44[,c("Family")]))) # 1126
HowManyGenus <- length(unique(unlist(taxonomyGTDB44[,c("Genus")]))) # 2608
HowManySpecies <- length(unique(unlist(taxonomyGTDB44[,c("Species")]))) # 3314


                           #################### Microbial diversity estimates ######################
##### (1) Alpha diversity ######
## create phyloseq file (long file 27872) for alpha diversity measure
asv_mat<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_Alpha27872.xlsx", sheet = "ASV_table"))
tax_mat<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_Alpha27872.xlsx", sheet = "taxa_table"))
sample_df<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_Alpha27872.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_mat)<-asv_mat$asv
asv_mat<-asv_mat%>% select(-asv)

row.names(tax_mat)<-tax_mat$asv
tax_mat<-tax_mat%>% select(-asv)

row.names(sample_df)<-sample_df$sample
sample_df<-sample_df%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_mat<-as.matrix(asv_mat)
tax_mat<-as.matrix(tax_mat)

## Create phyloseq objects.
(A)
ASV = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX  = tax_table(tax_mat)
samples = sample_data(sample_df)
metabarcoding <-phyloseq(ASV, TAX, samples)
metabarcoding
(B)
data <- phyloseq(otu_table(asv_mat, taxa_are_rows = TRUE), sample_data(sample_df))
tax_table(data) <- as.matrix(tax_mat)
data

## data preprocessing (sums, maximum, minimum)
reads <- sample_sums(metabarcoding)  # number of reads/seqs. per sample
max(sample_sums(metabarcoding))
min(sample_sums(metabarcoding))
length(which(reads<5000))
reads
## Normalize counts by the sample sums (Fractions).
otu <- data.frame(otu_table(data))
sum <- apply(otu, 2, sum)
norm_otu <- otu
for (i in 1:ncol(otu)){
  norm_otu[,i] <- otu[,i]/sum[i]
}
normalized_fraction <- phyloseq(otu_table(norm_otu, taxa_are_rows = TRUE), tax_table(as.matrix(tax_mat)), sample_data(sample_df))
normalized_fraction
 total = median(sample_sums(metabarcoding))  # find median sample read count (94861).
standf = function(x, t=total) round(t * (x / sum(x)))  # function to standardize to median sample read count.
normalized_median = transform_sample_counts(metabarcoding, standf)  # apply to phyloseq object.
normalized_median
ntaxa(normalized_median)
sample_sums(normalized_median)
## Plot Chao1 richness estimator and Shannon diversity estimator.
plot_richness(normalized_median, measures=c("Chao1", "Shannon"))
sample <- plot_richness(normalized_median, x="sample", color="layer", measures=c("Chao1", "Shannon"), 
              title = paste0("Standardized to total reads, N=", nsamples(normalized_median)))
sample 
area <- plot_richness(normalized_median, x="area", color="layer", measures=c("Observed", "Chao1", "Shannon"), 
              title = paste0("Standardized to total reads, N=", nsamples(normalized_median)))
area 
layer <- plot_richness(normalized_median, x="layer", color="layer", measures=c("Observed", "Chao1", "Shannon"), 
                       title = paste0("Standardized to total reads, N=", nsamples(normalized_median)))
layer
layer +  stat_compare_means(method = "t.test", label.y = 0.8)
kruskal.test(layer~Observed, data = div) # p-value= 0.5557
kruskal.test(layer~Chao1, data = div) # p-value= 0.5557
kruskal.test(layer~Shannon, data = div) # p-value= 0.4713

est <- estimate_richness(normalized_median, split = TRUE, measures = NULL)
tem <- cbind(sample_data(normalized_median), est)
aggregate(tem$Chao1, list(tem$location), mean)  # calculate the means
write.csv(tem, file = "diversity_normalized_richness GTDB44.csv")

# violin plot by ggstatsplot package
install.packages("tidyverse")
install.packages("ggstatsplot")
install.packages("palmerpenguins")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggpubr")
library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)

div <- read.csv("E:/sequencing data/microbiota/dada2 output/GTDB/diversity_normalized_richness GTDB44.csv")
div <- drop_na(div)
observed <- ggbetweenstats(div, x = layer, y = Observed)
observed
observed <- observed + labs(x = "sediment_layer", y = "observed_ASVs", title = paste0("Standardized to total reads, N=", nsamples(normalized_median))) +  
  theme( text = element_text(family = "Roboto", size = 8, color = "black"), plot.title = element_text(family = "Lobster Two", 
  size = 20, face = "bold", color = "#2a475e"), plot.subtitle = element_text(family = "Roboto", size = 15, face = "bold", color="#1b2838"),
  plot.title.position = "plot", axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 12))
observed <- observed + theme(axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
                       panel.grid = element_line(color = "#b4aea9"),
                       panel.grid.minor = element_blank(),
                       panel.grid.major.x = element_blank(),
                       panel.grid.major.y = element_line(linetype = "dashed"),
                       panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
                       plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))
observed

chao1 <- ggbetweenstats(div, x = layer, y = Chao1)
chao1
chao1 <- chao1 + labs(x = "sediment_layer", y = "Chao1_Richness", title = paste0("Standardized to total reads, N=", nsamples(normalized_median))) +  
  theme( text = element_text(family = "Roboto", size = 8, color = "black"), plot.title = element_text(family = "Lobster Two", 
  size = 20, face = "bold", color = "#2a475e"), plot.subtitle = element_text(family = "Roboto", size = 15, face = "bold", color="#1b2838"),
  plot.title.position = "plot", axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 12))
chao1 <- chao1 + theme(axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
                   panel.grid = element_line(color = "#b4aea9"),
                   panel.grid.minor = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.grid.major.y = element_line(linetype = "dashed"),
                   panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
                   plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))
chao1

shannon <- ggbetweenstats(div, x = area, y = Shannon)
shannon <- shannon + labs(x = "area", y = "Shannon_estimates", title = paste0("Standardized to total reads, N=", nsamples(normalized_median))) +  
  theme( text = element_text(family = "Roboto", size = 8, color = "black"), plot.title = element_text(family = "Lobster Two", 
  size = 20, face = "bold", color = "#2a475e"), plot.subtitle = element_text(family = "Roboto", size = 15, face = "bold", color="#1b2838"),
  plot.title.position = "plot", axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 12))
shannon <- shannon + theme(axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
                     panel.grid = element_line(color = "#b4aea9"),
                     panel.grid.minor = element_blank(),
                     panel.grid.major.x = element_blank(),
                     panel.grid.major.y = element_line(linetype = "dashed"),
                     panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
                     plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))
shannon

# Box plot
boxplot3 <- ggplot(div, aes(x=layer, y=Shannon)) + geom_boxplot() + geom_jitter(aes(color=layer)) +xlab("Sediment layer") + ylab("Shannon_estimates") + ggtitle("Standardized to total reads, N=44")
boxplot3
boxplot3 + stat_compare_means(method = "t.test", label.y = 0.8)

######## (2) Beta diversity ######
https://astrobiomike.github.io/amplicon/dada2_workflow_ex#beta-diversity

asv_mat<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_Alpha29473.xlsx", sheet = "ASV_table"))
tax_mat<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_Alpha29473.xlsx", sheet = "taxa_table"))
sample_df<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_Alpha29473.xlsx", sheet = "sample"))

count <- read.table("E:/sequencing data/microbiota/dada2 output/GTDB/shiangshan_16s_count_tab.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("E:/sequencing data/microbiota/dada2 output/GTDB/shiangshan_16s_tax_tab.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
sample_info_tab <- read.table("E:/sequencing data/microbiota/dada2 output/GTDB/shiangshan_16s_sample_info_tab.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")
sample_info_tab$color <- as.character(sample_info_tab$color)
 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
deseq_counts <- DESeqDataSetFromMatrix(count, colData = sample_info_tab, design = ~layer)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count <- assay(deseq_counts_vst)
euc_dist<-dist(t(vst_trans_count))

euc_clust <- hclust(euc_dist, method = "ward.D2")
plot(euc_clust, main="Cluster of Siangshan Samples")
rect.hclust(tree =euc_clust, k = 4, border = "red")

# running betadisper and permutational ANOVA test, first we need to make a phyloseq object with our transformed data
# --> the 'adonis' function is to perform a permutational ANOVA to see if there is a statistical difference in community structure across processing steps (sieved v nonsieved)
# --> but prior to running the adonis test, a betadisper test must be run to check that there is no difference between within-group dispersion
anova(betadisper(euc_dist, sample_info_tab$layer))
# if the anova test returns a non-significant p-value, then we can procedd and run the betadisper test 
adonis(euc_dist~sample_info_tab$layer)
# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
vst_physeq


# plot ordination graph
set.seed(20)
Siangshan_ord <- ordinate(normalized_median, "NMDS", "bray", k=4, trymax=100)
Siangshan_ord
# Call:
# metaMDS(comm = veganifyOTU(physeq), distance = distance, k = 2,      trymax = 100) 
# global Multidimensional Scaling using monoMDS
# Data:     wisconsin(sqrt(veganifyOTU(physeq))) 
# Distance: bray 
# Dimensions: 2 
# Stress:     0.150205 
# Stress type 1, weak ties
# Two convergent solutions found after 20 tries
# Scaling: centring, PC rotation, halfchange scaling 
# Species: expanded scores based on 'wisconsin(sqrt(veganifyOTU(physeq)))' 
stressplot(Siangshan_ord)
legend("topleft", "stress = 0.163", bty = "n", cex = 1)
goodness(Siangshan_ord)
[1] 0.016910206 0.019400639 0.016467265 0.020920543 0.014311920 0.024108916 0.011935666 0.024867342 0.009881586 0.020690863 0.010122143 0.011975328 0.027051805 0.012747413
[15] 0.023012771 0.010689225 0.022716389 0.013719627 0.048121753 0.010431576 0.013717733 0.011755105 0.020279143 0.019221610 0.032846691 0.017378956 0.028169973 0.031935007
[29] 0.022948560 0.033487378 0.026822089 0.028703615 0.021640781 0.016939743 0.018723744 0.020139444 0.028189512 0.032594688 0.024340644 0.014654760 0.014741414 0.013009011
[43] 0.013613416 0.014617894 0.013873367 0.012745840 0.028663113 0.012019052 0.020733452
color = c("layer")
shape = c("layer")
title = c("NMDS of 16s rRNA microbiome, Bray-Curtis distance, K=2")
Ordination <- plot_ordination(normalized_median, Siangshan_ord, type="sample", color=color, 
                shape=shape, title=title) 
ordination = Ordination + theme(axis.text = element_text(size = 16, face = "bold"), axis.title = element_text(size = 18, face = "bold"), 
             legend.title = element_text(size = 14)) + theme_bw() + labs(color = color, shape = shape) + geom_point(size=2) + 
geom_text(aes(label=rownames(sample_df), hjust=0.4, vjust=-0.5)) 
ordination 

pc = read.csv("E:/sequencing data/microbiota/dada2 output/GTDB/ASVs_countsGTDB44.csv", header= TRUE)
com = pc[,1:ncol(pc)]
m_com = as.matrix(com)
ano = anosim(m_com, pc$sample, distance = "bray", permutations = 9999)
ano
####### (3) Rarefaction curve ######
# three separate files are required for this pipeline, an ASV count table, a sample info table and a taxonomy table (these can be made in Excel from teh insect output ASV table, similar to the phyloseq pipeline requirements)
asv_mat<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_Alpha29473.xlsx", sheet = "ASV_table"))
tax_mat<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_Alpha29473.xlsx", sheet = "taxa_table"))
sample_df<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_Alpha29473.xlsx", sheet = "sample"))

install.packages("permute")
install.packages("lattice")
library("vegan")
library("permute")
library("lattice")

# import the count table
Data=read.table("E:/sequencing data/microbiota/dada2 output/GTDB/ASVs_countsGTDB.csv")
Data_sub=as.data.frame(Data[,c(1:49)])
Data_t=(t(Data_sub))
# count the number of species
S <- specnumber(Data_t)
raremax<-min(rowSums(Data_t)) # the minimum sample count achieved over the 26 samples. We will rarefy the sample counts to this value
# Rarefaction of the samples
Srare <- rarefy(Data_t, raremax)
# plotting the rarefaction curves
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(Data_t, step =1, sample = raremax, col = "blue", cex = 0.6,)
title("Rarefaction curve of Siangshan Wetland")
(1)
Data=read.table("E:/sequencing data/microbiota/dada2 output/GTDB/ASVs_countsGTDB.csv")
Data_sub=as.data.frame(Data[,c(1:49)])
Data_t=(t(Data_sub))
raremax <- min(rowSums(Data_t))
raremax # 3405
col <- rainbow(49)
lty <- c("solid")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)
rare <- with(pars[1:49, ], rarecurve(Data_t, step = 20, sample = raremax, col = col, lty = lty, label = TRUE))
col <- rainbow(49)
lty <- c("solid")
lwd <- c(1)
pars <- expand.grid(col = col, lty = lty, lwd = lwd, stringsAsFactors = FALSE)
Nmax <- sapply(rare, function(x) max(attr(x, "Subsample")))
Smax <- sapply(rare, max)
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size", ylab = "ASVs", type = "n")
abline(v = raremax)
for (i in seq_along(rare)) {
  N <- attr(rare[[i]], "Subsample")
  with(pars, lines(N, rare[[i]], col = col[i], lty = lty[i], lwd = lwd[i]), lable = TRUE)
}
title("Rarefaction curve of Siangshan Wetland")
(2)
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)
rarefaction <- with(pars[1:49,], rarecurve(Data_t, step = 20, sample = raremax, col = col, lty = lty, label = FALSE))
col <- c("black", "darkred", "forestgreen", "hotpink", "blue")
lty <- c("solid", "dashed", "dotdash")
lwd <- c(1, 2)
pars <- expand.grid(col = col, lty = lty, lwd = lwd, stringsAsFactors = FALSE)
head(pars)
Nmax <- sapply(rarefaction, function(x) max(attr(x, "Subsample")))
Smax <- sapply(rarefaction, max)
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size", ylab = "ASVs", type = "n")
abline(v = raremax)
for (i in seq_along(rarefaction)) {
  N <- attr(rarefaction[[i]], "Subsample")
  with(pars, lines(N, rarefaction[[i]], col = col[i], lty = lty[i], lwd = lwd[i]))
}
(3)
rarefied_microbes <- rarefy_even_depth(data)
rarecurve(t(otu_table(normalized_median)), step = 50, cex = .6, col = "red")
rarecurve(data, step = 20, col = "blue", cex = 0.6, label = FALSE)

#################################################################################################################################################
                                      ################## Visualize the relative abundance #################

library(hrbrthemes)
library(gcookbook)
library(tidyverse)

(1) surface sediment
asv_mat_surf<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_surface176.xlsx", sheet = "ASV_table"))
tax_mat_surf<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_surface176.xlsx", sheet = "taxa_table"))
sample_df_surf<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_surface176.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_mat_surf)<-asv_mat_surf$asv
asv_mat_surf<-asv_mat_surf%>% select(-asv)

row.names(tax_mat_surf)<-tax_mat_surf$asv
tax_mat_surf<-tax_mat_surf%>% select(-asv)

row.names(sample_df_surf)<-sample_df_surf$sample
sample_df_surf<-sample_df_surf%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_mat_surf<-as.matrix(asv_mat_surf)
tax_mat_surf<-as.matrix(tax_mat_surf)

## Create phyloseq objects. 
(A)
ASV_surf = otu_table(asv_mat_surf, taxa_are_rows = TRUE)
TAX_surf  = tax_table(tax_mat_surf)
samples_surf = sample_data(sample_df_surf)
metabarcoding_surf <-phyloseq(ASV_surf, TAX_surf, samples_surf)
metabarcoding_surf
(B)
data_surf <- phyloseq(otu_table(asv_mat_surf, taxa_are_rows = TRUE), sample_data(sample_df_surf))
tax_table(data_surf) <- as.matrix(tax_mat_surf)
data_surf

## Normalize counts by the sample sums (Fractions).
otuF <- data.frame(otu_table(data_surf))
sumF <- apply(otuF, 2, sum)
norm_otuF <- otuF
for (i in 1:ncol(otuF)){
  norm_otuF[,i] <- otuF[,i]/sumF[i]
}
normalized_surf <- phyloseq(otu_table(norm_otuF, taxa_are_rows = TRUE), tax_table(as.matrix(tax_mat_surf)), sample_data(sample_df_surf))
normalized_surf

plot_bar(normalized_surf, x = 'location', fill = "Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity", position = "stack") +
  labs(x='location', y='Relative abundance') + theme_bw() + scale_y_percent() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))


(1) surface sediment
asv_mat_surf<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_most abund_surf.xlsx", sheet = "ASV_table"))
tax_mat_surf<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_most abund_surf.xlsx", sheet = "taxa_table"))
sample_df_surf<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_most abund_surf.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_mat_surf)<-asv_mat_surf$asv
asv_mat_surf<-asv_mat_surf%>% select(-asv)

row.names(tax_mat_surf)<-tax_mat_surf$asv
tax_mat_surf<-tax_mat_surf%>% select(-asv)

row.names(sample_df_surf)<-sample_df_surf$sample
sample_df_surf<-sample_df_surf%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_mat_surf<-as.matrix(asv_mat_surf)
tax_mat_surf<-as.matrix(tax_mat_surf)

## Create phyloseq objects. 
(A)
ASV_surf = otu_table(asv_mat_surf, taxa_are_rows = TRUE)
TAX_surf  = tax_table(tax_mat_surf)
samples_surf = sample_data(sample_df_surf)
metabarcoding_surf_abun <-phyloseq(ASV_surf, TAX_surf, samples_surf)
metabarcoding_surf_abun
(B)
data_surf_abun <- phyloseq(otu_table(asv_mat_surf, taxa_are_rows = TRUE), sample_data(sample_df_surf))
tax_table(data_surf_abun) <- as.matrix(tax_mat_surf)
data_surf_abun

## Normalize counts by the sample sums (Fractions).
otuF <- data.frame(otu_table(data_surf_abun))
sumF <- apply(otuF, 2, sum)
norm_otuF <- otuF
for (i in 1:ncol(otuF)){
  norm_otuF[,i] <- otuF[,i]/sumF[i]
}
normalized_surf_abun <- phyloseq(otu_table(norm_otuF, taxa_are_rows = TRUE), tax_table(as.matrix(tax_mat_surf)), sample_data(sample_df_surf))
normalized_surf_abun

plot_bar(normalized_surf_abun, x = 'location', fill = "Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity", position = "stack") +
  labs(x='location', y='Relative abundance') +
theme_bw() + scale_y_percent() + ggtitle("Surface samples (Normalized by sum, N = 26)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))

plot_bar(normalized_surf_abun, x = 'location', fill = "Class") + geom_bar(aes(color=Class, fill=Class), stat = "identity", position = "stack") +
  labs(x='location', y='Relative abundance') + theme_bw() + scale_y_percent() + ggtitle("Surface samples (Normalized by sum, N = 26)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))

plot_bar(normalized_surf_abun, x = 'location', fill = "Genus") + geom_bar(aes(color=Genus, fill=Genus), stat = "identity", position = "stack") +
  labs(x='location', y='Relative abundance') + theme_bw() + scale_y_percent() + ggtitle("Surface samples (Normalized by sum, N = 26)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))

plot_bar(normalized_surf_abun, x="Phylum", fill = "Phylum") + ggtitle("Surface samples (N = 26)") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack", na.rm = TRUE) +
  labs(x='Phyla', y='Relative abundance') 

plot_bar(normalized_surf_abun, x="Class", fill = "Class") + ggtitle("Surface samples (N = 26)") +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  labs(x='Class', y='Relative abundance') 

(2) subsurface sediment
asv_mat_sub<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_most abund_sub.xlsx", sheet = "ASV_table"))
tax_mat_sub<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_most abund_sub.xlsx", sheet = "taxa_table"))
sample_df_sub<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_most abund_sub.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_mat_sub)<-asv_mat_sub$asv
asv_mat_sub<-asv_mat_sub%>% select(-asv)

row.names(tax_mat_sub)<-tax_mat_sub$asv
tax_mat_sub<-tax_mat_sub%>% select(-asv)

row.names(sample_df_sub)<-sample_df_sub$sample
sample_df_sub<-sample_df_sub%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_mat_sub<-as.matrix(asv_mat_sub)
tax_mat_sub<-as.matrix(tax_mat_sub)

## Create phyloseq objects.
(A)
ASV_sub = otu_table(asv_mat_sub, taxa_are_rows = TRUE)
TAX_sub  = tax_table(tax_mat_sub)
samples_sub = sample_data(sample_df_sub)
metabarcoding_sub_abun <-phyloseq(ASV_sub, TAX_sub, samples_sub)
metabarcoding_sub_abun
(B)
data_sub_abun <- phyloseq(otu_table(asv_mat_sub, taxa_are_rows = TRUE), sample_data(sample_df_sub))
tax_table(data_sub_abun) <- as.matrix(tax_mat_sub)
data_sub_abun

## Normalize counts by the sample sums (Fractions).
otuS <- data.frame(otu_table(data_sub_abun))
sumS <- apply(otuS, 2, sum)
norm_otuS <- otuS
for (i in 1:ncol(otuS)){
  norm_otuS[,i] <- otuS[,i]/sumS[i]
}
normalized_sub_abun <- phyloseq(otu_table(norm_otuS, taxa_are_rows = TRUE), tax_table(as.matrix(tax_mat_sub)), sample_data(sample_df_sub))
normalized_sub_abun

plot_bar(normalized_sub_abun, x = 'location', fill = "Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity", position = "stack") +
  labs(x='location', y='Relative abundance') + theme_bw() + scale_y_percent() + ggtitle("Sub-surface samples (Normalized by sum, N = 23)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))

plot_bar(normalized_sub_abun, x = 'location', fill = "Class") + geom_bar(aes(color=Class, fill=Class), stat = "identity", position = "stack") +
  labs(x='location', y='Relative abundance') + theme_bw() + scale_y_percent() + ggtitle("Sub-surface samples (Normalized by sum, N = 23)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))

plot_bar(normalized_sub_abun, x = 'location', fill = "Genus") + geom_bar(aes(color=palette.colors("alphabet"), fill=Genus), stat = "identity", position = "stack") +
  labs(x='location', y='Relative abundance') + theme_bw() + scale_y_percent() + ggtitle("Sub-surface samples (Normalized by sum, N = 23))") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))

plot_bar(normalized_sub_abun, x="Phylum", fill = "Phylum") + ggtitle("Sub-surface samples (N = 23)") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack", na.rm = TRUE) +
  labs(x='Phyla', y='Relative abundance') 

plot_bar(normalized_sub_abun, x="Class", fill = "Class") + ggtitle("Sub-surface samples (N = 23)") +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  labs(x='Class', y='Relative abundance') 



# Install the released version from CRAN:
install.packages("pals")

# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("kwstat/pals")

install.packages("tidyverse")
install.packages("dplyr")
install.packages("hrbrthemes")
install.packages("gcookbook")
library(microbiome)
library(dplyr)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)

########################################################################################################################################################################################################################################################
                                               ## Visualize the relative abundance for taxa higher than 500 ##
# create phyloseq file (short file after discard the taxa lower than 100). three tables are required (ASV table, Taxonomy table and Sample info table - these can be all in single excel file on three separate sheets) #####
(1) surface sediment
asv_mat_surf<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_surface176.xlsx", sheet = "ASV_table"))
tax_mat_surf<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_surface176.xlsx", sheet = "taxa_table"))
sample_df_surf<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB/phyloseq_surface176.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_mat_surf)<-asv_mat_surf$asv
asv_mat_surf<-asv_mat_surf%>% select(-asv)

row.names(tax_mat_surf)<-tax_mat_surf$asv
tax_mat_surf<-tax_mat_surf%>% select(-asv)

row.names(sample_df_surf)<-sample_df_surf$sample
sample_df_surf<-sample_df_surf%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_mat_surf<-as.matrix(asv_mat_surf)
tax_mat_surf<-as.matrix(tax_mat_surf)

## Create phyloseq objects. 
(A)
ASV_surf = otu_table(asv_mat_surf, taxa_are_rows = TRUE)
TAX_surf  = tax_table(tax_mat_surf)
samples_surf = sample_data(sample_df_surf)
metabarcoding_surf <-phyloseq(ASV_surf, TAX_surf, samples_surf)
metabarcoding_surf
(B)
data_surf <- phyloseq(otu_table(asv_mat_surf, taxa_are_rows = TRUE), sample_data(sample_df_surf))
tax_table(data_surf) <- as.matrix(tax_mat_surf)
data_surf

(2) subsurface sediment
asv_mat_sub<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/decipher/phyloseq_table109_subsurface.xlsx", sheet = "ASV_table"))
tax_mat_sub<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/decipher/phyloseq_table109_subsurface.xlsx", sheet = "taxa_table"))
sample_df_sub<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/decipher/phyloseq_table109_subsurface.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_mat_sub)<-asv_mat_sub$asv
asv_mat_sub<-asv_mat_sub%>% select(-asv)

row.names(tax_mat_sub)<-tax_mat_sub$asv
tax_mat_sub<-tax_mat_sub%>% select(-asv)

row.names(sample_df_sub)<-sample_df_sub$sample
sample_df_sub<-sample_df_sub%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_mat_sub<-as.matrix(asv_mat_sub)
tax_mat_sub<-as.matrix(tax_mat_sub)

## Create phyloseq objects.
(A)
ASV_sub = otu_table(asv_mat_sub, taxa_are_rows = TRUE)
TAX_sub  = tax_table(tax_mat_sub)
samples_sub = sample_data(sample_df_sub)
metabarcoding_sub <-phyloseq(ASV_sub, TAX_sub, samples_sub)
metabarcoding_sub
(B)
data_sub <- phyloseq(otu_table(asv_mat_sub, taxa_are_rows = TRUE), sample_data(sample_df_sub))
tax_table(data_sub) <- as.matrix(tax_mat_sub)
data_sub

## data preprocessing (sums, maximum, minimum)
sample_sums(data_sub)
max(sample_sums(data_sub))
min(sample_sums(data_sub))

## Normalize counts by the sample sums (Fractions).
otus_surf <- data.frame(otu_table(data_surf))
sums_surf <- apply(otus_surf, 2, sum)
norm_otus_surf <- otus_surf
for (i in 1:ncol(otus_surf)){
  norm_otus_surf[,i] <- otus_surf[,i]/sums_surf[i]
}
normalized_data_surf <- phyloseq(otu_table(norm_otus_surf, taxa_are_rows = TRUE), tax_table(as.matrix(tax_mat_surf)), sample_data(sample_df_surf))

## Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(data_surf))
standf = function(x, t=total) round(t * (x / sum(x)))
normalized_surf = transform_sample_counts(data_surf, standf)
normalized_surf

## Community composition 
# plot bar plot graphs by phyla across samples and processing step (making sure order of samples is correct across x axis) #####
cols <- c("#01B3ED", "#FFF8B8", "#B34E7E")
pal<-colorRampPalette(cols)

a<-brewer.pal(9, "Blues")
b<-brewer.pal(9, "Oranges")
c <- brewer.pal(9, "Greens")
d <- brewer.pal(9, "Purples")
e <- brewer.pal(12, "Set3")
cols<-c(a[c(9,7,6,5,3)],b[c(9,7,6,5,3)],c[c(9,7,6,5,3)],d[c(9,7,6,5,3)],e)

F = plot_bar(normalized_data_sub, fill="phylum")
F

F1 = plot_bar(normalized_surf, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity", position = "stack", levels = c("KY1.1F", "KY2.2F", "KY3.3F", "KW1.4F", "KW2.5F", "KW3.6F", "DJ1.7F", "DJ2.8F", "DJ3.9F", "HM1.10F", "HM2.11F", "FC1.12F", "FC2.13F", "FC3.14F", "HS1.15F", "HS2.16F", "HS3.17F", "HS4.18F", "OB1.19F", "OB2.20F", "OB3.21F", "YK1.22F", "YK2.23F", "MA1.24F", "MA2.25F", "MA3.26F")) + labs(x='Sites', y='Relative abundance') + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=270)) + ggtitle("Surface sediments")
F1

p2 = plot_bar(normalized_data, x = "layer", fill = "phylum") + 
  geom_bar(aes(color=phylum, fill=phylum), stat = "identity", position = "stack") +
  labs(x='Layers', y='Relative abundance') +
  theme(axis.text.x = element_text(size = 10, angle=0))
p2

p3 = plot_bar(normalized_data, x = 'area', fill = "phylum") +
  geom_bar(aes(color=phylum, fill=phylum), stat = "identity", position = "stack") +
  labs(x='Regions', y='Relative abundance') + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))
p3

p7 = plot_bar(normalized_data, x="phylum", fill = "phylum") +
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack", na.rm = TRUE) +
  labs(x='Phyla', y='Relative abundance')
p7

p8 = plot_bar(normalized_data, x="class", fill = "class") +
  geom_bar(aes(color=class, fill=class), stat="identity", position="stack") +
  labs(x='Class', y='Relative abundance')
p8
#######################################################################################################################################
## Heatmaps for most abundant OTUs >500 (use normalized reads by median sequencing depth not fraction).
(1) surface
## Normalize number of reads in each sample using median sequencing depth.
totalmap = median(sample_sums(metabarcoding_surf_abun))  # find median sample read count (98396).
standfmap = function(x, t=totalmap) round(t * (x / sum(x)))  # function to standardize to median sample read count.
normalized_medianmap = transform_sample_counts(metabarcoding_surf_abun, standfmap)  # apply to phyloseq object.
normalized_medianmap
ntaxa(normalized_medianmap)
sample_sums(normalized_medianmap)

metabarcoding_sub_abun

metabarcoding_surf_abun

plot_heatmap(normalized_medianmap, method = "NMDS", distance = "bray", low="#000033", high="#CCFF66")
plot_heatmap(normalized_medianmap, "NMDS", "bray", "location", low="#000033", high="#CCFF66", na.value="white")
plot_heatmap(normalized_medianmap, method = "NMDS", distance = phyloseq::distance(metabarcoding_surf_abun, method = "bray", binary = TRUE),
             taxa.label = "Class", taxa.order = "Class", trans=NULL, low="#000033", high="#CCFF66", na.value="white")

(1) sub-surface
## Normalize number of reads in each sample using median sequencing depth.
totalmap = median(sample_sums(metabarcoding_sub_abun))  # find median sample read count (73272).
standfmap = function(x, t=totalmap) round(t * (x / sum(x)))  # function to standardize to median sample read count.
normalized_medianmap = transform_sample_counts(metabarcoding_sub_abun, standfmap)  # apply to phyloseq object.
normalized_medianmap
ntaxa(normalized_medianmap)
sample_sums(normalized_medianmap)
plot_heatmap(normalized_medianmap, "NMDS", "bray", "location", low="#000033", high="#CCFF66", na.value="white")
plot_heatmap(normalized_medianmap, method = "NMDS", distance = phyloseq::distance(metabarcoding_sub_abun, method = "bray", binary = TRUE),
             taxa.label = "Class", trans=NULL, low="#000033", high="#CCFF66", na.value="white")
















































asv_matS<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB//phyloseq_surface 226.xlsx", sheet = "ASV_table"))
tax_matS<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB//phyloseq_surface 226.xlsx", sheet = "taxa_table"))
sample_dfS<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/GTDB//phyloseq_surface 226.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_matS)<-asv_matS$asv
asv_matS<-asv_matS%>% select(-asv)

row.names(tax_matS)<-tax_matS$asv
tax_matS<-tax_matS%>% select(-asv)

row.names(sample_dfS)<-sample_dfS$sample
sample_dfS<-sample_dfS%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_matS<-as.matrix(asv_matS)
tax_matS<-as.matrix(tax_matS)

## Create phyloseq objects.
(A)
ASVS = otu_table(asv_matS, taxa_are_rows = TRUE)
TAXS  = tax_table(tax_matS)
sampleS = sample_data(sample_dfS)
metabarcodingS <-phyloseq(ASVS, TAXS, sampleS)
metabarcodingS

## Normalize counts by the sample sums (Fractions).
otusS <- data.frame(otu_table(metabarcodingS))
sumsS <- apply(otusS, 2, sum)
norm_otus_S <- otusS
for (i in 1:ncol(otusS)){
  norm_otus_S[,i] <- otusS[,i]/sumsS[i]
}
normalized_data_ <- phyloseq(otu_table(norm_otus_surf, taxa_are_rows = TRUE), tax_table(as.matrix(tax_mat_surf)), sample_data(sample_df_surf))


## Normalize number of reads in each sample using median sequencing depth (The number of reads used for normalization is 65424).
total = median(sample_sums(metabarcodingS))
standf = function(x, t=total) round(t * (x / sum(x)))
normalized_median = transform_sample_counts(metabarcodingS, standf)
normalized_median

## A basic heatmap using the default parameters.
plot_heatmap(normalized_median, method = "NMDS", distance = "bray")

plot_heatmap(normalized_median, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "phylum", taxa.order = "phylum", 
             trans=NULL, low="beige", high="red", na.value="beige")


plot_heatmap(normalized_median, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "class", taxa.order = "class", 
             trans=NULL, low="beige", high="red", na.value="beige")

# Many different built-in distances can be used
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

plot_heatmap(normalized_median, method = "MDS", distance = "euclidean", 
             taxa.label = "phylum", taxa.order = "phylum", trans=NULL, low="beige", high="red", na.value="beige")
# longfile
plot_heatmap(normalized_medianL, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "class", taxa.order = "class", trans=NULL, low="beige", high="red", na.value="beige")
# 109 
plot_heatmap(normalized_median109, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "class", taxa.order = "class", trans=NULL, low="beige", high="red", na.value="beige")
plot_heatmap(normalized_median109, method = "NMDS", distance = "bray")
# 157
plot_heatmap(normalized_median157, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "class", taxa.order = "class", trans=NULL, low="beige", high="red", na.value="beige")
plot_heatmap(normalized_median157, method = "NMDS", distance = "bray")

########################################################################################################################################

## create phyloseq file (short file 109).
asv_mat109<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/decipher/phyloseq_table109.xlsx", sheet = "ASV_table"))
tax_mat109<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/decipher/phyloseq_table109.xlsx", sheet = "taxa_table"))
sample_df109<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/decipher/phyloseq_table109.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_mat109)<-asv_mat109$asv
asv_mat109<-asv_mat109%>% select(-asv)

row.names(tax_mat109)<-tax_mat109$asv
tax_mat109<-tax_mat109%>% select(-asv)

row.names(sample_df109)<-sample_df109$sample
sample_df109<-sample_df109%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_mat109<-as.matrix(asv_mat109)
tax_mat109<-as.matrix(tax_mat109)

## Create phyloseq objects.
(A)
ASV109 = otu_table(asv_mat109, taxa_are_rows = TRUE)
TAX109  = tax_table(tax_mat109)
samples109 = sample_data(sample_df109)
metabarcoding109 <-phyloseq(ASV109, TAX109, samples109)
metabarcoding109
(B)
data109 <- phyloseq(otu_table(asv_mat109, taxa_are_rows = TRUE), sample_data(sample_df109))
tax_table(data109) <- as.matrix(tax_mat109)
data109

## Normalize number of reads in each sample using median sequencing depth.
total109 = median(sample_sums(metabarcoding109))
standf109 = function(x, t=total109) round(t * (x / sum(x)))
normalized_median109 = transform_sample_counts(metabarcoding109, standf109)
normalized_median109

## create phyloseq file (short file).
asv_mat157<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/decipher/phyloseq_table157.xlsx", sheet = "ASV_table"))
tax_mat157<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/decipher/phyloseq_table157.xlsx", sheet = "taxa_table"))
sample_df157<-as.data.frame(read_excel("E:/sequencing data/microbiota/dada2 output/decipher/phyloseq_table157.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_mat157)<-asv_mat157$asv
asv_mat157<-asv_mat157%>% select(-asv)

row.names(tax_mat157)<-tax_mat157$asv
tax_mat157<-tax_mat157%>% select(-asv)

row.names(sample_df157)<-sample_df157$sample
sample_df157<-sample_df157%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_mat157<-as.matrix(asv_mat157)
tax_mat157<-as.matrix(tax_mat157)

## Create phyloseq objects.
(A)
ASV157 = otu_table(asv_mat157, taxa_are_rows = TRUE)
TAX157  = tax_table(tax_mat157)
samples157 = sample_data(sample_df157)
metabarcoding157 <-phyloseq(ASV157, TAX157, samples157)
metabarcoding157
(B)
data157 <- phyloseq(otu_table(asv_mat157, taxa_are_rows = TRUE), sample_data(sample_df157))
tax_table(data157) <- as.matrix(tax_mat157)
data157

## Normalize number of reads in each sample using median sequencing depth.
total157 = median(sample_sums(data157))
standf157 = function(x, t=total157) round(t * (x / sum(x)))
normalized_median157 = transform_sample_counts(data157, standf157)
normalized_median157

## Ordination
# 157
# Do multivariate analysis based on Bray-Curtis distance and NMDS ordination.
microbiota.ord <- ordinate(normalized_median157, "NMDS", "bray")
# Plot ASVs by breaking according to taxonomic phyla (most abundant file 31 taxa).
plot_ordination(normalized_median157, microbiota.ord, type="taxa", color="class", 
                title="ASVs", label="Class") + facet_wrap(~phylum, 3)
# Display samples and enlarge the points to make it more easy to read.
plot_ordination(normalized_median157, microbiota.ord, type="sample", color="layer", 
                shape="layer", title="Samples") + geom_point(size=3)
# Display both samples and ASVs but in 2 different panels.
plot_ordination(normalized_median157, microbiota.ord, type="split", color="class", 
                shape="layer", title="biplot", label = "station") + geom_point(size=3)

#109
microbiota.ord <- ordinate(normalized_median109, "NMDS", "bray")
# Plot ASVs by breaking according to taxonomic phyla (most abundant file 31 taxa).
plot_ordination(normalized_median109, microbiota.ord, type="taxa", color="class", 
                title="ASVs", label="Class") + facet_wrap(~phylum, 3)
# Display samples and enlarge the points to make it more easy to read.
plot_ordination(normalized_median109, microbiota.ord, type="sample", color="layer", 
                shape="layer", title="Samples") + geom_point(size=3)
# Display both samples and ASVs but in 2 different panels.
plot_ordination(normalized_median109, microbiota.ord, type="split", color="class", 
                shape="layer", title="biplot", label = "station") + geom_point(size=3)

#26549
microbiota.ord <- ordinate(normalized_medianL, "NMDS", "bray")
# Plot ASVs by breaking according to taxonomic phyla (most abundant file 31 taxa).
plot_ordination(normalized_medianL, microbiota.ord, type="taxa", color="class", 
                title="ASVs", label="Class") + facet_wrap(~phylum, 3)
# Display samples and enlarge the points to make it more easy to read.
plot_ordination(normalized_medianL, microbiota.ord, type="sample", color="layer", 
                shape="layer", title="Samples") + geom_point(size=3)
# Display both samples and ASVs but in 2 different panels.
plot_ordination(normalized_medianL, microbiota.ord, type="split", color="class", 
                shape="layer", title="biplot", label = "station") + geom_point(size=3)

# 31 taxa
microbiota.ord <- ordinate(normalized_median, "NMDS", "bray")
# Plot ASVs by breaking according to taxonomic phyla (most abundant file 31 taxa).
plot_ordination(normalized_median, microbiota.ord, type="taxa", color="class", 
                title="ASVs", label="Class") + facet_wrap(~phylum, 3)
# Display samples and enlarge the points to make it more easy to read.
plot_ordination(normalized_median, microbiota.ord, type="sample", color="layer", 
                shape="layer", title="Samples") + geom_point(size=3)
# Display both samples and ASVs but in 2 different panels.
plot_ordination(normalized_median, microbiota.ord, type="split", color="class", 
                shape="layer", title="biplot", label = "station") + geom_point(size=3)

# Network analysis
# 157
plot_net(normalized_median157, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="class", point_label="phylum")
plot_net(normalized_median157, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="class", point_label="class") 

# 109
plot_net(normalized_median109, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="class", point_label="phylum")
plot_net(normalized_median109, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="class", point_label="class")

# 31
plot_net(normalized_median, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="class", point_label="phylum")
plot_net(normalized_median, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="class", point_label="class")












######## DESeq2/vegan packages for beta-diversity analysis ##########
install.packages("DESeq2")
library(vegan)
library(DESeq2)
library(dendextend)

##### three separate files are required for this pipeline, an ASV count table, a sample info table and a taxonomy table (these can be made in Excel from the insect output ASV table, similar to the phyloseq pipeline requirements) #####
count <- read.table("harwich_count.txt", header = TRUE, row.names = 1, check.names = F, sep = "\t")
sample_info <- read.table("harwich_sample_info.txt", header = TRUE, row.names = 1, check.names = F, sep = "\t")
tax <- as.matrix(read.table("harwich_tax.tsv", header = TRUE, row.names = 1, check.names = F, sep = "\t"))

count <- asv_mat
tax <- tax_mat
sample_info <- sample_df

head(count)
head(sample_info)
head(tax)

sample_info$color <- rainbow(49)
rarecurve(t(count), col = sample_info$color, step = 100 , lwd=3, ylab = "ASVs", label = F)
abline(v=(min(rowSums(t(count)))), col = "black")
title("Rarefraction of Shiangshan Samples")



















# communty composition 
tax <- data.frame(tax_table(normalized_data))
otuN <- data.frame(otu_table(normalized_data))
long_dataN <- cbind(tax, otuN)
long_dataN <- reshape2::melt(long_dataN)
#for every sample in long_dataN, add the sample description
#use as.character to remove factors
long_dataN$Sample <- as.character(long_dataN$variable)
for (i in 1:length(long_dataN$Sample)) {
  sample_name <- long_dataN$variable[[i]]
  long_dataN$Sample[[i]] <- as.character(sample_df$location[rownames(sample_df) == sample_name])
}
p=plot_bar(normalized_data, aes(x=variable, y=value, fill = "phylum")) +
  geom_bar(position='stack', stat = 'identity') +
  labs(x='Sample', y='Relative abundance') + theme_minimal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=270))