## this script oversees the full bioinformatic pipeline required to complete the analysis of Illumina MiSeq sequence data of COI amplicons (mlCOInintF/jgh20198 313bp amplicon region)
## this pipeline uses DADA2 for read processing/filtering to produce Amplicon Sequence Variants (asv), the insect package to assign taxonomy and LULU (+VSEARCH) for asv table curation and phyloseq & vegan packages for asv table visualisation and statistical tests 
## the following online tutorials were used to create the following script
# 1) DADA2 online tutorial (v 1.12) for read processing/filtering: http://benjjneb.github.io/dada2/tutorial.html
# 2) INSECT R package online tutorial for ASV taxonomy assignement: https://cran.r-project.org/web/packages/insect/vignettes/insect-vignette.html
# 3) LULU curation using vsearch option: https://github.com/tobiasgf/lulu


################ DADA2 pipeline ###################


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18")


install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") 
library(dada2)

BiocManager::install("DECIPHER")
BiocManager::install("phangorn")
BiocManager::install("phyloseq")
BiocManager::install("ggplot2")
library(DECIPHER)
library(phangorn)
library(phyloseq)
library(ggplot2)


file_path_COI <- "D:/E/sequencing data/Invertebrates (2 seasons in water)/raw data"
list.files(file_path_COI)
file.exists(file_path_COI)
fwdsF<-sort(list.files(file_path_COI, pattern="R1.clean.fastq", full.names = TRUE))
revsR<-sort(list.files(file_path_COI, pattern="R2.clean.fastq", full.names = TRUE))
fwdsF
revsR

## Inspect read quality profiles & visualize it.
#visualise forward read quality profiles
plotQualityProfile(fwdsF[1:40])
plotQualityProfile(fwdsF[1])
#visualise reverse read quality profiles
plotQualityProfile(revsR[1:40])
plotQualityProfile(revsR[1])


sample.names <- sapply(strsplit(basename(fwdsF), "_"), '[', 1)
sample.names
sample.names <- sapply(strsplit(basename(revsR), "_"), '[', 1)
sample.names
sample.names <- sapply(strsplit(basename(revsR), "R"), '[', 1)
sample.names

fnF <- file.path(file_path_COI, fwdsF)
fnR <- file.path(file_path_COI, revsR)
filtFs_COI <-file.path(file_path_COI, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs_COI <-file.path(file_path_COI, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs_COI) <- sample.names
names(filtRs_COI) <- sample.names
filtFs_COI
filtRs_COI


##filter and trim
#remove primers (25bp long from beginning of each sequence, remove phix control sequence, trim reads to 150 for both forward and reverse files as quality drops below a score of 30 around 150bp)
outS<- filterAndTrim(fwdsF, filtFs_COI, revsR,filtRs_COI, maxN=0, maxEE = c(2,2), truncQ=2, rm.phix = TRUE, compress = TRUE, multithread = T, trimLeft = 63, trimRight = 50)
head(outS)
dim(outS)
outS

## Now let's take a look at our filtered reads (filtFs_COI & filtRs_COI):
plotQualityProfile(filtFs_COI[1])
plotQualityProfile(filtFs_COI[1:40])
plotQualityProfile(filtRs_COI[1])
plotQualityProfile(filtRs_COI[1:40])

#learn error rates, using a machine learning algorhithm
errF<- learnErrors(filtFs_COI, multithread = TRUE)
errR<- learnErrors(filtRs_COI, multithread = TRUE)
#visualise estimated error rates with
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
gc()

#sample inference - where the core sample inference algorithm is applied
dadaFs_COI <- dada(filtFs_COI, err = errF, multithread = TRUE)
dadaRs_COI <- dada(filtRs_COI, err = errR, multithread = TRUE)
dadaFs_COI[[1]]
dadaRs_COI[[1]]

#merge paired reads aka align forward and reverse complement reads
mergers_COI <- mergePairs(dadaFs_COI, filtFs_COI, dadaRs_COI, filtRs_COI, verbose = TRUE)
head(mergers_COI[[1]])
mergers_COI[[1]]

#construct sequence table (construct ASV table, higher-resolution version than traditional OTU table)
seqtab_COI <- makeSequenceTable(mergers_COI)
dim(seqtab_COI)
#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_COI)))

#remove sequences that are shorter than 200bp or longer than 300bp as they are likely results of non-specific priming
seqtab2_COI <- seqtab_COI[,nchar(colnames(seqtab_COI)) %in% 200:350]
dim(seqtab2_COI)

#remove chimeras
seqtab.nochim_COI <- removeBimeraDenovo(seqtab_COI, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab2.nochim_COI <- removeBimeraDenovo(seqtab2_COI, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim_COI)
dim(seqtab2.nochim_COI)
sum(seqtab2.nochim_COI)/sum(seqtab2_COI)
sum(seqtab.nochim_COI)/sum(seqtab_COI)
table(nchar(getSequences(seqtab2.nochim_COI)))
hist(nchar(getSequences(seqtab2.nochim_COI)), main="Distribution of sequence lengths")


# track reads throughout the above pipeline and print out number of reads that make it through each step
getN_COI <- function(x) sum(getUniques(x))
track_COI <- cbind(outS, sapply(dadaFs_COI, getN_COI), sapply(dadaRs_COI, getN_COI), sapply(mergers_COI, getN_COI), rowSums(seqtab2.nochim_COI))
colnames(track_COI) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_COI) <- sample.names
head(track_COI)
track_COI
write.csv(track_COI, "D:/E/sequencing data/Invertebrates (2 seasons in water)/DADA output new.csv")

summary_tab <- data.frame(dada2_input=outS[,1],
                          filtered=outS[,2], denoisedF=sapply(dadaFs_COI, getN_COI),
                          denoisedR=sapply(dadaRs_COI, getN_COI), merged=sapply(mergers_COI, getN_COI),
                          nonchim=rowSums(seqtab2.nochim_COI),
                          final_perc_reads_retained=round(rowSums(seqtab2.nochim_COI)/outS[,1]*100, 1))
write.csv(summary_tab, "D:/E/sequencing data/Invertebrates (2 seasons in water)/DADA output new/summary_trim.csv")


# make output fasta file with sequences and associated ASV ID strings
library(ShortRead)
sq <- getSequences(seqtab2.nochim_COI)
writeFasta(sq, file="D:/E/sequencing data/Invertebrates (2 seasons in water)/DADA output new.fasta")

#giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqS <- colnames(seqtab2.nochim_COI)
asv_headerS <- vector(dim(seqtab2.nochim_COI)[2], mode="character")
for (i in 1:dim(seqtab2.nochim_COI)[2]) {
  asv_headerS[i] <- paste(">ASV", i, sep="_")
}

#making and writing out a fasta of our final ASV seqs:
asv_fastaS <- c(rbind(asv_headerS, asv_seqS))
write(asv_fastaS, "ASVs.fa")
write(asv_fastaS, "ASVs.csv")

#count table:
asv_tabS <- t(seqtab2.nochim_COI)
row.names(asv_tabS) <- sub(">", "", asv_headerS)
write.table(asv_tabS, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tabS, "ASVs_counts.csv", quote=F, col.names=NA)


                 ########## Assigning Taxonomy  of output ASVs from the DADA2 pipeline, using different databases ############

## (3) Genome Taxonomy Database (GTDB): The new version of our dada2
taxa <- assignTaxonomy(seqtab2.nochim_16rRNA, "D:/E/sequencing data/GTDB database/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz", multithread=TRUE)
Genus <- assignTaxonomy(taxa, "D:/E/sequencing data/GTDB database/GTDB_bac120_arc122_ssu_r202_Genus.fa.gz")
taxa.print <- taxa 
rownames(taxa.print) <- NULL  # Remove sequence rownames for display only
head(taxa.print)
write.csv(taxa.print, "D:/Dr. walaa/New folder/Taxonomy.GTDB.csv")


## tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_taxS <- t(sapply(taxa, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_taxS) <- ranks
rownames(asv_taxS) <- gsub(pattern=">", replacement="", x=asv_headerS)
write.table(asv_taxS, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

#write.table(longDF and shortDF)
longDF <- cbind(taxa.print, t(seqtab2.nochim_16rRNA))
write.csv(longDF, "D:/Dr. walaa/New folder/16srRNA_asv_table_long", sep = "\t")
write.csv(longDF, "D:/Dr. walaa/New folder/16srRNA_asv_table_long.csv")

longDF <- read.csv("D:/Dr. walaa/New folder/16srRNA_asv_table_long.csv",header=T)
taxa <- aggregate(longDF[2:8], longDF["Class"], head, 1)
taxa <- taxa[-1]
counts <- aggregate(longDF[9:ncol(longDF)], longDF["Class"], sum)
shortDF <- merge(taxa, counts, by = "Class")
write.csv(shortDF, "D:/Dr. walaa/New folder/16srRNA_asv_table_short", sep = "\t")
write.csv(shortDF, "D:/Dr. walaa/New folder/16srRNA_asv_table_short.csv")
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
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
taxonomyCOI <- read.csv("D:/E/sequencing data/invertebrates paper/dada outputs/longfile siangshan.csv")
HowManyKingdom <- length(unique(unlist(taxonomyCOI[,c("kingdom")]))) #1
HowManyPhyla <- length(unique(unlist(taxonomyCOI[,c("phylum")]))) #46
HowManyClass <- length(unique(unlist(taxonomyCOI[,c("class")]))) #133
HowManyOrder <- length(unique(unlist(taxonomyCOI[,c("order")]))) #450
HowManyFamily <- length(unique(unlist(taxonomyCOI[,c("family")]))) #1278
HowManyGenus <- length(unique(unlist(taxonomyCOI[,c("genus")]))) #2168
HowManySpecies <- length(unique(unlist(taxonomyCOI[,c("species")]))) #2689


################## create phyloseq file (long file 5,475) for alpha diversity measure######################

asv_mat<-as.data.frame(read_excel("D:/E/sequencing data/invertebrates paper/dada outputs/phyloseq.xlsx", sheet = "ASV_table"))
tax_mat<-as.data.frame(read_excel("D:/E/sequencing data/invertebrates paper/dada outputs/phyloseq.xlsx", sheet = "taxa_table"))
sample_df<-as.data.frame(read_excel("D:/E/sequencing data/invertebrates paper/dada outputs/phyloseq.xlsx", sheet = "sample"))

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
ASV = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX  = tax_table(tax_mat)
samples = sample_data(sample_df)
metabarcoding <-phyloseq(ASV, TAX, samples)
metabarcoding


## data preprocessing (sums, maximum, minimum)
reads <- sample_sums(metabarcoding)  # number of reads/seqs. per sample
max(sample_sums(metabarcoding))
min(sample_sums(metabarcoding))
length(which(reads<5000))
reads
DJ1_S  DJ1_W  DJ2_S  DJ2_W  FC1_S  FC1_W  FC2_S  FC2_W  FC3_S  FC3_W  HM1_S  HM1_W  HM2_S  HM2_W  HS1_S  HS1_W  HS2_S 
213622  36775 233456  41282  93236  28489  97301  15851  59177  23793 109055  30477 111665  36912 151019  58730  71554 
HS2_W  HS3_S  HS3_W  HS4_S  HS4_W  KW1_S  KW1_W  KW2_S  KW2_W  KW3_S  KW3_W  KY1_S  KY1_W  KY2_S  KY2_W   MA_S   MA_W 
18866  63197  55835  62327  80443 108969  44657  73742  40253  87528  29059  68079  21852  71646  34317 109294  81839 
OB_S   OB_W  YK1_S  YK1_W  YK2_S  YK2_W 
58696  66980 105275 107803 207265  47561 

## Normalize number of reads in each sample using median sequencing depth (total = 65088.5).
total = median(sample_sums(metabarcoding))
stand = function(x, t=total) round(t * (x / sum(x)))
normalized_COI = transform_sample_counts(metabarcoding, stand)
normalized_COI
ntaxa(normalized_COI)
sample_sums(normalized_COI)
DJ1_S DJ1_W DJ2_S DJ2_W FC1_S FC1_W FC2_S FC2_W FC3_S FC3_W HM1_S HM1_W HM2_S HM2_W HS1_S HS1_W HS2_S HS2_W HS3_S HS3_W 
65092 65086 65094 65089 65072 65101 65069 65087 65060 65089 65090 65092 65072 65086 65076 65092 65083 65092 65064 65070 
HS4_S HS4_W KW1_S KW1_W KW2_S KW2_W KW3_S KW3_W KY1_S KY1_W KY2_S KY2_W  MA_S  MA_W  OB_S  OB_W YK1_S YK1_W YK2_S YK2_W 
65060 65099 65091 65080 65086 65085 65078 65089 65137 65098 65116 65090 65089 65091 65081 65111 65080 65073 65119 65091

## Plot Chao1 richness estimator and Shannon diversity estimator.
plot_richness(normalized_COI, measures=c("Chao1", "Shannon"))
sample <- plot_richness(normalized_COI, x="sample", color="season", measures=c("Chao1", "Shannon"), 
                        title = paste0("Standardized to total reads, N=", nsamples(normalized_COI)))
sample 
area <- plot_richness(normalized_COI, x="area", color="season", measures=c("Observed", "Chao1", "Shannon"), 
                      title = paste0("Standardized to total reads, N=", nsamples(normalized_COI)))
area 
season <- plot_richness(normalized_COI, x="season", color="season", measures=c("Observed", "Chao1", "Shannon"), 
                       title = paste0("Standardized to total reads, N=", nsamples(normalized_COI)))
season
area +  stat_compare_means(method = "t.test", label.y = 0.8)
kruskal.test(area~Observed, data = div) # p-value= 0.5557
kruskal.test(area~Chao1, data = div) # p-value= 0.5557
kruskal.test(area~Shannon, data = div) # p-value= 0.4713

est <- estimate_richness(normalized_COI, split = TRUE, measures = NULL)
tem <- cbind(sample_data(normalized_COI), est)
aggregate(tem$Chao1, list(tem$location), mean)  # calculate the means
write.csv(tem, file = "diversity_normalized_richness COI.csv")


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

div <- read.csv("D:/E/sequencing data/invertebrates paper/dada outputs/diversity_normalized_richness COI.csv")
div <- drop_na(div)
observed <- ggbetweenstats(div, x = season, y = Observed)
observed
observed <- observed + labs(x = "season", y = "observed_ASVs", title = paste0("Standardized to total reads, N=", nsamples(normalized_COI))) +  
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

chao1 <- ggbetweenstats(div, x = season, y = Chao1)
chao1
chao1 <- chao1 + labs(x = "season", y = "Chao1_Richness", title = paste0("Standardized to total reads, N=", nsamples(normalized_COI))) +  
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

shannon <- ggbetweenstats(div, x = season, y = Shannon)
shannon <- shannon + labs(x = "season", y = "Shannon_estimates", title = paste0("Standardized to total reads, N=", nsamples(normalized_COI))) +  
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
boxplot3 <- ggplot(div, aes(x=season, y=Observed)) + geom_boxplot() + geom_jitter(aes(color=season)) +xlab("") + ylab("observed_ASVs") + ggtitle("Standardized to total reads, N=40")
boxplot3
boxplot3 + stat_compare_means(method = "t.test", label.y = 0.8)


# plot ordination graph (Beta-diversity, NMDS)
set.seed(20)
COI_ord <- ordinate(normalized_COI, "NMDS", "bray", k=2, trymax=100)
COI_ord

Call:
  metaMDS(comm = veganifyOTU(physeq), distance = distance, k = 2,      trymax = 100) 
global Multidimensional Scaling using monoMDS
Data:     wisconsin(sqrt(veganifyOTU(physeq))) 
Distance: bray 
Dimensions: 2 
Stress:     0.1277176 
Stress type 1, weak ties
Best solution was repeated 3 times in 20 tries
The best solution was from try 8 (random start)
Scaling: centring, PC rotation, halfchange scaling 
Species: expanded scores based on ‘wisconsin(sqrt(veganifyOTU(physeq)))’  

stressplot(COI_ord)
legend("top", "stress = 0.127", bty = "n", cex = 1)
goodness(COI_ord)
[1] 0.020308004 0.021559250 0.017307778 0.018921383 0.011876912 0.012237106 0.013660214 0.017976616 0.016960063
[10] 0.012880048 0.023605864 0.021869325 0.015257994 0.013368414 0.035348867 0.015746787 0.014051332 0.012179246
[19] 0.013104766 0.009504751 0.018140449 0.012045057 0.026313068 0.015331323 0.014410784 0.010987397 0.037197690
[28] 0.025097936 0.026093002 0.022463193 0.023690121 0.017831311 0.021118975 0.029800029 0.021829317 0.024301864
[37] 0.022837735 0.019950341 0.023880670 0.016471559

color = c("season")
shape = c("season")
title = c("NMDS of COI, Bray-Curtis distance, K=2")
Ordination <- plot_ordination(normalized_COI, COI_ord, type="area", color=color, 
                              shape=shape, title=title) 
ordination = Ordination + theme(axis.text = element_text(size = 16, face = "bold"), axis.title = element_text(size = 18, face = "bold"), 
             legend.title = element_text(size = 14)) + theme_bw() + labs(color = color, shape = shape) + geom_point(size=2) + 
             geom_text(aes(label=rownames(sample_df), hjust=0.4, vjust=-0.5)) 
ordination

#ANOSIM
presabs = get_variable(normalized_COI)
presabs
region_ano = anosim(phyloseq::distance(normalized_COI, "bray"), presabs$season)
region_ano
summary(region_ano)

Call:
  anosim(x = phyloseq::distance(normalized_COI, "bray"), grouping = presabs$season) 
Dissimilarity: bray 
ANOSIM statistic R: 0.5975 
Significance: 0.001 
Permutation: free
Number of permutations: 999
Upper quantiles of permutations (null model):
  90%    95%  97.5%    99% 
  0.0374 0.0545 0.0732 0.0983 
Dissimilarity ranks between and within classes:
  0%    25%   50%    75% 100%   N
Between  8 379.25 536.0 657.25  780 400
Spring  41 258.25 354.0 515.00  764 190
Winter   1  61.50 124.5 192.75  530 190


##ANOSIM analysis
pc = read.csv("D:/E/sequencing data/invertebrates paper/dada outputs/ASV_counts.csv", header= TRUE)
com = pc[,1:ncol(pc)]
m_com = as.matrix(com)
ano = anosim(m_com, pc$season, distance = "bray", permutations = 9999)
ano


####### Rarefaction curve################
install.packages("permute")
install.packages("lattice")
library("vegan")
library("permute")
library("lattice")

# import the count table
Data=read.table("D:/E/sequencing data/invertebrates paper/dada outputs/ASV_counts.csv")
Data_sub=as.data.frame(Data[,c(1:40)])
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
title("Rarefaction curve")
(1)
Data=read.table("D:/Dr. walaa/New folder/ASVs_counts.csv")
Data_sub=as.data.frame(Data[,c(1:8)])
Data_t=(t(Data_sub))
raremax <- min(rowSums(Data_t))
raremax # 29375
col <- rainbow(8)
lty <- c("solid")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)
rare <- with(pars[1:8, ], rarecurve(Data_t, step = 20, sample = raremax, col = col, lty = lty, label = TRUE))
col <- rainbow(8)
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
title("Rarefaction curve")
(2)
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)
rarefaction <- with(pars[1:8,], rarecurve(Data_t, step = 20, sample = raremax, col = col, lty = lty, label = FALSE))
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
rarefied_COI <- rarefy_even_depth(metabarcoding)
rarecurve(t(otu_table(normalized_COI)), step = 8, cex = .6, col = "red")
rarecurve(metabarcoding, step = 20, col = "blue", cex = 0.6, label = FALSE)


#######################################################################################################################################################################################################################################################
## Visualize the relative abundance for taxa > .1% ##
# create phyloseq file (short file after discard the taxa lower than 0.1%). three tables are required (ASV table, Taxonomy table and Sample info table - these can be all in single excel file on three separate sheets) #####

asv_mat_COI<-as.data.frame(read_excel("D:/E/sequencing data/invertebrates paper/dada outputs/RA/phyloseq short invertebrates.xlsx", sheet = "ASV_table"))
tax_mat_COI<-as.data.frame(read_excel("D:/E/sequencing data/invertebrates paper/dada outputs/RA/phyloseq short invertebrates.xlsx", sheet = "taxa_table"))
sample_df_COI<-as.data.frame(read_excel("D:/E/sequencing data/invertebrates paper/dada outputs/RA/phyloseq short invertebrates.xlsx", sheet = "sample"))

## Transform into matrices.
row.names(asv_mat_COI)<-asv_mat_COI$asv
asv_mat_COI<-asv_mat_COI%>% select(-asv)

row.names(tax_mat_COI)<-tax_mat_COI$asv
tax_mat_COI<-tax_mat_COI%>% select(-asv)

row.names(sample_df_COI)<-sample_df_COI$sample
sample_df_COI<-sample_df_COI%>% select(-sample)

# Transform into matrices OTU and Tax tables (sample table can be left as data frame).
asv_mat_COI<-as.matrix(asv_mat_COI)
tax_mat_COI<-as.matrix(tax_mat_COI)

## Create phyloseq objects. 
(A)
ASV_COI = otu_table(asv_mat_COI, taxa_are_rows = TRUE)
TAX_COI  = tax_table(tax_mat_COI)
samples_COI = sample_data(sample_df_COI)
metabarcoding_COI <-phyloseq(ASV_COI, TAX_COI, samples_COI)
metabarcoding_COI
(B)
data_COI <- phyloseq(otu_table(asv_mat_COI, taxa_are_rows = TRUE), sample_data(sample_df_COI))
tax_table(data_COI) <- as.matrix(tax_mat_COI)
data_COI

## data preprocessing (sums, maximum, minimum)
sample_sums(data_COI)
max(sample_sums(data_COI))
min(sample_sums(data_COI))
DJ1_S  DJ1_W  DJ2_S  DJ2_W  FC1_S  FC1_W  FC2_S  FC2_W  FC3_S  FC3_W  HM1_S  HM1_W  HM2_S  HM2_W  HS1_S  HS1_W  HS2_S 
213622  36775 233456  41282  93236  28489  97301  15851  59177  23793 109055  30477 111665  36912 151019  58730  71554 
HS2_W  HS3_S  HS3_W  HS4_S  HS4_W  KW1_S  KW1_W  KW2_S  KW2_W  KW3_S  KW3_W  KY1_S  KY1_W  KY2_S  KY2_W   MA_S   MA_W 
18866  63197  55835  62327  80443 108969  44657  73742  40253  87528  29059  68079  21852  71646  34317 109294  81839 
OB_S   OB_W  YK1_S  YK1_W  YK2_S  YK2_W 
58696  66980 105275 107803 207265  47561 

[1] 233456 (Max)
[1] 15851 (Min)


## Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(data_COI))
standCOI = function(x, t=total) round(t * (x / sum(x)))
normalized_COI = transform_sample_counts(data_COI, standCOI)
normalized_COI

## Normalize counts by the sample sums (Fractions).
otu <- data.frame(otu_table(data_COI))
sum <- apply(otu, 2, sum)
norm_otu <- otu
for (i in 1:ncol(otu)){
  norm_otu[,i] <- otu[,i]/sum[i]
}
normalized_fraction <- phyloseq(otu_table(norm_otu, taxa_are_rows = TRUE), tax_table(as.matrix(tax_mat_COI)), sample_data(sample_df_COI))
normalized_fraction

plot_bar(normalized_fraction, x = 'season', fill = "phylum") + geom_bar(aes(color=phylum, fill=phylum), stat = "identity", position = "stack") +
  labs(x='season', y='Relative abundance') + theme_bw() + scale_y_percent() +ggtitle("COI distribution, N = 40") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))

plot_bar(normalized_surf, x = 'station', fill = "Class") + geom_bar(aes(color=Class, fill=Class), stat = "identity", position = "stack") +
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



library(hrbrthemes)
library(gcookbook)
library(tidyverse)


#https://joey711.github.io/phyloseq-demo/Restroom-Biogeography
wl_tr = transform_sample_counts(data_COI, function(x) 100 * x/sum(x))
preod=c("Bacillariophyta","Arthropoda","Mollusca","Cnidaria","Chlorophyta","Annelida","Chordata","Ascomycota","Nemertea","Rhodophyta","Porifera","Echinodermata","Platyhelminthes","Gastrotricha","Nematoda","Others","NA") 
col <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "deeppink4", "darksalmon", "yellow4","royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "orange4", "lightblue4", "maroon2", "tan", "azure3", "lightsalmon", "darkgrey", "gold", "beige", "azure1", "slategray1", "yellow", "plum", "pink", "green", "black", "lavender", "blue")

##### Phyla level #####
plot_bar(wl_tr, x = "season", fill=factor("phylum", preod)) +
  geom_bar(aes(fill = factor(phylum, preod)), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4","darksalmon","gold","pink","green"))+
  theme_bw() + ylab("Reletive abundance (%)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=0))

plot_bar(wl_tr, x = "area", fill=factor("phylum", preod)) +
  geom_bar(aes(fill = factor(phylum, preod)), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4","darksalmon","gold","pink","green"))+
  theme_bw() + ylab("Reletive abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=0))

##color palette (https://bookdown.org/hneth/ds4psy/D.3-apx-colors-basics.html)
preod=c("Bacillariophyceae","Insecta","Gastropoda","Hydrozoa","Scyphozoa","Mamiellophyceae","Xanthophyceae","Arachnida","Polychaeta","Aves","Lepidosauria","Sordariomycetes","Actinopteri","Enopla","Mediophyceae","Coscinodiscophyceae","Bangiophyceae","Malacostraca","Demospongiae","Ophiuroidea","Trematoda","Copepoda","Anthozoa","Amphibia","Cestoda","Chaetonotida","Ascidiacea","Chromadorea","Mammalia","Collembola","Asteroidea","Clitellata","Dothideomycetes","Others","NA") 
col <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "deeppink4", "darksalmon", "yellow4","royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "orange4", "lightblue4", "maroon2", "tan", "azure3", "lightsalmon", "darkgrey", "gold", "beige", "azure1", "slategray1", "yellow", "plum", "pink", "green", "black", "lavender", "blue")

##### class level #####
plot_bar(wl_tr, x = "season", fill=factor("class", preod)) +
  geom_bar(aes(fill = factor(class, preod)), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "deeppink4", "darksalmon", "yellow4","royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "orange4", "lightblue4", "maroon2", "tan", "azure3", "lightsalmon", "darkgrey", "gold", "beige", "azure1", "slategray1", "yellow", "plum", "pink", "green", "black", "lavender", "blue"))+
  theme_bw() + ylab("Reletive abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=0))

plot_bar(wl_tr, x = "area", fill=factor("class", preod)) +
  geom_bar(aes(fill = factor(class, preod)), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "deeppink4", "darksalmon", "yellow4","royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "orange4", "lightblue4", "maroon2", "tan", "azure3", "lightsalmon", "darkgrey", "gold", "beige", "azure1", "slategray1", "yellow", "plum", "pink", "green", "black", "lavender", "blue"))+
  theme_bw() + ylab("Reletive abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=0))


packageVersion("dada2")





plot_bar(data_COI, x="phylum", fill = "phylum") + ggtitle("Invertebrates (40)") + 
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack", na.rm = TRUE) +
  labs(x='Phyla', y='Abundance (Number of reads)')



