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


file_path_16srRNA <- "D:/E/sequencing data/Invertebrates (2 seasons in water)/raw data"
list.files(file_path_16srRNA)
file.exists(file_path_16srRNA)
fwdsF<-sort(list.files(file_path_16srRNA, pattern="R1.clean.fastq", full.names = TRUE))
revsR<-sort(list.files(file_path_16srRNA, pattern="R2.clean.fastq", full.names = TRUE))
fwdsF
revsR
#visualise forward read quality profiles
plotQualityProfile(fwdsF[1:8])
plotQualityProfile(fwdsF[1])
#visualise reverse read quality profiles
plotQualityProfile(revsR[1:8])
plotQualityProfile(revsR[1])


sample.names <- sapply(strsplit(basename(fwdsF), "_"), '[', 1)
sample.names
sample.names <- sapply(strsplit(basename(revsR), "_"), '[', 1)
sample.names
sample.names <- sapply(strsplit(basename(revsR), "R"), '[', 1)
sample.names

fnFsub <- file.path(file_path_16srRNA, fwdsF)
fnRsub <- file.path(file_path_16srRNA, revsR)
filtFs_16srRNA <-file.path(file_path_16srRNA, paste0(sample.names, "_F_filt.fastq"))
filtRs_16srRNA <-file.path(file_path_16srRNA, paste0(sample.names, "_R_filt.fastq"))
names(filtFs_16srRNA) <- sample.names
names(filtRs_16srRNA) <- sample.names


##filter and trim
#remove primers (25bp long from beginning of each sequence, remove phix control sequence, trim reads to 150 for both forward and reverse files as quality drops below a score of 30 around 150bp)
outS<- filterAndTrim(fwdsF, filtFs_16srRNA, revsR,filtRs_16srRNA,  maxN= 0, maxEE = c(2,2), truncQ = 2, truncLen = c(240, 230), rm.phix = TRUE, compress = TRUE, multithread = F)
head(outS)
dim(outS)
outS

#learn error rates, using a machine learning algorhithm
errF<- learnErrors(filtFs_16srRNA, multithread = TRUE)
errR<- learnErrors(filtRs_16srRNA, multithread = TRUE)
#visualise estimated error rates with
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
gc()

#sample inference - where the core sample inference algorithm is applied
dadaFs_16srRNA <- dada(filtFs_16srRNA, err = errF, multithread = TRUE)
dadaRs_16srRNA <- dada(filtRs_16srRNA, err = errR, multithread = TRUE)
dadaFs_16srRNA[[1]]
dadaRs_16srRNA[[1]]

#merge paired reads aka align forward and reverse complement reads
mergers_16srRNA <- mergePairs(dadaFs_16srRNA, filtFs_16srRNA, dadaRs_16srRNA, filtRs_16srRNA, verbose = TRUE)
head(mergers_16srRNA[[1]])
mergers_16srRNA[[1]]

#construct sequence table (construct ASV table, higher-resolution version than traditional OTU table)
seqtab_16srRNA <- makeSequenceTable(mergers_16srRNA)
dim(seqtab_16srRNA)
#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_16srRNA)))

#remove sequences that are shorter than 200bp or longer than 300bp as they are likely results of non-specific priming
seqtab2_16srRNA <- seqtab_16srRNA[,nchar(colnames(seqtab_16srRNA)) %in% 200:500]
dim(seqtab2_16srRNA)

#remove chimeras
seqtab.nochim_16srRNA <- removeBimeraDenovo(seqtab_16srRNA, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab2.nochim_16rRNA <- removeBimeraDenovo(seqtab2_16srRNA, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim_16srRNA)
dim(seqtab2.nochim_16rRNA)
sum(seqtab2.nochim_16rRNA)/sum(seqtab2_16srRNA)
sum(seqtab.nochim_16srRNA)/sum(seqtab_16srRNA)

# track reads throughout the above pipeline and print out number of reads that make it through each step
getN_16srRNA <- function(x) sum(getUniques(x))
track_16srRNA <- cbind(outS, sapply(dadaFs_16srRNA, getN_16srRNA), sapply(dadaRs_16srRNA, getN_16srRNA), sapply(mergers_16srRNA, getN_16srRNA), rowSums(seqtab2.nochim_16rRNA))
colnames(track_16srRNA) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_16srRNA) <- sample.names
head(track_16srRNA)
track_16srRNA
write.csv(track_16srRNA, "D:/Dr. walaa/New folder/UAE_16srRNA_trim_info.csv")

# make output fasta file with sequences and associated ASV ID strings
library(ShortRead)
sq <- getSequences(seqtab2.nochim_16rRNA)
writeFasta(sq, file="D:/Dr. walaa/New folder/UAE_16srRNA.fasta")


                 ########## Assigning Taxonomy  of output ASVs from the DADA2 pipeline, using different databases ############

## (3) Genome Taxonomy Database (GTDB): The new version of our dada2
taxa <- assignTaxonomy(seqtab2.nochim_16rRNA, "D:/E/sequencing data/GTDB database/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz", multithread=TRUE)
Genus <- assignTaxonomy(taxa, "D:/E/sequencing data/GTDB database/GTDB_bac120_arc122_ssu_r202_Genus.fa.gz")
taxa.print <- taxa 
rownames(taxa.print) <- NULL  # Remove sequence rownames for display only
head(taxa.print)
write.csv(taxa.print, "D:/Dr. walaa/New folder/Taxonomy.GTDB.csv")


#giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqS <- colnames(seqtab2.nochim_16rRNA)
asv_headerS <- vector(dim(seqtab2.nochim_16rRNA)[2], mode="character")
for (i in 1:dim(seqtab2.nochim_16rRNA)[2]) {
  asv_headerS[i] <- paste(">ASV", i, sep="_")
}

#making and writing out a fasta of our final ASV seqs:
asv_fastaS <- c(rbind(asv_headerS, asv_seqS))
write(asv_fastaS, "ASVs.fa")
write(asv_fastaS, "ASVs.csv")

#count table:
asv_tabS <- t(seqtab2.nochim_16rRNA)
row.names(asv_tabS) <- sub(">", "", asv_headerS)
write.table(asv_tabS, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tabS, "ASVs_counts.csv", quote=F, col.names=NA)

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
taxonomyGTDB <- read.csv("D:/Dr. walaa/New folder/16srRNA_asv_table_long.csv")
HowManyKingdom <- length(unique(unlist(taxonomyGTDB[,c("Kingdom")]))) # 4
HowManyPhyla <- length(unique(unlist(taxonomyGTDB[,c("Phylum")]))) # 39
HowManyClass <- length(unique(unlist(taxonomyGTDB[,c("Class")]))) # 80
HowManyOrder <- length(unique(unlist(taxonomyGTDB[,c("Order")]))) # 189
HowManyFamily <- length(unique(unlist(taxonomyGTDB[,c("Family")]))) # 343
HowManyGenus <- length(unique(unlist(taxonomyGTDB[,c("Genus")]))) # 718
HowManySpecies <- length(unique(unlist(taxonomyGTDB[,c("Species")]))) # 806


################## create phyloseq file (long file 4243) for alpha diversity measure######################

asv_mat<-as.data.frame(read_excel("D:/Dr. walaa/New folder/phyloseq.xlsx", sheet = "ASV_table"))
tax_mat<-as.data.frame(read_excel("D:/Dr. walaa/New folder/phyloseq.xlsx", sheet = "taxa_table"))
sample_df<-as.data.frame(read_excel("D:/Dr. walaa/New folder/phyloseq.xlsx", sheet = "sample"))

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
AH-E  AH-R  HP-E  HP-R  HS-E  HS-R  PT-E  PT-R 
41003 43447 38183 29742 47352 36984 44029 29375 

## Normalize number of reads in each sample using median sequencing depth (total = 39593).
total = median(sample_sums(metabarcoding))
stand = function(x, t=total) round(t * (x / sum(x)))
normalized_walaa = transform_sample_counts(metabarcoding, stand)
normalized_walaa
ntaxa(normalized_walaa)
sample_sums(normalized_walaa)
AH-E  AH-R  HP-E  HP-R  HS-E  HS-R  PT-E  PT-R 
39615 39601 39564 39614 39608 39520 39592 39568 

## Plot Chao1 richness estimator and Shannon diversity estimator.
plot_richness(normalized_walaa, measures=c("Chao1", "Shannon"))
sample <- plot_richness(normalized_walaa, x="sample", color="type", measures=c("Chao1", "Shannon"), 
                        title = paste0("Standardized to total reads, N=", nsamples(normalized_walaa)))
sample 
area <- plot_richness(normalized_walaa, x="station", color="type", measures=c("Observed", "Chao1", "Shannon"), 
                      title = paste0("Standardized to total reads, N=", nsamples(normalized_walaa)))
area 
layer <- plot_richness(normalized_walaa, x="layer", color="layer", measures=c("Observed", "Chao1", "Shannon"), 
                       title = paste0("Standardized to total reads, N=", nsamples(normalized_walaa)))
layer
area +  stat_compare_means(method = "t.test", label.y = 0.8)
kruskal.test(area~Observed, data = div) # p-value= 0.5557
kruskal.test(area~Chao1, data = div) # p-value= 0.5557
kruskal.test(area~Shannon, data = div) # p-value= 0.4713

est <- estimate_richness(normalized_walaa, split = TRUE, measures = NULL)
tem <- cbind(sample_data(normalized_walaa), est)
aggregate(tem$Chao1, list(tem$location), mean)  # calculate the means
write.csv(tem, file = "diversity_normalized_richness GTDB.csv")


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

div <- read.csv("D:/Dr. walaa/New folder/diversity_normalized_richness GTDB.csv")
div <- drop_na(div)
observed <- ggbetweenstats(div, x = type, y = Observed)
observed
observed <- observed + labs(x = "type", y = "observed_ASVs", title = paste0("Standardized to total reads, N=", nsamples(normalized_walaa))) +  
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

chao1 <- ggbetweenstats(div, x = type, y = Chao1)
chao1
chao1 <- chao1 + labs(x = "type", y = "Chao1_Richness", title = paste0("Standardized to total reads, N=", nsamples(normalized_walaa))) +  
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

shannon <- ggbetweenstats(div, x = type, y = Shannon)
shannon <- shannon + labs(x = "type", y = "Shannon_estimates", title = paste0("Standardized to total reads, N=", nsamples(normalized_walaa))) +  
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
boxplot3 <- ggplot(div, aes(x=type, y=Chao1)) + geom_boxplot() + geom_jitter(aes(color=sample)) +xlab("") + ylab("Chao1_Richness") + ggtitle("Standardized to total reads, N=8")
boxplot3
boxplot3 + stat_compare_means(method = "t.test", label.y = 0.8)


# plot ordination graph (Beta-diversity, NMDS)
set.seed(20)
UAE_ord <- ordinate(normalized_walaa, "NMDS", "bray", k=2, trymax=100)
UAE_ord

Call:
  metaMDS(comm = veganifyOTU(physeq), distance = distance, k = 2,      trymax = 100) 
global Multidimensional Scaling using monoMDS
Data:     wisconsin(sqrt(veganifyOTU(physeq))) 
Distance: bray 
Dimensions: 2 
Stress:     0.0287 
Stress type 1, weak ties
Best solution was repeated 5 times in 20 tries
The best solution was from try 13 (random start)
Scaling: centring, PC rotation 
Species: expanded scores based on ‘wisconsin(sqrt(veganifyOTU(physeq)))’  

stressplot(UAE_ord)
legend("top", "stress = 0.0287", bty = "n", cex = 1)
goodness(UAE_ord)
[1] 0.00856 0.00594 0.01200 0.01292 0.00477 0.01020 0.01364 0.00971

color = c("type")
shape = c("sample")
title = c("NMDS of 16s rRNA microbiome, Bray-Curtis distance, K=2")
Ordination <- plot_ordination(normalized_walaa, UAE_ord, type="sample", color=color, 
                              shape=shape, title=title) 
ordination = Ordination + theme(axis.text = element_text(size = 16, face = "bold"), axis.title = element_text(size = 18, face = "bold"), 
             legend.title = element_text(size = 14)) + theme_bw() + labs(color = color, shape = shape) + geom_point(size=2) + 
             geom_text(aes(label=rownames(sample_df), hjust=0.4, vjust=-0.5)) 
ordination

#ANOSIM
presabs = get_variable(normalized_walaa)
presabs
region_ano = anosim(phyloseq::distance(normalized_walaa, "bray"), presabs$type)
region_ano
summary(region_ano)

##ANOSIM analysis
pc = read.csv("D:/Dr. walaa/New folder/ASVs_counts.csv", header= TRUE)
com = pc[,1:ncol(pc)]
m_com = as.matrix(com)
ano = anosim(m_com, pc$sample, distance = "bray", permutations = 9999)
ano


####### Rarefaction curve################
install.packages("permute")
install.packages("lattice")
library("vegan")
library("permute")
library("lattice")

# import the count table
Data=read.table("D:/Dr. walaa/New folder/ASVs_counts.csv")
Data_sub=as.data.frame(Data[,c(1:8)])
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
rarefied_microbes <- rarefy_even_depth(metabarcoding)
rarecurve(t(otu_table(normalized_walaa)), step = 8, cex = .6, col = "red")
rarecurve(metabarcoding, step = 20, col = "blue", cex = 0.6, label = FALSE)


#######################################################################################################################################################################################################################################################
## Visualize the relative abundance for taxa > 100 ##
# create phyloseq file (short file after discard the taxa lower than 100). three tables are required (ASV table, Taxonomy table and Sample info table - these can be all in single excel file on three separate sheets) #####

asv_mat_surf<-as.data.frame(read_excel("D:/Dr. walaa/New folder/phyloseq _relative.xlsx", sheet = "ASV_table"))
tax_mat_surf<-as.data.frame(read_excel("D:/Dr. walaa/New folder/phyloseq _relative.xlsx", sheet = "taxa_table"))
sample_df_surf<-as.data.frame(read_excel("D:/Dr. walaa/New folder/phyloseq _relative.xlsx", sheet = "sample"))

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

## data preprocessing (sums, maximum, minimum)
sample_sums(data_surf)
max(sample_sums(data_surf))
min(sample_sums(data_surf))

## Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(data_surf))
standf = function(x, t=total) round(t * (x / sum(x)))
normalized_surf = transform_sample_counts(data_surf, standf)
normalized_surf


plot_bar= plot_bar(data_surf, x = 'type', fill = "Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity", position = "stack") +
  labs(x='sample', y='Relative abundance') + theme_bw() + scale_y_percent() +ggtitle("Composition of bacteria, N = 8") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10))
plot_bar
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
wl_tr = transform_sample_counts(data_surf, function(x) 100 * x/sum(x)/4)
preod=c("Acidimicrobiia","Actinomycetia","Alphaproteobacteria","Anaerolineae","Bacilli","Bacteroidia","Binatia","Blastocatellia","Chloroflexia","Deinococci","Gammaproteobacteria","Ellin6529","Clostridia","Dehalococcoidia","Methylomirabilia","Gemmatimonadetes","Thermoleophilia","Polyangia","UBA4738","Rhodothermia","Vicinamibacteria","Rubrobacteria","Sericytochromatia","UBA1144","Saccharimonadia","Others","NA") 
col <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
         "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "gold", "beige", "azure1", "slategray1", "yellow", "plum", "pink", "green", "black", "lavender", "blue")

##### Phyla level #####
plot_bar(wl_tr, x = "type", fill=factor("Phylum", preod)) +
  geom_bar(aes(fill = factor(Phylum, preod)), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4","darksalmon"))+
  theme_bw() + ylab("Reletive abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=0))
plot_bar(wl_tr, x = "area", fill=factor("Phylum", preod)) +
  geom_bar(aes(fill = factor(Phylum, preod)), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4","darksalmon"))+
  theme_bw() + ylab("Reletive abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=0))

##### class level #####
plot_bar(wl_tr, x = "type", fill=factor("Class", preod)) +
  geom_bar(aes(fill = factor(Class, preod)), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4","darksalmon","gold", "beige", "azure1", "slategray1", "yellow", "plum", "pink", "green", "black", "lavender", "blue","darkgray"))+
  theme_bw() + ylab("Reletive abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=0))
plot_bar(wl_tr, x = "area", fill=factor("Class", preod)) +
  geom_bar(aes(fill = factor(Class, preod)), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4","darksalmon","gold", "beige", "azure1", "slategray1", "yellow", "plum", "pink", "green", "black", "lavender", "blue","darkgray"))+
  theme_bw() + ylab("Reletive abundance (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 10, angle=0))


packageVersion("dada2")

