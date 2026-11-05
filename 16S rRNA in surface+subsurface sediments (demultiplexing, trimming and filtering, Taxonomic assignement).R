                                                   ################ DADA2 pipeline ###################
                                                   
## install DADA2 package                                                   
(1)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.11")
(2)
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

## Define the path variable so that it points to the extracted directory on my PC:
data_path <- "E:/sequencing data/microbiota"
list.files(data_path) # forward (R1) and reverse (R2) fastq filenames have format: samplename_R1_clean.fastq and samplename_R2_clean.fastq.

## read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.
file.exists(data_path)
Fwds<-sort(list.files(data_path, pattern="R1.clean.fastq", full.names = TRUE))
Revs<-sort(list.files(data_path, pattern="R2.clean.fastq", full.names = TRUE))
Fwds
Revs

## Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.
sample.names <- sapply(strsplit(basename(Fwds), "_"), '[', 1)
sample.names
sample.names <- sapply(strsplit(basename(Revs), "_"), '[', 1)
sample.names
sample.names <- sapply(strsplit(basename(Fwds), ".R"), '[', 1)
sample.names

## Inspect read quality profiles & visualize it.
plotQualityProfile(Fwds[1])
plotQualityProfile(Fwds[1:49])
plotQualityProfile(Revs[1])
plotQualityProfile(Revs[1:49])

## data filtering and  Place the filtered files in filtered/ subdirectory.
fnF <- file.path(data_path, Fwds)
fnR <- file.path(data_path, Revs)
filtFs <-file.path(data_path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <-file.path(data_path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
filtFs
filtRs

## Filter and trim (remove primers (trimleft= remove the reads from beginning of reads, trunclen= setting the minimum size to trim the forward and reverse reads to in order to keep the quality scores roughly above 30 overall, rm.phix= removes any reads that match phix bacteriophage genome that added to illumina sequencing runs for quality monitoring).
out<- filterAndTrim(Fwds, filtFs, Revs, filtRs, maxN= 0, maxEE = c(2,2), truncQ = 2, minLen=170, truncLen = c(240, 210), rm.phix = TRUE, compress = TRUE, multithread = F)
head(out)
dim(out)
out
class(out) # matrix

## Now let's take a look at our filtered reads (filtFs & filtRs):
plotQualityProfile(filtFs[1])
plotQualityProfile(filtFs[1:49])
plotQualityProfile(filtRs[1])
plotQualityProfile(filtRs[1:49])

## learn and plot error rates, using a machine learning algorithm.
errF <- learnErrors(filtFs, multithread = T)
errR <- learnErrors(filtRs, multithread = T)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
gc()

## Dereplication
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- samples

## sample inference (apply the core sample inference algorithm to the filtered and trimmed sequence data)- where the core sample inference algorithm is applied.
dadaFs <- dada(derep_forward, err = errF, multithread=TRUE)
dadaRs <- dada(derep_reverse, err = errR, multithread=TRUE)
dadaFs[[1]]
?dada
## Merge paired reads aka align forward and reverse complement reads.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
head(mergers[[1]])

## Construct amplicon sequence variants (ASVs) table (a higher-resolution version of the OTU table produced by traditional methods).
seqtable <- makeSequenceTable(mergers)
dim(seqtable)
table(nchar(getSequences(seqtable)))

## Remove sequences that are shorter than 200bp or longer than 300bp as they are likely results of non-specific priming.
seqtab2 <- seqtable[,nchar(colnames(seqtable)) %in% 200:300]
dim(seqtab2)

## Remove chimeras
seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab2.nochim)
sum(seqtab2.nochim)/sum(seqtab2)

## track reads throughout the above pipeline and print out number of reads that make it through each step.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "E:/sequencing data/microbiota/microbiota_trim_new.csv")
track
## make output fasta file with sequences and associated ASV ID strings.
library(ShortRead)
seq <- getSequences(seqtab2.nochim)
writeFasta(seq, file="E:/sequencing data/microbiota/dada2 output/microbiota.fasta")

## taxonomic assignment using DECIPHER Bioconductor package (recent developed IdTaxa taxonomic classification method).
# downloading DECIPHER-formatted SILVA v138 reference.
download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData") 
# loading reference taxonomy object.
load("SILVA_SSU_r138_2019.RData")
#loading DECIPHER.
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")
library(DECIPHER)
packageVersion("DECIPHER")
# creating DNAStringSet object of our ASVs.
dna <- DNAStringSet(getSequences(seqtab2.nochim))
# classifying.
tax_info <- IdTaxa(dna, trainingSet=trainingSet, strand="both", processors=NULL)
# giving our seq headers more manageable names (ASV_1, ASV_2...).
asv_seq <- colnames(seqtab2.nochim)
asv_header <- vector(dim(seqtab2.nochim)[2], mode="character")
for (i in 1:dim(seqtab2.nochim)[2]) {
  asv_header[i] <- paste(">ASV", i, sep="_")
}
# making and writing out a fasta of our final ASVs.
asv_fasta <- c(rbind(asv_header, asv_seq))
write(asv_fasta, "ASVs.fa")
write(asv_fasta, "ASVs.csv")

## Count table.
asv_tab <- t(seqtab2.nochim)
row.names(asv_tab) <- sub(">", "", asv_header)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tab, "ASVs_counts.csv")

## tax table.
# creating table of taxonomy and setting any that are unclassified as "NA".
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_header)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

## write.table(longDF and shortDF).
longDF <- cbind(asv_tax, t(seqtab2.nochim))
write.csv(longDF, "E:/sequencing data/microbiota/dada2 output/microbiota_asv_table_long", sep = "\t")
write.csv(longDF, "E:/sequencing data/microbiota/dada2 output/decipher/microbiota_asv_table_long.csv")

longDF <- read.csv("E:/sequencing data/microbiota/dada2 output/decipher/microbiota_asv_table_long.csv",header=T)
taxa <- aggregate(longDF[2:8], longDF["class"], head, 1)
taxa <- taxa[-1]
counts <- aggregate(longDF[9:ncol(longDF)], longDF["class"], sum)
shortDF <- merge(taxa, counts, by = "class")
write.csv(shortDF, "E:/sequencing data/microbiota/dada2 output/decipher/microbiota_asv_table_short",
          sep = "\t")
write.csv(shortDF, "E:/sequencing data/microbiota/dada2 output/decipher/microbiota_asv_table_short.csv")