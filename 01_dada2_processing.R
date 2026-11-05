# 16S rRNA paired-end processing with DADA2 and DECIPHER
# Adapted and documented from the original Xiang-Shan Wetland workflow.
# Review all settings against read quality and target amplicon before use.

suppressPackageStartupMessages({
  library(dada2)
  library(DECIPHER)
  library(Biostrings)
  library(ShortRead)
})

# --------------------------- USER CONFIGURATION ---------------------------

input_dir <- "data/fastq"       # directory containing paired FASTQ files
output_dir <- "results/dada2"
silva_training_set <- "references/SILVA_SSU_r138_2019.RData"

# Original-study starting parameters; not universal defaults.
trunc_len <- c(240, 210)
max_ee <- c(2, 2)
trunc_q <- 2
min_read_length <- 170
min_amplicon_length <- 200
max_amplicon_length <- 300
remove_phix <- TRUE
multithread <- FALSE            # set TRUE only if supported by your system

# Expected filename pattern: SAMPLE_R1.clean.fastq[.gz] and SAMPLE_R2...
forward_pattern <- "_R1\\.clean\\.fastq(\\.gz)?$"
reverse_pattern <- "_R2\\.clean\\.fastq(\\.gz)?$"

# IdTaxa confidence threshold (DECIPHER default is 60; state chosen value).
taxonomy_bootstrap <- 60

# ------------------------------- VALIDATION -------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
filtered_dir <- file.path(output_dir, "filtered")
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_dir)) stop("FASTQ input directory does not exist: ", input_dir)
if (!file.exists(silva_training_set)) {
  stop("SILVA training-set file not found: ", silva_training_set,
       "\nDownload the DECIPHER-formatted SILVA SSU r138 training set and update the path.")
}

fastq_files <- list.files(input_dir, pattern = "\\.fastq(\\.gz)?$", full.names = TRUE)
forward_files <- sort(grep(forward_pattern, fastq_files, value = TRUE))
reverse_files <- sort(grep(reverse_pattern, fastq_files, value = TRUE))

if (length(forward_files) == 0L || length(reverse_files) == 0L) {
  stop("No paired FASTQ files found. Check the filename convention and patterns.")
}

sample_id <- function(x, pattern) sub(pattern, "", basename(x))
forward_ids <- sample_id(forward_files, forward_pattern)
reverse_ids <- sample_id(reverse_files, reverse_pattern)

if (anyDuplicated(forward_ids) || anyDuplicated(reverse_ids)) {
  stop("Duplicate sample identifiers detected in FASTQ filenames.")
}
if (!setequal(forward_ids, reverse_ids)) {
  stop("Forward/reverse sample IDs do not match.\nForward only: ",
       paste(setdiff(forward_ids, reverse_ids), collapse = ", "),
       "\nReverse only: ", paste(setdiff(reverse_ids, forward_ids), collapse = ", "))
}

# Reorder reverse reads to match forward-read order.
reverse_files <- reverse_files[match(forward_ids, reverse_ids)]
names(forward_files) <- forward_ids
names(reverse_files) <- forward_ids
samples <- forward_ids

message("Matched ", length(samples), " paired-end samples.")

# ---------------------------- QUALITY INSPECTION ---------------------------

pdf(file.path(output_dir, "quality_profiles_raw.pdf"), width = 10, height = 7)
print(plotQualityProfile(forward_files))
print(plotQualityProfile(reverse_files))
dev.off()

# --------------------------- FILTER AND TRIM -------------------------------

filtered_forward <- file.path(filtered_dir, paste0(samples, "_F_filt.fastq.gz"))
filtered_reverse <- file.path(filtered_dir, paste0(samples, "_R_filt.fastq.gz"))
names(filtered_forward) <- samples
names(filtered_reverse) <- samples

filter_stats <- filterAndTrim(
  fwd = forward_files,
  filt = filtered_forward,
  rev = reverse_files,
  filt.rev = filtered_reverse,
  maxN = 0,
  maxEE = max_ee,
  truncQ = trunc_q,
  minLen = min_read_length,
  truncLen = trunc_len,
  rm.phix = remove_phix,
  compress = TRUE,
  multithread = multithread
)

write.csv(filter_stats, file.path(output_dir, "filtering_summary.csv"))

pdf(file.path(output_dir, "quality_profiles_filtered.pdf"), width = 10, height = 7)
print(plotQualityProfile(filtered_forward))
print(plotQualityProfile(filtered_reverse))
dev.off()

# -------------------------- ERROR LEARNING / DENOISING ---------------------

err_forward <- learnErrors(filtered_forward, multithread = multithread)
err_reverse <- learnErrors(filtered_reverse, multithread = multithread)

pdf(file.path(output_dir, "learned_error_profiles.pdf"), width = 10, height = 7)
print(plotErrors(err_forward, nominalQ = TRUE))
print(plotErrors(err_reverse, nominalQ = TRUE))
dev.off()

derep_forward <- derepFastq(filtered_forward, verbose = TRUE)
derep_reverse <- derepFastq(filtered_reverse, verbose = TRUE)
names(derep_forward) <- samples
names(derep_reverse) <- samples

dada_forward <- dada(derep_forward, err = err_forward, multithread = multithread)
dada_reverse <- dada(derep_reverse, err = err_reverse, multithread = multithread)

# ------------------------------- MERGE READS -------------------------------

merged <- mergePairs(
  dada_forward, filtered_forward,
  dada_reverse, filtered_reverse,
  verbose = TRUE
)

sequence_table <- makeSequenceTable(merged)
write.csv(as.data.frame(table(nchar(getSequences(sequence_table)))),
          file.path(output_dir, "merged_sequence_length_distribution.csv"),
          row.names = FALSE)

# Retain the amplicon length range used in the original script.
sequence_table_filtered <- sequence_table[
  , nchar(colnames(sequence_table)) >= min_amplicon_length &
    nchar(colnames(sequence_table)) <= max_amplicon_length,
  drop = FALSE
]

sequence_table_nochim <- removeBimeraDenovo(
  sequence_table_filtered,
  method = "consensus",
  multithread = multithread,
  verbose = TRUE
)

# ------------------------------- READ TRACKING -----------------------------

get_n <- function(x) sum(getUniques(x))
track <- cbind(
  input = filter_stats[, "reads.in"],
  filtered = filter_stats[, "reads.out"],
  denoised_forward = vapply(dada_forward, get_n, numeric(1)),
  denoised_reverse = vapply(dada_reverse, get_n, numeric(1)),
  merged = vapply(merged, get_n, numeric(1)),
  nonchim = rowSums(sequence_table_nochim)
)
rownames(track) <- samples
write.csv(track, file.path(output_dir, "read_tracking.csv"))

# -------------------------- TAXONOMIC ASSIGNMENT ---------------------------

load(silva_training_set) # Expected to provide DECIPHER object `trainingSet`.
if (!exists("trainingSet")) {
  stop("The loaded reference file did not create an object named 'trainingSet'.")
}

asv_sequences <- colnames(sequence_table_nochim)
dna <- DNAStringSet(asv_sequences)
names(dna) <- paste0("ASV_", seq_along(dna))

taxonomy <- IdTaxa(
  dna,
  trainingSet = trainingSet,
  strand = "both",
  processors = NULL,
  threshold = taxonomy_bootstrap
)

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxonomy_matrix <- t(vapply(taxonomy, function(x) {
  idx <- match(tolower(ranks), tolower(x$rank))
  taxon <- x$taxon[idx]
  taxon[is.na(taxon) | grepl("^unclassified_", taxon)] <- NA_character_
  taxon
}, character(length(ranks))))

colnames(taxonomy_matrix) <- ranks
rownames(taxonomy_matrix) <- names(dna)

# ------------------------------- EXPORTS -----------------------------------

asv_counts <- t(sequence_table_nochim)
rownames(asv_counts) <- names(dna)
colnames(asv_counts) <- samples

write.table(
  asv_counts,
  file = file.path(output_dir, "asv_counts.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)
write.table(
  taxonomy_matrix,
  file = file.path(output_dir, "taxonomy.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)
writeXStringSet(dna, filepath = file.path(output_dir, "asv_sequences.fasta"))
write.csv(as.data.frame(track), file.path(output_dir, "read_tracking.csv"))

saveRDS(
  list(
    sequence_table = sequence_table_nochim,
    taxonomy = taxonomy_matrix,
    read_tracking = track,
    parameters = list(
      trunc_len = trunc_len, max_ee = max_ee, trunc_q = trunc_q,
      min_read_length = min_read_length,
      amplicon_length = c(min_amplicon_length, max_amplicon_length),
      taxonomy_bootstrap = taxonomy_bootstrap
    )
  ),
  file = file.path(output_dir, "dada2_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("DADA2 processing complete. Outputs saved to: ", normalizePath(output_dir))
