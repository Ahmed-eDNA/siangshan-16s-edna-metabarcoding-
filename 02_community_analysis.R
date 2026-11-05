# 16S rRNA microbial community analysis
# Input: ASV count table, taxonomy table, and sample metadata.
# Rows in counts/taxonomy = ASVs; columns in counts = samples.

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(vegan)
  library(tibble)
})

# --------------------------- USER CONFIGURATION ---------------------------

counts_file <- "data/example/asv_counts.tsv"
taxonomy_file <- "data/example/taxonomy.tsv"
metadata_file <- "data/example/metadata.tsv"
output_dir <- "results/community_analysis"

# Set to a taxonomy rank present in taxonomy.tsv.
composition_rank <- "Phylum"

# Metadata variable used to group alpha-diversity plots and color NMDS points.
group_variable <- "layer"

# Number of permutations for PERMANOVA and dispersion test.
permutations <- 999

set.seed(2026)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------- IMPORT -----------------------------------

read_tsv_with_ids <- function(path, id_name) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                  row.names = 1, na.strings = c("", "NA"))
  if (nrow(x) == 0L || ncol(x) == 0L) stop("Empty input table: ", path)
  x
}

counts <- read_tsv_with_ids(counts_file, "ASV_ID")
taxonomy <- read_tsv_with_ids(taxonomy_file, "ASV_ID")
metadata <- read_tsv_with_ids(metadata_file, "SampleID")

counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))
count_matrix <- as.matrix(counts)
storage.mode(count_matrix) <- "numeric"

if (anyNA(count_matrix) || any(count_matrix < 0) ||
    any(count_matrix != floor(count_matrix))) {
  stop("ASV counts must be non-missing, non-negative integers.")
}
if (is.null(rownames(count_matrix)) || is.null(colnames(count_matrix))) {
  stop("Count table requires ASV row names and sample column names.")
}
if (!setequal(rownames(count_matrix), rownames(taxonomy))) {
  stop("ASV identifiers differ between count and taxonomy tables.")
}
if (!setequal(colnames(count_matrix), rownames(metadata))) {
  stop("Sample identifiers differ between count table and metadata.")
}
if (!group_variable %in% colnames(metadata)) {
  stop("Grouping variable not found in metadata: ", group_variable)
}
if (!composition_rank %in% colnames(taxonomy)) {
  stop("Taxonomic rank not found in taxonomy table: ", composition_rank)
}

# Reorder rows/columns to ensure exact matching.
taxonomy <- taxonomy[rownames(count_matrix), , drop = FALSE]
metadata <- metadata[colnames(count_matrix), , drop = FALSE]
taxonomy_matrix <- as.matrix(taxonomy)
storage.mode(count_matrix) <- "integer"

# ---------------------------- PHYLOSEQ OBJECT ------------------------------

OTU <- otu_table(count_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy_matrix)
SAM <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, SAM)

saveRDS(ps, file.path(output_dir, "phyloseq_object.rds"))

# ------------------------------ ALPHA DIVERSITY ----------------------------

alpha <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
alpha$SampleID <- rownames(alpha)
alpha <- alpha %>%
  left_join(
    metadata %>% rownames_to_column("SampleID"),
    by = "SampleID"
  )
write.csv(alpha, file.path(output_dir, "alpha_diversity.csv"), row.names = FALSE)

alpha_long <- alpha %>%
  pivot_longer(
    cols = c("Observed", "Chao1", "Shannon"),
    names_to = "Metric", values_to = "Value"
  )

p_alpha <- ggplot(alpha_long, aes(x = .data[[group_variable]], y = Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(width = 0.12, alpha = 0.75, size = 2) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(x = group_variable, y = "Diversity estimate",
       title = "Alpha diversity of sediment-associated microbial communities") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(output_dir, "alpha_diversity.png"), p_alpha,
       width = 9, height = 5.5, dpi = 300)

# -------------------------- RELATIVE COMPOSITION ---------------------------

ps_relative <- transform_sample_counts(ps, function(x) x / sum(x))
composition <- psmelt(ps_relative)
composition[[composition_rank]] <- as.character(composition[[composition_rank]])
composition[[composition_rank]][is.na(composition[[composition_rank]]) |
                                  composition[[composition_rank]] == ""] <- "Unclassified"

composition_summary <- composition %>%
  group_by(.data[[group_variable]], .data[[composition_rank]]) %>%
  summarise(mean_relative_abundance = mean(Abundance), .groups = "drop")
write.csv(composition_summary,
          file.path(output_dir, paste0("relative_abundance_", composition_rank, ".csv")),
          row.names = FALSE)

p_composition <- ggplot(
  composition_summary,
  aes(x = .data[[group_variable]], y = mean_relative_abundance,
      fill = .data[[composition_rank]])
) +
  geom_col(position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = group_variable, y = "Mean relative abundance",
       fill = composition_rank,
       title = paste("Taxonomic composition at", composition_rank, "level")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(output_dir, paste0("composition_", composition_rank, ".png")),
       p_composition, width = 10, height = 6, dpi = 300)

# -------------------------- BETA DIVERSITY / NMDS --------------------------

# Bray-Curtis on within-sample relative abundance.
relative_matrix <- as(otu_table(ps_relative), "matrix")
if (!taxa_are_rows(ps_relative)) relative_matrix <- t(relative_matrix)
bray <- vegdist(t(relative_matrix), method = "bray")
write.csv(as.matrix(bray), file.path(output_dir, "bray_curtis_distance.csv"))

nmds <- metaMDS(
  t(relative_matrix),
  distance = "bray",
  k = 2,
  trymax = 100,
  autotransform = FALSE,
  trace = FALSE
)
writeLines(capture.output(nmds), file.path(output_dir, "nmds_summary.txt"))

ordination <- as.data.frame(scores(nmds, display = "sites"))
ordination$SampleID <- rownames(ordination)
ordination <- left_join(
  ordination,
  metadata %>% rownames_to_column("SampleID"),
  by = "SampleID"
)
p_nmds <- ggplot(
  ordination,
  aes(x = NMDS1, y = NMDS2, color = .data[[group_variable]])
) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "NMDS ordination based on Bray–Curtis dissimilarity",
    subtitle = paste("Stress =", signif(nmds$stress, 3)),
    color = group_variable
  ) +
  coord_equal() +
  theme_bw()
ggsave(file.path(output_dir, "nmds_bray_curtis.png"), p_nmds,
       width = 7, height = 6, dpi = 300)

# PERMANOVA and multivariate-dispersion check.
# Add strata=... to adonis2 if the design requires restricted permutations.
group <- droplevels(factor(metadata[[group_variable]]))
if (nlevels(group) > 1L) {
  permanova <- adonis2(bray ~ group, data = metadata,
                       permutations = permutations, by = "terms")
  dispersion <- betadisper(bray, group)
  dispersion_test <- permutest(dispersion, permutations = permutations)

  capture.output(permanova,
                 file = file.path(output_dir, "permanova_results.txt"))
  capture.output(dispersion_test,
                 file = file.path(output_dir, "dispersion_test_results.txt"))
} else {
  warning("PERMANOVA skipped: grouping variable has fewer than two levels.")
}

# ------------------------------- RAREFACTION -------------------------------

# vegan::rarecurve expects samples in rows and taxa in columns.
sample_by_asv <- t(count_matrix)
pdf(file.path(output_dir, "rarefaction_curves.pdf"), width = 9, height = 7)
rarecurve(sample_by_asv, step = 100, label = FALSE,
          xlab = "Sequencing depth (reads)", ylab = "Observed ASVs")
dev.off()

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))
message("Community analysis complete. Outputs saved to: ",
        normalizePath(output_dir))
