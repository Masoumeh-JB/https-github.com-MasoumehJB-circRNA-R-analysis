# ============================================================
# Project: circRNA network analysis
# Script : circRNA.R
# Purpose: CircRNA abundance in breast cancer subtypes by CIRIquant
# Triple Negative (GSE142731): TNBC
# HER2 Positive (GSE191230): HER2
# Inflammatory Breast Cancer (GSE207248): IBC
# Author : Masoumeh Jomhori Baloch
# Date   : <2026>
# ============================================================

###############################################################################
# Manuscript section: 3.1 — circRNA abundance in human Breast Cancer tissues
#
# This block contains the complete analysis pipeline for Section 3.1
# (Stages 1–12). The next manuscript section (e.g., 3.2) starts after the
# "END OF SECTION 3.1" delimiter below to avoid mixing outputs or logic.
#
# Stages included:
#   Stage 1  — Filter circRNA counts (>=2 reads in >=1 sample)
#   Stage 2  — SRPBM normalization + filters
#   Stage 3  — QC (depth/burden, PCA/MDS)
#   Stage 4  — BED generation (hg38)
#   Stage 5  — Database overlap (circBase/circAtlas/MiOncoCirc)
#   Stage 6  — Read support distribution + donut chart
#   Stage 7  — Exon structure + chi-square test
#   Stage 8  — Chromosomal SRPBM density
#   Stage 9  — Venn diagram (shared high-confidence circRNAs)
#   Stage 10 — SRPBM/FPKM integration + heatmap + PCC
#   Stage 11 — circRNA–miRNA–mRNA network + Cytoscape export
#   Stage 12 — GEO/DAVID enrichment plotting
###############################################################################

# ==============================
# START OF SECTION 3.1 (Stages 1–12)
# ==============================

# ============================================================
# Stage 1 — Filter circRNA counts (>=2 reads in >=1 sample)
# Inputs :
#   data/raw/<GSE>/expressionsfinal.txt
# Outputs:
#   data/processed/stage01/<GSE>_filtered_counts_ge2.csv
#   data/processed/stage01/<GSE>_filtered_counts_ge2_preview.tsv
#   outputs/tables/stage01_circRNA_counts_before_after_ge2_filter.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(here)
})

dir.create(here("data", "processed", "stage01"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)

# ---- Inputs (EDIT these paths to match your repo raw-data layout) ----
paths_counts <- list(
  GSE142731 = here("data", "raw", "GSE142731", "expressionsfinal.txt"),
  GSE191230 = here("data", "raw", "GSE191230", "expressionsfinal.txt"),
  GSE207248 = here("data", "raw", "GSE207248", "expressionsfinal.txt")
)

filter_counts_ge2 <- function(file_path, label) {
  cat("\n=====================================\n")
  cat("Stage1 | Processing:", label, "\n")
  cat("File:", file_path, "\n")
  cat("=====================================\n")
  
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  
  Counts <- read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(Counts) == 0) stop("Counts table is empty: ", file_path)
  
  total_raw <- nrow(Counts)
  
  # Ensure circ_id exists
  if (!"circ_id" %in% colnames(Counts)) {
    colnames(Counts)[1] <- "circ_id"
    cat("circ_id not found. First column assumed circ_id.\n")
  }
  
  # Clean circ_id
  Counts <- Counts %>% mutate(circ_id = gsub("\\|", "-", as.character(circ_id)))
  
  numeric_cols <- setdiff(colnames(Counts), "circ_id")
  if (length(numeric_cols) == 0) stop("No sample columns found in: ", file_path)
  
  Counts[, numeric_cols] <- lapply(Counts[, numeric_cols, drop = FALSE], function(x) {
    suppressWarnings(as.numeric(x))
  })
  
  ge2_flag <- rowSums(Counts[, numeric_cols, drop = FALSE] >= 2, na.rm = TRUE) > 0
  filtered_counts <- Counts[ge2_flag, , drop = FALSE]
  total_ge2 <- nrow(filtered_counts)
  
  out_csv <- here("data", "processed", "stage01", paste0(label, "_filtered_counts_ge2.csv"))
  out_tsv <- here("data", "processed", "stage01", paste0(label, "_filtered_counts_ge2_preview.tsv"))
  
  write.csv(filtered_counts, out_csv, row.names = FALSE)
  write_tsv(head(filtered_counts, 20), out_tsv)
  
  cat("Saved:", out_csv, " | kept:", total_ge2, "/", total_raw, "\n")
  
  tibble(
    GSE = label,
    Total_raw = total_raw,
    Total_ge2 = total_ge2,
    Percent_kept = ifelse(total_raw > 0, round(100 * total_ge2 / total_raw, 2), NA_real_)
  )
}

summary_counts <- bind_rows(lapply(names(paths_counts), function(label) {
  filter_counts_ge2(paths_counts[[label]], label)
}))

out_summary <- here("outputs", "tables", "stage01_circRNA_counts_before_after_ge2_filter.csv")
write.csv(summary_counts, out_summary, row.names = FALSE)

cat("\n===== Stage1 Summary (>=2 reads filter) =====\n")
print(summary_counts)
cat(" Saved summary:", out_summary, "\n")

# ============================================================
# Stage 2 — SRPBM normalization + filters
# Filters:
#   1) mean_SRPBM >= 1
#   2) present in >=2 samples (SRPBM > 0)
# Inputs :
#   data/processed/stage01/<GSE>_filtered_counts_ge2.csv
#   data/raw/<GSE>/mapped_reads*.csv  (2 columns: Sample, Total_Mapped)
# Outputs:
#   data/processed/stage02/<GSE>_filtered_counts_ge2_SRPBM_ge1.csv
#   data/processed/stage02/<GSE>_filtered_counts_ge2_SRPBM_ge1_preview.tsv
#   outputs/tables/stage02_circRNA_counts_after_SRPBM_filter_hg38_only.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(here)
})

dir.create(here("data", "processed", "stage02"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)

paths_counts_ge2 <- list(
  GSE142731 = here("data", "processed", "stage01", "GSE142731_filtered_counts_ge2.csv"),
  GSE191230 = here("data", "processed", "stage01", "GSE191230_filtered_counts_ge2.csv"),
  GSE207248 = here("data", "processed", "stage01", "GSE207248_filtered_counts_ge2.csv")
)

# ---- mapped reads inputs (EDIT to match your raw-data layout) ----
paths_mapped_reads <- list(
  GSE142731 = here("data", "raw", "GSE142731", "mapped_reads.csv"),
  GSE191230 = here("data", "raw", "GSE191230", "mapped_reads.csv"),
  GSE207248 = here("data", "raw", "GSE207248", "mapped_reads_Samples.csv")
)

subtype_map <- c(GSE142731 = "TNBC", GSE191230 = "HER2+", GSE207248 = "IBC")

quick_file_check <- function(fp, label = "") {
  if (!file.exists(fp)) stop("File not found: ", fp)
  sz <- file.info(fp)$size
  if (is.na(sz) || sz < 5) stop("File exists but looks empty: ", fp)
  cat("  - OK:", label, "|", fp, "| size =", sz, "bytes\n")
}

clean_sample <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace(x, "\\.sorted\\.bam$", "")
  x <- str_replace(x, "\\.bam$", "")
  x <- str_replace(x, "\\.fastq\\.gz$", "")
  x <- str_replace(x, "\\.fq\\.gz$", "")
  x <- str_replace_all(x, "[[:space:]]+", "")
  x
}

read_mapped_reads <- function(fp) {
  # try csv, then tsv, then semicolon
  mr <- suppressWarnings(tryCatch(read_csv(fp, show_col_types = FALSE), error = function(e) NULL))
  if (is.null(mr) || ncol(mr) < 2) {
    mr <- suppressWarnings(tryCatch(read_delim(fp, delim = "\t", show_col_types = FALSE), error = function(e) NULL))
  }
  if (is.null(mr) || ncol(mr) < 2) {
    mr <- suppressWarnings(tryCatch(read_delim(fp, delim = ";", show_col_types = FALSE), error = function(e) NULL))
  }
  if (is.null(mr) || ncol(mr) < 2) stop("Could not read mapped_reads as 2-column table: ", fp)
  
  mr <- mr[, 1:2]
  colnames(mr) <- c("Sample", "Total_Mapped")
  mr %>%
    mutate(
      Sample = clean_sample(Sample),
      Total_Mapped = suppressWarnings(as.numeric(Total_Mapped))
    )
}

normalize_srpbm_one <- function(label) {
  cat("\n=====================================\n")
  cat("Stage2 | Processing:", label, "(", subtype_map[[label]], ")\n")
  cat("=====================================\n")
  
  counts_fp <- paths_counts_ge2[[label]]
  mapped_fp <- paths_mapped_reads[[label]]
  
  quick_file_check(counts_fp, "Stage1 counts_ge2")
  quick_file_check(mapped_fp, "mapped_reads")
  
  Counts <- read.csv(counts_fp, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"circ_id" %in% colnames(Counts)) colnames(Counts)[1] <- "circ_id"
  if (nrow(Counts) == 0) stop("Stage1 output is empty: ", counts_fp)
  
  sample_cols <- setdiff(colnames(Counts), "circ_id")
  Counts[, sample_cols] <- lapply(Counts[, sample_cols, drop = FALSE], function(x) suppressWarnings(as.numeric(x)))
  
  n_counts_ge2 <- nrow(Counts)
  mapped_reads <- read_mapped_reads(mapped_fp)
  
  counts_samples_clean <- clean_sample(sample_cols)
  matched <- counts_samples_clean %in% mapped_reads$Sample
  n_matched <- sum(matched)
  
  cat("  Matched samples:", n_matched, "/", length(sample_cols), "\n")
  if (n_matched == 0) {
    cat("Example counts samples:\n"); print(head(sample_cols, 15))
    cat("Example mapped_reads samples:\n"); print(head(mapped_reads$Sample, 15))
    stop("No samples matched. Fix sample names between counts and mapped_reads.")
  }
  
  matched_orig_cols  <- sample_cols[matched]
  matched_clean_cols <- counts_samples_clean[matched]
  
  Counts_norm <- Counts
  
  for (i in seq_along(matched_orig_cols)) {
    s_orig  <- matched_orig_cols[i]
    s_clean <- matched_clean_cols[i]
    total_reads <- mapped_reads$Total_Mapped[mapped_reads$Sample == s_clean][1]
    
    if (is.na(total_reads) || total_reads <= 0) {
      Counts_norm[[s_orig]] <- NA_real_
    } else {
      Counts_norm[[s_orig]] <- (Counts_norm[[s_orig]] * 1e9) / total_reads
    }
  }
  
  # Unmatched -> NA so they don't affect means
  if (any(!matched)) {
    for (s_orig in sample_cols[!matched]) Counts_norm[[s_orig]] <- NA_real_
  }
  
  mean_SRPBM <- rowMeans(Counts_norm[, matched_orig_cols, drop = FALSE], na.rm = TRUE)
  present_in_n_samples <- rowSums(Counts_norm[, matched_orig_cols, drop = FALSE] > 0, na.rm = TRUE)
  
  Counts_norm$mean_SRPBM <- mean_SRPBM
  Counts_norm$present_in_n_samples <- as.integer(present_in_n_samples)
  Counts_norm$matched_samples_n <- as.integer(n_matched)
  
  keep1 <- !is.na(Counts_norm$mean_SRPBM) & Counts_norm$mean_SRPBM >= 1
  Counts_f1 <- Counts_norm[keep1, , drop = FALSE]
  keep2 <- Counts_f1$present_in_n_samples >= 2
  Counts_f2 <- Counts_f1[keep2, , drop = FALSE]
  
  cat("  Input (GE2):", n_counts_ge2, "\n")
  cat("  After mean_SRPBM>=1:", nrow(Counts_f1), "\n")
  cat("  After present>=2:", nrow(Counts_f2), "\n")
  
  out_csv <- here("data", "processed", "stage02", paste0(label, "_filtered_counts_ge2_SRPBM_ge1.csv"))
  out_tsv <- here("data", "processed", "stage02", paste0(label, "_filtered_counts_ge2_SRPBM_ge1_preview.tsv"))
  
  if (nrow(Counts_f2) == 0) {
    dbg <- tibble(
      label = label,
      subtype = subtype_map[[label]],
      n_counts_ge2 = n_counts_ge2,
      n_matched_samples = n_matched,
      n_after_mean_SRPBM_ge1 = nrow(Counts_f1),
      n_after_present_ge2 = 0L,
      note = "EMPTY after filters. Check thresholds OR sample match OR mapped reads."
    )
    write.csv(dbg, here("outputs", "tables", paste0("stage02_", label, "_DEBUG_empty_report.csv")), row.names = FALSE)
    stop("Stage2 produced EMPTY output for ", label, ". See debug report in outputs/tables.")
  }
  
  write.csv(Counts_f2, out_csv, row.names = FALSE)
  write_tsv(head(Counts_f2, 20), out_tsv)
  
  cat("Saved:", out_csv, "\n")
  
  tibble(
    GSE = label,
    Subtype = subtype_map[[label]],
    n_counts_ge2 = n_counts_ge2,
    n_matched_samples = n_matched,
    n_after_mean_SRPBM_ge1 = nrow(Counts_f1),
    n_after_present_ge2 = nrow(Counts_f2)
  )
}

summary_stage2 <- bind_rows(lapply(names(paths_counts_ge2), normalize_srpbm_one))

out_summary <- here("outputs", "tables", "stage02_circRNA_counts_after_SRPBM_filter_hg38_only.csv")
write.csv(summary_stage2, out_summary, row.names = FALSE)

cat("\n===== Stage2 Summary (SRPBM) =====\n")
print(summary_stage2)
cat(" Saved summary:", out_summary, "\n")

# ============================================================
# Stage 3 — QC (depth, burden, depth-vs-detected) + PCA/MDS
# Inputs :
#   data/processed/stage01/<GSE>_filtered_counts_ge2.csv
#   data/raw/<GSE>/mapped_reads*.csv
# Outputs:
#   outputs/figures/QC_*.pdf
#   outputs/tables/QC_*.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(ggplot2)
  library(here)
})

dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)

subtype_colors <- c(
  "HER2+" = "#868686FF",
  "IBC"   = "#CD534CFF",
  "TNBC"  = "#0073C2FF"
)

base_dir_processed <- here("data", "processed", "stage01")

paths_counts_ge2 <- list(
  GSE142731 = here(base_dir_processed, "GSE142731_filtered_counts_ge2.csv"),
  GSE191230 = here(base_dir_processed, "GSE191230_filtered_counts_ge2.csv"),
  GSE207248 = here(base_dir_processed, "GSE207248_filtered_counts_ge2.csv")
)

paths_mapped_reads <- list(
  GSE142731 = here("data", "raw", "GSE142731", "mapped_reads.csv"),
  GSE191230 = here("data", "raw", "GSE191230", "mapped_reads.csv"),
  GSE207248 = here("data", "raw", "GSE207248", "mapped_reads_Samples.csv")
)

read_ge2_counts <- function(path) {
  df <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"circ_id" %in% colnames(df)) stop("circ_id not found: ", path)
  rownames(df) <- df$circ_id
  df$circ_id <- NULL
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
  as.matrix(df)
}

counts_list <- lapply(paths_counts_ge2, read_ge2_counts)
all_circs <- Reduce(union, lapply(counts_list, rownames))

counts_all <- do.call(cbind, lapply(names(counts_list), function(gse) {
  mat <- counts_list[[gse]]
  out <- matrix(0, nrow = length(all_circs), ncol = ncol(mat),
                dimnames = list(all_circs, colnames(mat)))
  out[rownames(mat), colnames(mat)] <- mat
  out
}))

make_meta_for_gse <- function(gse, subtype_label) {
  samples <- colnames(counts_list[[gse]])
  depth <- read.csv(paths_mapped_reads[[gse]], header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  depth <- depth[, 1:2]
  colnames(depth) <- c("Sample", "MappedReads")
  depth$Sample <- as.character(depth$Sample)
  
  tibble(
    sample = samples,
    dataset = gse,
    subtype = subtype_label,
    mapped_reads = depth$MappedReads[match(samples, depth$Sample)]
  )
}

meta <- bind_rows(
  make_meta_for_gse("GSE142731", "TNBC"),
  make_meta_for_gse("GSE191230", "HER2+"),
  make_meta_for_gse("GSE207248", "IBC")
) %>%
  mutate(
    subtype = factor(subtype, levels = c("TNBC", "HER2+", "IBC")),
    dataset = factor(dataset)
  )

qc_df <- tibble(
  sample = colnames(counts_all),
  detected_circs = colSums(counts_all > 0),
  total_counts = colSums(counts_all)
) %>%
  left_join(meta, by = c("sample")) %>%
  mutate(
    depth_million = mapped_reads / 1e6,
    circRNAs_per_10M = detected_circs / (mapped_reads / 1e7)
  )

# 1) Depth boxplot
p_depth <- ggplot(qc_df, aes(x = subtype, y = depth_million, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  scale_fill_manual(values = subtype_colors) +
  labs(x = "Subtype", y = "Mapped reads (millions)", title = "Sequencing depth per sample") +
  theme_classic() + theme(legend.position = "none")

ggsave(here("outputs", "figures", "QC_depth_per_subtype.pdf"), p_depth, width = 5, height = 4, dpi = 1200)

# 2) Burden per 10M
p_burden <- ggplot(qc_df, aes(x = subtype, y = circRNAs_per_10M, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  scale_fill_manual(values = subtype_colors) +
  labs(x = "Subtype", y = "Detected circRNAs per 10M mapped reads", title = "circRNA burden adjusted for depth") +
  theme_classic() + theme(legend.position = "none")

ggsave(here("outputs", "figures", "QC_circRNA_burden_per10M.pdf"), p_burden, width = 5, height = 4, dpi = 1200)

# 3) Scatter depth vs detected
p_scatter <- ggplot(qc_df, aes(x = depth_million, y = detected_circs, color = subtype)) +
  geom_point(alpha = 0.75, size = 2.2) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  scale_color_manual(values = subtype_colors) +
  labs(x = "Mapped reads (millions)", y = "Number of detected circRNAs (>=1 read after GE2 list)",
       title = "Relationship between sequencing depth and detected circRNAs") +
  theme_classic()

ggsave(here("outputs", "figures", "QC_depth_vs_detected_circs.pdf"), p_scatter, width = 5.5, height = 4.5, dpi = 1200)

# Summary table
qc_summary <- qc_df %>%
  group_by(subtype) %>%
  summarise(
    n_samples = n(),
    median_depth_M = median(depth_million, na.rm = TRUE),
    median_circ_per10M = median(circRNAs_per_10M, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(qc_summary, here("outputs", "tables", "QC_summary_depth_and_burden_per_subtype.csv"), row.names = FALSE)

# Correlations
cor_all <- cor(qc_df$depth_million, qc_df$detected_circs, use = "complete.obs", method = "spearman")
cor_sub <- qc_df %>%
  group_by(subtype) %>%
  summarise(rho = cor(depth_million, detected_circs, use = "complete.obs", method = "spearman"),
            .groups = "drop")

write.csv(cor_sub, here("outputs", "tables", "QC_correlation_depth_vs_circs_by_subtype.csv"), row.names = FALSE)

# PCA/MDS (TMM)
lib_sizes <- colSums(counts_all)
keep_samples <- lib_sizes > 0
counts_all_f <- counts_all[, keep_samples, drop = FALSE]
meta_f <- meta[match(colnames(counts_all_f), meta$sample), ]

dge <- DGEList(counts = counts_all_f)
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

pca <- prcomp(t(logCPM))
pca_df <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2]) %>%
  left_join(meta_f, by = "sample")

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = subtype, shape = dataset)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = subtype_colors) +
  labs(
    title = "PCA of TMM-normalized circRNA counts (QC only)",
    x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")
  ) +
  theme_classic()

ggsave(here("outputs", "figures", "QC_PCA_TMM_logCPM_circRNA.pdf"), p_pca, width = 6, height = 5, dpi = 1200)

pdf(here("outputs", "figures", "QC_MDS_TMM_logCPM_circRNA.pdf"), width = 6, height = 5)
plotMDS(dge, labels = meta_f$subtype, col = as.numeric(meta_f$subtype))
dev.off()

cat("Stage3 QC completed.\n")

# ============================================================
# Stage 4 — Final BED generation (hg38) from Stage2 SRPBM-filtered sets
# Inputs :
#   data/processed/stage02/<GSE>_filtered_counts_ge2_SRPBM_ge1.csv
# Outputs:
#   data/processed/stage04/<GSE>_final_hg38.bed
#   data/processed/stage04/<GSE>_final_hg38_with_header.tsv
#   outputs/tables/stage04_BED_summary_hg38_only.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(here)
})

dir.create(here("data", "processed", "stage04"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)

paths_final <- list(
  GSE142731 = here("data", "processed", "stage02", "GSE142731_filtered_counts_ge2_SRPBM_ge1.csv"),
  GSE191230 = here("data", "processed", "stage02", "GSE191230_filtered_counts_ge2_SRPBM_ge1.csv"),
  GSE207248 = here("data", "processed", "stage02", "GSE207248_filtered_counts_ge2_SRPBM_ge1.csv")
)

subtype_names <- c(GSE142731 = "TNBC", GSE191230 = "HER2+", GSE207248 = "IBC")

quick_file_check <- function(fp, label = "") {
  if (!file.exists(fp)) stop("File not found: ", fp)
  sz <- file.info(fp)$size
  if (is.na(sz) || sz < 5) stop("File looks empty: ", fp)
  cat("  - OK:", label, "|", fp, "| size =", sz, "bytes\n")
}

parse_circ_id_to_coords <- function(df) {
  if (!"circ_id" %in% colnames(df)) colnames(df)[1] <- "circ_id"
  
  df %>%
    mutate(
      circ_id = gsub("\\|", "-", as.character(circ_id)),
      chr = sub(":.*", "", circ_id),
      pos = sub(".*:", "", circ_id),
      start = suppressWarnings(as.integer(sub("-.*", "", pos))),
      end = suppressWarnings(as.integer(sub(".*-", "", pos)))
    ) %>%
    select(circ_id, chr, start, end)
}

bed_summary_list <- list()

for (label in names(paths_final)) {
  cat("\n=====================================\n")
  cat("Stage4 | Processing:", label, "(", subtype_names[[label]], ")\n")
  cat("=====================================\n")
  
  fp <- paths_final[[label]]
  quick_file_check(fp, "Stage2 output")
  
  df <- read.csv(fp, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) stop("Input has 0 rows: ", fp)
  
  coords <- parse_circ_id_to_coords(df)
  
  bad <- is.na(coords$chr) | is.na(coords$start) | is.na(coords$end)
  if (any(bad)) {
    cat(" Dropping NA coords:", sum(bad), "\n")
    coords <- coords[!bad, , drop = FALSE]
  }
  if (nrow(coords) == 0) stop("No valid coords after parsing circ_id for ", label)
  
  bed_df <- coords %>%
    mutate(
      start2 = pmin(start, end),
      end2 = pmax(start, end)
    ) %>%
    transmute(
      chr = chr,
      start = start2,
      end = end2,
      circ_id = paste0(chr, ":", start2, "-", end2),
      strand = ".",
      Subtype = subtype_names[[label]]
    ) %>%
    distinct(chr, start, end, .keep_all = TRUE)
  
  out_bed <- here("data", "processed", "stage04", paste0(label, "_final_hg38.bed"))
  out_tsv <- here("data", "processed", "stage04", paste0(label, "_final_hg38_with_header.tsv"))
  
  write_tsv(bed_df %>% select(chr, start, end, circ_id, strand), out_bed, col_names = FALSE)
  write_tsv(bed_df, out_tsv)
  
  cat("Saved:", out_bed, "\n")
  
  bed_summary_list[[label]] <- tibble(
    GSE = label,
    Subtype = subtype_names[[label]],
    n_stage2_rows = nrow(df),
    n_in_BED = nrow(bed_df),
    Genome = "hg38"
  )
}

bed_summary <- bind_rows(bed_summary_list)
out_summary <- here("outputs", "tables", "stage04_BED_summary_hg38_only.csv")
write.csv(bed_summary, out_summary, row.names = FALSE)

cat("\n===== Stage4 BED summary =====\n")
print(bed_summary)
cat(" Saved:", out_summary, "\n")
cat("Stage4 done.\n")

# ============================================================
# Stage 5 — Compare GSE BEDs (TNBC / HER2+ / IBC) vs databases
# Databases: circBase(hg38) + circAtlas(hg38) + MiOncoCirc(hg38)
#
# Inputs:
#   - data/processed/stage04/<GSE>_final_hg38.bed
#   - data/reference/circBase_Homo_sapiens_hg38_lifted.bed
#   - data/reference/MiOncoCirc_all_hg38.bed
#   - data/reference/CircAtlas_human_bed_v3.0.zip
#
# Outputs:
#   - outputs/tables/stage05_overlap_summary_OVERLAP.csv
#   - data/processed/stage05/<GSE>_overlap_DIAGNOSTIC_exact_vs_overlap.csv
#   - data/processed/stage05/<GSE>_NOVEL_circRNAs.csv
#   - data/processed/stage05/NOVEL_circRNAs_all.csv
#   - data/processed/stage05/NOVEL_circRNAs_canonical.csv
#   - data/processed/stage05/NOVEL_circRNAs_altcontigs.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(GenomicRanges)
  library(IRanges)
  library(here)
})

# -----------------------------
# 0) Config
# -----------------------------
RUN_SANITY_CHECKS <- TRUE
MIN_BP_OVERLAP <- 10  # increase to 50 for stricter overlaps

dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "processed", "stage05"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Inputs
# -----------------------------
paths_gse_bed <- list(
  GSE142731 = here("data", "processed", "stage04", "GSE142731_final_hg38.bed"),
  GSE191230 = here("data", "processed", "stage04", "GSE191230_final_hg38.bed"),
  GSE207248 = here("data", "processed", "stage04", "GSE207248_final_hg38.bed")
)

subtype_names <- c(GSE142731 = "TNBC", GSE191230 = "HER2+", GSE207248 = "IBC")

# Put reference BEDs/zips here in your repo
circbase_bed_path <- here("data", "reference", "circBase_Homo_sapiens_hg38_lifted.bed")
mion_bed_path     <- here("data", "reference", "MiOncoCirc_all_hg38.bed")

circatlas_zip <- here("data", "reference", "CircAtlas_human_bed_v3.0.zip")
circatlas_dir <- here("data", "reference", "CircAtlas_unzip")
circatlas_txt <- file.path(circatlas_dir, "human_bed_v3.0.txt")

# -----------------------------
# 2) Helpers
# -----------------------------
quick_file_check <- function(fp, label = "") {
  if (!file.exists(fp)) stop("File not found: ", fp)
  sz <- file.info(fp)$size
  if (is.na(sz) || sz < 5) stop("File exists but looks empty: ", fp)
  cat("  - OK:", label, "|", fp, "| size =", sz, "bytes\n")
}

clean_strand <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "*"
  ifelse(x %in% c("+", "-", "*"), x, "*")
}

read_bed5 <- function(fp) {
  df <- read_tsv(fp, col_names = FALSE, show_col_types = FALSE)
  if (ncol(df) < 3) stop("BED has <3 columns. Corrupted? ", fp)
  
  if (ncol(df) >= 5) {
    df <- df[, 1:5]
    colnames(df) <- c("chr", "start", "end", "circ_id", "strand")
  } else {
    colnames(df)[1:3] <- c("chr", "start", "end")
    if (ncol(df) == 3) {
      df$circ_id <- paste0(df$chr, ":", df$start, "-", df$end)
      df$strand <- "*"
    } else if (ncol(df) == 4) {
      colnames(df)[4] <- "circ_id"
      df$strand <- "*"
    }
    df <- df %>% select(chr, start, end, circ_id, strand)
  }
  
  df %>%
    mutate(
      chr   = as.character(chr),
      start = suppressWarnings(as.integer(start)),
      end   = suppressWarnings(as.integer(end)),
      start0 = start,
      end0   = end,
      start = pmin(start0, end0),
      end   = pmax(start0, end0),
      strand = clean_strand(strand)
    ) %>%
    select(chr, start, end, circ_id, strand) %>%
    filter(!is.na(chr), !is.na(start), !is.na(end)) %>%
    distinct(chr, start, end, .keep_all = TRUE)
}

make_keys <- function(df) paste(df$chr, df$start, df$end)

df_to_gr <- function(df) {
  GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start, end = df$end),
    strand   = clean_strand(df$strand),
    circ_id  = df$circ_id
  )
}

inspect_bed <- function(fp, label) {
  quick_file_check(fp, label)
  df <- read_bed5(fp)
  width <- df$end - df$start
  cat("\n---", label, "---\n")
  cat("n =", nrow(df), "\n")
  cat("start min/median/max:", min(df$start), median(df$start), max(df$start), "\n")
  cat("end   min/median/max:", min(df$end), median(df$end), max(df$end), "\n")
  cat("width min/median/max:", min(width), median(width), max(width), "\n")
  cat("Any start==0 ?", any(df$start == 0), "\n")
  cat("First 5 rows:\n")
  print(head(df, 5))
  invisible(df)
}

overlap_flag <- function(query_gr, subject_gr, min_bp = 10) {
  hits <- findOverlaps(query_gr, subject_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(rep(FALSE, length(query_gr)))
  
  ints <- pintersect(
    ranges(query_gr[queryHits(hits)]),
    ranges(subject_gr[subjectHits(hits)])
  )
  bp <- width(ints)
  
  keep_hit <- bp >= min_bp
  qhit <- queryHits(hits)[keep_hit]
  
  flag <- rep(FALSE, length(query_gr))
  flag[unique(qhit)] <- TRUE
  flag
}

# -----------------------------
# 3) Load databases (hg38)
# -----------------------------
cat("\n=== Loading databases (hg38) ===\n")

quick_file_check(circbase_bed_path, "circBase hg38")
cb_df <- read_bed5(circbase_bed_path)
cb_gr <- df_to_gr(cb_df)
cb_keys <- make_keys(cb_df)
cat("circBase unique:", nrow(cb_df), "\n")

quick_file_check(mion_bed_path, "MiOncoCirc hg38")
mion_df <- read_bed5(mion_bed_path)
mion_gr <- df_to_gr(mion_df)
mion_keys <- make_keys(mion_df)
cat("MiOncoCirc unique:", nrow(mion_df), "\n")

quick_file_check(circatlas_zip, "circAtlas zip")
if (!dir.exists(circatlas_dir)) dir.create(circatlas_dir, recursive = TRUE)
unzip(circatlas_zip, exdir = circatlas_dir)

quick_file_check(circatlas_txt, "circAtlas txt")
ca_raw <- read_tsv(circatlas_txt, col_names = TRUE, show_col_types = FALSE)
stopifnot(all(c("Chro", "Start", "End") %in% colnames(ca_raw)))

ca_df <- ca_raw %>%
  mutate(
    Chro  = as.character(Chro),
    Start = suppressWarnings(as.integer(Start)),
    End   = suppressWarnings(as.integer(End))
  ) %>%
  filter(!is.na(Chro), !is.na(Start), !is.na(End)) %>%
  transmute(
    chr = Chro,
    start0 = Start,
    end0   = End,
    start  = pmin(start0, end0),
    end    = pmax(start0, end0),
    circ_id = paste0(chr, ":", start, "-", end),
    strand  = "*"
  ) %>%
  distinct(chr, start, end, .keep_all = TRUE)

ca_gr <- df_to_gr(ca_df)
ca_keys <- make_keys(ca_df)
cat("circAtlas unique:", nrow(ca_df), "\n")

cat("\n Databases loaded successfully.\n")

# -----------------------------
# 4) Optional sanity checks
# -----------------------------
if (isTRUE(RUN_SANITY_CHECKS)) {
  cat("\n=== Sanity checks (coordinate conventions) ===\n")
  gse_df  <- inspect_bed(paths_gse_bed$GSE142731, "GSE142731 BED")
  cb_df2  <- inspect_bed(circbase_bed_path, "circBase hg38")
  mion_df2 <- inspect_bed(mion_bed_path, "MiOncoCirc hg38")
  
  cat("\nExact-match baseline (GSE142731 vs circBase):",
      sum(make_keys(gse_df) %in% make_keys(cb_df2)), "\n")
  cat("Exact-match baseline (GSE142731 vs MiOncoCirc):",
      sum(make_keys(gse_df) %in% make_keys(mion_df2)), "\n")
  
  gse_gr0  <- df_to_gr(gse_df)
  cb_gr0   <- df_to_gr(cb_df2)
  mion_gr0 <- df_to_gr(mion_df2)
  
  cat("Overlap (GSE142731 vs circBase):",
      length(unique(queryHits(findOverlaps(gse_gr0, cb_gr0, ignore.strand = TRUE)))),
      "of", length(gse_gr0), "\n")
  
  cat("Overlap (GSE142731 vs MiOncoCirc):",
      length(unique(queryHits(findOverlaps(gse_gr0, mion_gr0, ignore.strand = TRUE)))),
      "of", length(gse_gr0), "\n")
}

# -----------------------------
# 5) Core compare per GSE
# -----------------------------
compare_one_gse <- function(label) {
  cat("\n=====================================\n")
  cat("Stage5 | Processing:", label, "(", subtype_names[[label]], ")\n")
  cat("=====================================\n")
  
  bed_path <- paths_gse_bed[[label]]
  quick_file_check(bed_path, paste0(label, " BED"))
  
  g_df <- read_bed5(bed_path)
  if (nrow(g_df) == 0) stop("GSE BED empty after cleaning: ", label)
  
  g_keys <- make_keys(g_df)
  
  # Exact diagnostics
  in_cb_exact   <- g_keys %in% cb_keys
  in_ca_exact   <- g_keys %in% ca_keys
  in_mion_exact <- g_keys %in% mion_keys
  
  # Start-1 exact diagnostics (0/1-based hint)
  g_keys_m1 <- paste(g_df$chr, g_df$start - 1L, g_df$end)
  in_cb_exact_m1   <- g_keys_m1 %in% cb_keys
  in_ca_exact_m1   <- g_keys_m1 %in% ca_keys
  in_mion_exact_m1 <- g_keys_m1 %in% mion_keys
  
  # Overlap-based (main)
  g_gr <- df_to_gr(g_df)
  in_cb_ov   <- overlap_flag(g_gr, cb_gr,   min_bp = MIN_BP_OVERLAP)
  in_ca_ov   <- overlap_flag(g_gr, ca_gr,   min_bp = MIN_BP_OVERLAP)
  in_mion_ov <- overlap_flag(g_gr, mion_gr, min_bp = MIN_BP_OVERLAP)
  
  annot_df <- g_df %>%
    mutate(
      GSE = label,
      Subtype = subtype_names[[label]],
      
      in_circBase_exact   = in_cb_exact,
      in_circAtlas_exact  = in_ca_exact,
      in_MiOncoCirc_exact = in_mion_exact,
      
      in_circBase_exact_m1   = in_cb_exact_m1,
      in_circAtlas_exact_m1  = in_ca_exact_m1,
      in_MiOncoCirc_exact_m1 = in_mion_exact_m1,
      
      # MAIN overlap annotations
      in_circBase   = in_cb_ov,
      in_circAtlas  = in_ca_ov,
      in_MiOncoCirc = in_mion_ov,
      
      n_databases = as.integer(in_circBase) + as.integer(in_circAtlas) + as.integer(in_MiOncoCirc),
      novel_all = (n_databases == 0)
    )
  
  # Save per-GSE outputs
  out_diag <- here("data", "processed", "stage05", paste0(label, "_overlap_DIAGNOSTIC_exact_vs_overlap.csv"))
  write.csv(annot_df, out_diag, row.names = FALSE)
  
  novel_df <- annot_df %>%
    filter(novel_all) %>%
    transmute(
      GSE, Subtype, chr, start, end, circ_id, strand,
      n_databases, in_circBase, in_circAtlas, in_MiOncoCirc
    ) %>%
    arrange(chr, start, end)
  
  out_novel <- here("data", "processed", "stage05", paste0(label, "_NOVEL_circRNAs.csv"))
  write.csv(novel_df, out_novel, row.names = FALSE)
  
  cat("Saved diagnostic:", out_diag, "\n")
  cat("Saved NOVEL list:", out_novel, "| n =", nrow(novel_df), "\n")
  
  summary_row <- tibble(
    GSE = label,
    Subtype = subtype_names[[label]],
    total_circRNAs = nrow(annot_df),
    
    # MAIN overlap counts
    in_circBase   = sum(annot_df$in_circBase),
    in_circAtlas  = sum(annot_df$in_circAtlas),
    in_MiOncoCirc = sum(annot_df$in_MiOncoCirc),
    novel_all     = sum(annot_df$novel_all),
    
    # Diagnostics
    exact_circBase = sum(annot_df$in_circBase_exact),
    exact_MiOncoCirc = sum(annot_df$in_MiOncoCirc_exact),
    exactm1_circBase = sum(annot_df$in_circBase_exact_m1),
    exactm1_MiOncoCirc = sum(annot_df$in_MiOncoCirc_exact_m1)
  )
  
  list(summary = summary_row, novel = novel_df)
}

res_list <- lapply(names(paths_gse_bed), compare_one_gse)
overlap_summary <- bind_rows(lapply(res_list, `[[`, "summary"))

out_sum <- here("outputs", "tables", "stage05_overlap_summary_three_databases_hg38_only_OVERLAP.csv")
write.csv(overlap_summary, out_sum, row.names = FALSE)

cat("\n===== Stage5 Overall overlap summary (OVERLAP-based; hg38 only) =====\n")
print(overlap_summary)
cat(" Saved:", out_sum, "\n")

# -----------------------------
# 6) Novel classification (canonical vs alt contigs)
# -----------------------------
novel_all <- bind_rows(lapply(res_list, `[[`, "novel"))

out_all <- here("data", "processed", "stage05", "NOVEL_circRNAs_all.csv")
write_csv(novel_all, out_all)

canonical_chr <- c(paste0("chr", 1:22), "chrX", "chrY")

novel2 <- novel_all %>%
  mutate(
    chr = as.character(chr),
    Novel_class = ifelse(chr %in% canonical_chr, "novel_canonical", "novel_alt")
  )

novel_canonical <- novel2 %>%
  filter(Novel_class == "novel_canonical") %>%
  arrange(chr, start, end)

novel_alt <- novel2 %>%
  filter(Novel_class == "novel_alt") %>%
  arrange(chr, start, end)

out_can <- here("data", "processed", "stage05", "NOVEL_circRNAs_canonical.csv")
out_alt <- here("data", "processed", "stage05", "NOVEL_circRNAs_altcontigs.csv")

write_csv(novel_canonical, out_can)
write_csv(novel_alt, out_alt)

cat("\n===== Stage5 Novel circRNA classification summary =====\n")
cat("Total novel:", nrow(novel2), "\n")
cat("Canonical (chr1-22,X,Y):", nrow(novel_canonical), "\n")
cat("Alt contigs / non-canonical:", nrow(novel_alt), "\n\n")

cat("Canonical breakdown by subtype:\n")
print(novel_canonical %>% count(Subtype, name = "n_canonical"))

cat("\nAlt-contig breakdown by subtype:\n")
print(novel_alt %>% count(Subtype, name = "n_alt"))

cat("\nSaved:\n",
    "- ", out_all, "\n",
    "- ", out_can, "\n",
    "- ", out_alt, "\n", sep = "")
cat(" Stage5 DONE.\n")

# ============================================================
# Stage 6 — Read support distribution + Donut chart (>=2 reads distribution)
# Inputs:
#   - outputs/tables/stage01_circRNA_counts_before_after_ge2_filter.csv
#   - outputs/tables/stage02_circRNA_counts_after_SRPBM_filter_hg38_only.csv
#   - outputs/tables/stage05_overlap_summary_three_databases_hg38_only_OVERLAP.csv
#   - data/processed/stage01/<GSE>_filtered_counts_ge2.csv
# Outputs:
#   - outputs/tables/Table2_circRNA_summary_hg38_only_OVERLAP.csv
#   - outputs/figures/Fig3a_read_support_distribution_hg38_only.png
#   - outputs/figures/Fig3b_circRNAs_ge2_donut_hg38_only.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(here)
})

dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Config
# -----------------------------
gse_ids <- c("GSE142731", "GSE191230", "GSE207248")

subtype_names <- c(
  GSE142731 = "TNBC",
  GSE191230 = "HER2+",
  GSE207248 = "IBC"
)

pal <- c(
  "HER2+" = "#868686FF",
  "IBC"   = "#CD534CFF",
  "TNBC"  = "#0073C2FF"
)

# -----------------------------
# A) Table 2
# -----------------------------
counts_ge2_path <- here("outputs", "tables", "stage01_circRNA_counts_before_after_ge2_filter.csv")
counts_srp_path <- here("outputs", "tables", "stage02_circRNA_counts_after_SRPBM_filter_hg38_only.csv")
overlap_path    <- here("outputs", "tables", "stage05_overlap_summary_three_databases_hg38_only_OVERLAP.csv")

if (!file.exists(counts_ge2_path)) stop("Missing input: ", counts_ge2_path)
if (!file.exists(counts_srp_path)) stop("Missing input: ", counts_srp_path)
if (!file.exists(overlap_path)) stop("Missing input: ", overlap_path)

counts_ge2 <- read_csv(counts_ge2_path, show_col_types = FALSE) %>%
  filter(GSE %in% gse_ids)

counts_srp <- read_csv(counts_srp_path, show_col_types = FALSE) %>%
  filter(GSE %in% gse_ids)

hc_candidates <- c("n_after_present_ge2", "n_final_both", "n_final", "High_confidence")
hc_col <- hc_candidates[hc_candidates %in% colnames(counts_srp)][1]
if (is.na(hc_col)) {
  stop(
    "Could not find a high-confidence column in Stage02 summary.\n",
    "Present columns:\n- ", paste(colnames(counts_srp), collapse = "\n- "), "\n",
    "Expected one of: ", paste(hc_candidates, collapse = ", ")
  )
}

overlap_summary <- read_csv(overlap_path, show_col_types = FALSE) %>%
  filter(GSE %in% gse_ids)

tab2 <- counts_ge2 %>%
  inner_join(
    counts_srp %>% select(GSE, Subtype, High_confidence = all_of(hc_col)),
    by = "GSE"
  ) %>%
  inner_join(
    overlap_summary %>% select(GSE, total_circRNAs, novel_all),
    by = "GSE"
  ) %>%
  mutate(
    Known_any_db  = total_circRNAs - novel_all,
    Percent_novel = round(100 * novel_all / total_circRNAs, 2)
  ) %>%
  select(
    Subtype, GSE,
    Total_raw, Total_ge2,
    High_confidence,
    total_circRNAs,
    Known_any_db,
    novel_all,
    Percent_novel
  ) %>%
  arrange(match(Subtype, c("TNBC", "HER2+", "IBC")))

tab2_all <- tab2 %>%
  summarise(
    Subtype = "All subtypes",
    GSE = "All",
    Total_raw = sum(Total_raw, na.rm = TRUE),
    Total_ge2 = sum(Total_ge2, na.rm = TRUE),
    High_confidence = sum(High_confidence, na.rm = TRUE),
    total_circRNAs = sum(total_circRNAs, na.rm = TRUE),
    Known_any_db = sum(Known_any_db, na.rm = TRUE),
    novel_all = sum(novel_all, na.rm = TRUE)
  ) %>%
  mutate(Percent_novel = round(100 * novel_all / total_circRNAs, 2))

Table2_final <- bind_rows(tab2, tab2_all)

out_table2 <- here("outputs", "tables", "Table2_circRNA_summary_hg38_only_OVERLAP.csv")
write.csv(Table2_final, out_table2, row.names = FALSE)

cat("Saved:", out_table2, "\n")

# -----------------------------
# B) Read support distribution
# -----------------------------
file_paths_fig3a <- list(
  "TNBC"  = here("data", "processed", "stage01", "GSE142731_filtered_counts_ge2.csv"),
  "HER2+" = here("data", "processed", "stage01", "GSE191230_filtered_counts_ge2.csv"),
  "IBC"   = here("data", "processed", "stage01", "GSE207248_filtered_counts_ge2.csv")
)

processed_list <- lapply(names(file_paths_fig3a), function(subtype) {
  fp <- file_paths_fig3a[[subtype]]
  if (!file.exists(fp)) stop("Missing Fig3a input: ", fp)
  
  df <- read_csv(fp, show_col_types = FALSE)
  if (!"circ_id" %in% colnames(df)) colnames(df)[1] <- "circ_id"
  
  numeric_cols <- setdiff(colnames(df), "circ_id")
  df[, numeric_cols] <- lapply(df[, numeric_cols, drop = FALSE],
                               function(x) suppressWarnings(as.numeric(x)))
  
  df$Total_Reads <- rowSums(df[, numeric_cols, drop = FALSE], na.rm = TRUE)
  
  freq_df <- as.data.frame(table(df$Total_Reads), stringsAsFactors = FALSE)
  colnames(freq_df) <- c("Number_of_Reads", "Total_circRNAs")
  freq_df$Number_of_Reads <- as.integer(freq_df$Number_of_Reads)
  freq_df$Total_circRNAs  <- as.integer(freq_df$Total_circRNAs)
  freq_df$Subtype <- subtype
  freq_df
})

final_data <- bind_rows(processed_list) %>%
  mutate(Reads_Grouped = ifelse(Number_of_Reads > 150, 150, Number_of_Reads)) %>%
  group_by(Subtype, Reads_Grouped) %>%
  summarise(Total_circRNAs = sum(Total_circRNAs), .groups = "drop") %>%
  filter(Reads_Grouped >= 2)

final_data$Subtype <- factor(final_data$Subtype, levels = c("TNBC", "HER2+", "IBC"))

p_fig3a <- ggplot(final_data, aes(x = factor(Reads_Grouped), y = Total_circRNAs, fill = Subtype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  scale_x_discrete(breaks = as.character(seq(2, 150, by = 10))) +
  labs(
    x = "Number of supporting reads (summed per circRNA)",
    y = "Number of circRNAs (log10 scale)",
    fill = "Subtype"
  ) +
  scale_fill_manual(values = pal) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    panel.border = element_blank()
  )

out_fig3a <- here("outputs", "figures", "Fig3a_read_support_distribution_hg38_only.png")
ggsave(out_fig3a, plot = p_fig3a, width = 8, height = 4, dpi = 1200)
cat("Saved:", out_fig3a, "\n")

# -----------------------------
# C) Donut chart (>=2 reads distribution)
# -----------------------------
counts_df <- read_csv(counts_ge2_path, show_col_types = FALSE) %>%
  filter(GSE %in% gse_ids) %>%
  mutate(Subtype = subtype_names[GSE])

donut_df <- counts_df %>%
  transmute(GSE, Subtype, Count = Total_ge2) %>%
  mutate(
    Subtype = factor(Subtype, levels = c("IBC", "TNBC", "HER2+")),
    Percent = Count / sum(Count) * 100
  )

total_all <- sum(donut_df$Count)

p_fig3b <- ggplot(donut_df, aes(x = 2, y = Count, fill = Subtype)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = pal) +
  geom_text(
    aes(label = paste0(round(Percent, 2), "%")),
    position = position_stack(vjust = 0.5),
    size = 3.8,
    color = "black"
  ) +
  annotate("text", x = 0.5, y = 0, label = paste0("Total = ", total_all),
           size = 5, fontface = "bold") +
  labs(
    title = "Distribution of circRNAs (>=2 reads) in breast cancer subtypes",
    fill  = "Subtype"
  ) +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

out_fig3b <- here("outputs", "figures", "Fig3b_circRNAs_ge2_donut_hg38_only.png")
ggsave(out_fig3b, plot = p_fig3b, width = 7, height = 4.5, dpi = 1200)
cat("Saved:", out_fig3b, "\n")

cat("Stage 6 complete.\n")

# ============================================================
# Stage 7 — Exon structure + Chi-square test
# Inputs:
#   - data/processed/stage01/<GSE>_filtered_counts_ge2.csv
#   - data/reference/genes_hg38.gtf
# Outputs:
#   - outputs/tables/Exon_structure_summary_by_subtype.csv
#   - outputs/tables/circRNA_exon_counts_per_circRNA.csv
#   - outputs/figures/Fig3c_ExonCountDistribution_AllSubtypes.png
#   - outputs/tables/ExonCount_Groups_by_Subtype_for_ChiSquare.csv
#   - outputs/tables/ChiSquare_exon_groups_results.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(IRanges)
  library(tibble)
  library(here)
})

dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)

counts_files <- list(
  "TNBC"  = here("data", "processed", "stage01", "GSE142731_filtered_counts_ge2.csv"),
  "HER2+" = here("data", "processed", "stage01", "GSE191230_filtered_counts_ge2.csv"),
  "IBC"   = here("data", "processed", "stage01", "GSE207248_filtered_counts_ge2.csv")
)

gtf_hg38 <- here("data", "reference", "genes_hg38.gtf")
if (!file.exists(gtf_hg38)) stop("Missing GTF reference: ", gtf_hg38)

parse_circ_file <- function(file_path) {
  if (!file.exists(file_path)) stop("Counts file not found: ", file_path)
  df <- read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"circ_id" %in% colnames(df)) colnames(df)[1] <- "circ_id"
  
  df %>%
    mutate(
      circ_id = gsub("\\|", "-", as.character(circ_id)),
      chr     = sub(":.*", "", circ_id),
      pos     = sub(".*:", "", circ_id),
      start   = suppressWarnings(as.integer(sub("-.*", "", pos))),
      end     = suppressWarnings(as.integer(sub(".*-", "", pos)))
    ) %>%
    select(circ_id, chr, start, end) %>%
    filter(!is.na(start), !is.na(end))
}

cat("Building TxDb from GTF...\n")
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_hg38, format = "gtf")
exons_gr <- unique(GenomicFeatures::exons(txdb))
cat("Total unique exons:", length(exons_gr), "\n")

compute_exon_counts <- function(subtype, counts_file) {
  circ_df <- parse_circ_file(counts_file)
  
  circ_gr <- GRanges(
    seqnames = circ_df$chr,
    ranges   = IRanges(start = circ_df$start, end = circ_df$end),
    strand   = "*"
  )
  
  n_exons <- GenomicRanges::countOverlaps(circ_gr, exons_gr)
  
  detail <- circ_df %>%
    mutate(Subtype = subtype, n_exons = as.integer(n_exons))
  
  n_total    <- nrow(detail)
  n_exon1    <- sum(detail$n_exons == 1, na.rm = TRUE)
  n_exon_lt6 <- sum(detail$n_exons < 6,  na.rm = TRUE)
  
  summary <- tibble(
    Subtype      = subtype,
    n_circ_ge2   = n_total,
    n_exon1      = n_exon1,
    pct_exon1    = round(100 * n_exon1 / n_total, 2),
    n_exon_lt6   = n_exon_lt6,
    pct_exon_lt6 = round(100 * n_exon_lt6 / n_total, 2)
  )
  
  list(detail = detail, summary = summary)
}

res <- lapply(names(counts_files), function(st) compute_exon_counts(st, counts_files[[st]]))
names(res) <- names(counts_files)

exon_details <- bind_rows(lapply(res, `[[`, "detail"))
exon_summary <- bind_rows(lapply(res, `[[`, "summary"))

out_sum <- here("outputs", "tables", "Exon_structure_summary_by_subtype.csv")
out_det <- here("outputs", "tables", "circRNA_exon_counts_per_circRNA.csv")
write.csv(exon_summary, out_sum, row.names = FALSE)
write.csv(exon_details, out_det, row.names = FALSE)
cat("Saved:", out_sum, "\n")
cat("Saved:", out_det, "\n")

# Fig 3c (pooled)
exon_dist_all <- exon_details %>%
  mutate(Exon_group = ifelse(n_exons >= 10, ">=10", as.character(n_exons))) %>%
  count(Exon_group, name = "n_circRNAs") %>%
  mutate(Exon_group = factor(Exon_group, levels = c(as.character(0:9), ">=10"))) %>%
  filter(!is.na(Exon_group))

total_circ <- sum(exon_dist_all$n_circRNAs)
exon_dist_all <- exon_dist_all %>%
  mutate(Percent = 100 * n_circRNAs / total_circ)

max_count <- max(exon_dist_all$n_circRNAs)
max_pct   <- max(exon_dist_all$Percent)
scaleFactor <- max_count / max_pct

p_fig3c <- ggplot(exon_dist_all, aes(x = Exon_group)) +
  geom_col(aes(y = n_circRNAs), width = 0.7) +
  geom_line(aes(y = Percent * scaleFactor, group = 1), linetype = "dashed", linewidth = 1.2) +
  geom_point(aes(y = Percent * scaleFactor), size = 3) +
  scale_y_continuous(
    name = "Number of circRNAs",
    sec.axis = sec_axis(~ . / scaleFactor, name = "Percentage (%)")
  ) +
  labs(x = "Number of exons in circRNAs",
       title = "Exon count distribution in breast cancer circRNAs") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

out_fig <- here("outputs", "figures", "Fig3c_ExonCountDistribution_AllSubtypes.png")
ggsave(out_fig, plot = p_fig3c, width = 8, height = 4.5, dpi = 1200)
cat("Saved:", out_fig, "\n")

# Chi-square grouping
count_exon_groups <- function(df) {
  df %>%
    mutate(group = case_when(
      n_exons == 1 ~ "1_exon",
      n_exons >= 2 & n_exons <= 5 ~ "2_5_exons",
      n_exons >= 6 ~ "6plus_exons",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(group)) %>%
    count(group, name = "n")
}

exon_group_table <- exon_details %>%
  group_by(Subtype) %>%
  group_modify(~ count_exon_groups(.x)) %>%
  ungroup() %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)

out_group <- here("outputs", "tables", "ExonCount_Groups_by_Subtype_for_ChiSquare.csv")
write.csv(exon_group_table, out_group, row.names = FALSE)
cat("Saved:", out_group, "\n")

chisq_mat <- exon_group_table %>%
  select(1_exon, 2_5_exons, 6plus_exons) %>%
  as.matrix()
rownames(chisq_mat) <- exon_group_table$Subtype

chisq_res <- chisq.test(chisq_mat)

out_txt <- here("outputs", "tables", "ChiSquare_exon_groups_results.txt")
sink(out_txt)
cat("Chi-square test for exon group distributions\n\n")
print(chisq_res)
cat("\nExpected counts:\n")
print(chisq_res$expected)
cat("\nResiduals:\n")
print(chisq_res$residuals)
cat("\nRow-wise percentages (%):\n")
prop_table <- chisq_mat / rowSums(chisq_mat)
print(round(prop_table * 100, 2))
sink()

cat("Saved:", out_txt, "\n")
cat("Stage 7 complete.\n")

# ============================================================
# Stage 8 — Chromosomal SRPBM density
# Inputs:
#   - data/processed/stage02/<GSE>_filtered_counts_ge2_SRPBM_ge1.csv  (must contain mean_SRPBM)
#   - data/reference/genes_hg38.gtf
# Outputs:
#   - outputs/tables/Chromosomal_SRPBM_density_all_subtypes.csv
#   - outputs/tables/Table3_peak_chromosomal_density.csv
#   - outputs/tables/Chromosomal_density_summary_per_subtype.csv
#   - outputs/figures/Fig3d_chromosomal_SRPBM_density_bar.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(IRanges)
  library(here)
})

dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)

paths_final <- list(
  "TNBC"  = here("data", "processed", "stage02", "GSE142731_filtered_counts_ge2_SRPBM_ge1.csv"),
  "HER2+" = here("data", "processed", "stage02", "GSE191230_filtered_counts_ge2_SRPBM_ge1.csv"),
  "IBC"   = here("data", "processed", "stage02", "GSE207248_filtered_counts_ge2_SRPBM_ge1.csv")
)

subtype_order <- c("TNBC", "HER2+", "IBC")

gtf_hg38 <- here("data", "reference", "genes_hg38.gtf")
if (!file.exists(gtf_hg38)) stop("Missing GTF reference: ", gtf_hg38)

# Build TxDb and estimate chromosome lengths (Mb)
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_hg38, format = "gtf")
exons_df <- as.data.frame(GenomicFeatures::exons(txdb))

chrom_lengths <- exons_df %>%
  group_by(chr = as.character(seqnames)) %>%
  summarise(length_bp = max(end, na.rm = TRUE), .groups = "drop") %>%
  filter(chr %in% paste0("chr", c(1:22, "X"))) %>%
  mutate(
    length_Mb = length_bp / 1e6,
    Chr = sub("chr", "Chr", chr)
  ) %>%
  arrange(chr)

compute_chr_density <- function(subtype, final_file) {
  if (!file.exists(final_file)) stop("Missing Stage02 input: ", final_file)
  
  df <- read.csv(final_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"circ_id" %in% colnames(df)) colnames(df)[1] <- "circ_id"
  if (!"mean_SRPBM" %in% colnames(df)) stop("mean_SRPBM column missing in: ", final_file)
  
  df2 <- df %>%
    mutate(
      circ_id = gsub("\\|", "-", as.character(circ_id)),
      chr = sub(":.*", "", circ_id)
    )
  
  canonical <- paste0("chr", c(1:22, "X"))
  
  df2 %>%
    filter(chr %in% canonical) %>%
    group_by(chr) %>%
    summarise(
      Total_SRPBM = sum(mean_SRPBM, na.rm = TRUE),
      circRNA_count = n(),
      .groups = "drop"
    ) %>%
    mutate(Subtype = subtype)
}

chrom_density_all <- bind_rows(lapply(names(paths_final), function(st) {
  compute_chr_density(st, paths_final[[st]])
}))

chrom_density_all <- chrom_density_all %>%
  left_join(chrom_lengths, by = "chr") %>%
  filter(!is.na(length_Mb)) %>%
  mutate(
    Density = Total_SRPBM / length_Mb,
    Chr = factor(Chr, levels = paste0("Chr", c(1:22, "X"))),
    Subtype = factor(Subtype, levels = subtype_order)
  )

out_all <- here("outputs", "tables", "Chromosomal_SRPBM_density_all_subtypes.csv")
write.csv(chrom_density_all, out_all, row.names = FALSE)
cat("Saved:", out_all, "\n")

Table3 <- chrom_density_all %>%
  group_by(Subtype) %>%
  slice_max(order_by = Density, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    Group = Subtype,
    Chr,
    Total_SRPBM,
    length_Mb,
    Density,
    circRNA_count
  )

out_t3 <- here("outputs", "tables", "Table3_peak_chromosomal_density.csv")
write.csv(Table3, out_t3, row.names = FALSE)
cat("Saved:", out_t3, "\n")

density_summary <- chrom_density_all %>%
  group_by(Subtype) %>%
  summarise(
    mean_density = mean(Density, na.rm = TRUE),
    median_density = median(Density, na.rm = TRUE),
    .groups = "drop"
  )

out_sum <- here("outputs", "tables", "Chromosomal_density_summary_per_subtype.csv")
write.csv(density_summary, out_sum, row.names = FALSE)
cat("Saved:", out_sum, "\n")

pal <- c("HER2+" = "#868686FF", "IBC" = "#CD534CFF", "TNBC" = "#0073C2FF")

p_fig3d <- ggplot(chrom_density_all, aes(x = Chr, y = Density, fill = Subtype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = pal) +
  labs(
    x = "Chromosome",
    y = "circRNA SRPBM density (SRPBM / Mb)",
    fill = "Subtype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

out_fig <- here("outputs", "figures", "Fig3d_chromosomal_SRPBM_density_bar.png")
ggsave(out_fig, plot = p_fig3d, width = 9, height = 4.5, dpi = 1200)
cat("Saved:", out_fig, "\n")

cat("Stage 8 complete.\n")

# ============================================================
# Stage 9 — Venn diagram of shared high-confidence circRNAs
# Inputs:
#   - data/processed/stage04/GSE142731_final_hg38.bed
#   - data/processed/stage04/GSE191230_final_hg38.bed
#   - data/processed/stage04/GSE207248_final_hg38.bed
# Outputs:
#   - outputs/tables/Venn_set_sizes_high_confidence.csv
#   - outputs/figures/Fig3e_Venn_high_confidence_circRNAs.pdf
#   - outputs/tables/Shared_circRNAs_TNBC_HER2_IBC.csv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(VennDiagram)
  library(here)
})

dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Inputs (Stage4 BEDs)
# -----------------------------
bed_paths <- list(
  TNBC   = here("data", "processed", "stage04", "GSE142731_final_hg38.bed"),
  `HER2+`= here("data", "processed", "stage04", "GSE191230_final_hg38.bed"),
  IBC    = here("data", "processed", "stage04", "GSE207248_final_hg38.bed")
)

# consistent palette with other figures
pal <- c(
  "TNBC"  = "#0073C2FF",
  "HER2+" = "#868686FF",
  "IBC"   = "#CD534CFF"
)

read_circ_ids_from_bed <- function(path) {
  if (!file.exists(path)) stop("BED file not found: ", path)
  
  bed <- read_tsv(path, col_names = FALSE, show_col_types = FALSE)
  if (ncol(bed) < 4) stop("Unexpected BED format (need >=4 cols): ", path)
  
  circ_ids <- as.character(bed[[4]])
  circ_ids <- gsub("\\|", "-", circ_ids)
  unique(circ_ids)
}

# -----------------------------
# Build sets
# -----------------------------
set_list <- lapply(bed_paths, read_circ_ids_from_bed)

venn_sizes <- tibble(
  Subtype = names(set_list),
  n_circ  = sapply(set_list, length)
)

out_sizes <- here("outputs", "tables", "Venn_set_sizes_high_confidence.csv")
write.csv(venn_sizes, out_sizes, row.names = FALSE)
cat("Saved:", out_sizes, "\n")
print(venn_sizes)

# -----------------------------
# Draw Venn (3 sets)
# -----------------------------
venn_plot <- venn.diagram(
  x = set_list,
  filename = NULL,
  category.names = c("TNBC", "HER2+", "IBC"),
  fill = pal[c("TNBC", "HER2+", "IBC")],
  col = "white",
  lwd = 1.2,
  alpha = 0.7,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = "black",
  cat.pos = 0,
  cat.dist = 0.05,
  margin = 0.1
)

out_pdf <- here("outputs", "figures", "Fig3e_Venn_high_confidence_circRNAs.pdf")
pdf(out_pdf, width = 7, height = 7, useDingbats = FALSE)
grid::grid.newpage()
grid::grid.draw(venn_plot)
dev.off()

cat("Saved:", out_pdf, "\n")

# -----------------------------
# Shared circRNAs across all 3 subtypes
# -----------------------------
shared_all3 <- Reduce(intersect, set_list)

out_shared <- here("outputs", "tables", "Shared_circRNAs_TNBC_HER2_IBC.csv")
write.csv(data.frame(circ_id = shared_all3), out_shared, row.names = FALSE)

cat("Shared across TNBC + HER2+ + IBC:", length(shared_all3), "\n")
cat("Saved:", out_shared, "\n")

cat("Stage 9 complete.\n")


# ============================================================
# Stage 10 — SRPBM / FPKM integration for shared circRNAs
# Purpose:
#   - Merge SRPBM (circRNA abundance) with parental mRNA FPKM
#   - Prepare clean LONG and WIDE matrices for shared circRNAs
#   - Generate SRPBM heatmap and PCC analyses
# Notes on upstream inputs:
#   - *_combined_gene_abundance.tsv files are outputs of parental gene expression
#     quantification (FPKM) computed from TopHat2-aligned BAM files.
#   - *_circfinal_with_gene_names.csv files are derived from CIRIquant outputs
#     (circfinal.csv).
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
})

# ============================================================
# Configuration
# ============================================================
base_dir <- "F:/Me/Expresion/Articles/My_Data/Final_Words/R_Analysis/Redrawing_shapes/3.1"

stage5_dir <- file.path(base_dir, "Stage5_FPKM_HeatMap")
stage7_dir <- file.path(base_dir, "Stage7_Shared_Prepare")
stage8_dir <- file.path(base_dir, "Stage8_Heatmap_PCC")

dir.create(stage7_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stage8_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Input files
# ============================================================

shared_fp <- file.path(
  base_dir,
  "Stage4_Venn",
  "Shared_circRNAs_TNBC_HER2_IBC.csv"
)

merged_long_paths <- list(
  TNBC = file.path(stage5_dir, "GSE142731_SRPBM_FPKM_Merged_LONG.csv"),
  HER2 = file.path(stage5_dir, "GSE191230_SRPBM_FPKM_Merged_LONG.csv"),
  IBC  = file.path(stage5_dir, "GSE207248_SRPBM_FPKM_Merged_LONG.csv")
)

# ============================================================
# Helper functions
# ============================================================

normalize_circ_id <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\\|", "-") %>%
    str_replace_all("\u2013|\u2014", "-") %>%
    str_remove_all("\\s+") %>%
    str_trim()
}

# ============================================================
# Stage 8.1 — Build clean LONG files for shared circRNAs
# ============================================================

shared_df <- read.csv(shared_fp, check.names = FALSE)
circ_col <- if ("circ_id" %in% names(shared_df)) "circ_id" else names(shared_df)[1]

shared_circs <- shared_df[[circ_col]] %>%
  normalize_circ_id() %>%
  unique()

build_shared_long <- function(subtype, fp_long, shared_circs, out_dir) {
  
  df <- read.csv(fp_long, check.names = FALSE)
  
  stopifnot(all(c("circ_id", "sample", "SRPBM", "FPKM") %in% names(df)))
  
  df$circ_id <- normalize_circ_id(df$circ_id)
  
  out <- df %>%
    filter(circ_id %in% shared_circs) %>%
    mutate(
      Subtype = subtype,
      SRPBM = as.numeric(SRPBM),
      FPKM  = as.numeric(FPKM)
    ) %>%
    select(
      Subtype, circ_id,
      any_of(c("Gene.ID", "Gene.Name")),
      sample, SRPBM, FPKM
    ) %>%
    arrange(circ_id, sample)
  
  out_fp <- file.path(out_dir, paste0("SharedCircRNAs_", subtype, "_LONG.csv"))
  write.csv(out, out_fp, row.names = FALSE)
  
  out
}

shared_long_list <- lapply(names(merged_long_paths), function(st) {
  build_shared_long(st, merged_long_paths[[st]], shared_circs, stage7_dir)
})
names(shared_long_list) <- names(merged_long_paths)

shared_all <- bind_rows(shared_long_list)
write.csv(
  shared_all,
  file.path(stage7_dir, "SharedCircRNAs_ALL_LONG.csv"),
  row.names = FALSE
)

# ============================================================
# Stage 8.2 — Build WIDE matrices (SRPBM / FPKM)
# ============================================================

make_wide <- function(df, value_col, subtype, out_dir) {
  wide <- df %>%
    select(circ_id, sample, !!sym(value_col)) %>%
    pivot_wider(names_from = sample, values_from = !!sym(value_col))
  
  out_fp <- file.path(
    out_dir,
    paste0("SharedCircRNAs_", subtype, "_WIDE_", value_col, ".csv")
  )
  write.csv(wide, out_fp, row.names = FALSE)
}

for (st in names(shared_long_list)) {
  make_wide(shared_long_list[[st]], "SRPBM", st, stage7_dir)
  make_wide(shared_long_list[[st]], "FPKM",  st, stage7_dir)
}

# ============================================================
# Stage 11 — SRPBM heatmap of shared circRNAs
# ============================================================

read_wide <- function(fp) {
  df <- read.csv(fp, check.names = FALSE)
  stopifnot("circ_id" %in% names(df))
  df
}

TNBC_w <- read_wide(file.path(stage7_dir, "SharedCircRNAs_TNBC_WIDE_SRPBM.csv"))
HER2_w <- read_wide(file.path(stage7_dir, "SharedCircRNAs_HER2_WIDE_SRPBM.csv"))
IBC_w  <- read_wide(file.path(stage7_dir, "SharedCircRNAs_IBC_WIDE_SRPBM.csv"))

common_circ <- Reduce(intersect, list(TNBC_w$circ_id, HER2_w$circ_id, IBC_w$circ_id))
common_circ <- sort(common_circ)

align_rows <- function(df, circ_ids) {
  df %>%
    filter(circ_id %in% circ_ids) %>%
    arrange(match(circ_id, circ_ids))
}

TNBC_w <- align_rows(TNBC_w, common_circ)
HER2_w <- align_rows(HER2_w, common_circ)
IBC_w  <- align_rows(IBC_w,  common_circ)

srpbm_mat <- cbind(
  IBC_w[, -1, drop = FALSE],
  TNBC_w[, -1, drop = FALSE],
  HER2_w[, -1, drop = FALSE]
)

srpbm_mat <- apply(srpbm_mat, 2, as.numeric)
srpbm_mat <- as.matrix(srpbm_mat)

srpbm_log <- log2(srpbm_mat + 1)

annotation_col <- data.frame(
  Subtype = factor(
    c(
      rep("IBC",   ncol(IBC_w)  - 1),
      rep("TNBC",  ncol(TNBC_w) - 1),
      rep("HER2",  ncol(HER2_w) - 1)
    ),
    levels = c("IBC", "TNBC", "HER2")
  )
)
rownames(annotation_col) <- colnames(srpbm_log)

subtype_colors <- list(
  Subtype = c(
    TNBC = "#0073C2FF",
    HER2 = "#868686FF",
    IBC  = "#CD534CFF"
  )
)

pdf(
  file.path(stage8_dir, "Heatmap_SharedCircRNAs_SRPBM.pdf"),
  width = 18, height = 6, useDingbats = FALSE
)

pheatmap(
  srpbm_log,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  annotation_col = annotation_col,
  annotation_colors = subtype_colors,
  main = "Shared circRNAs — SRPBM heatmap (IBC | TNBC | HER2)"
)

dev.off()

# ============================================================
# Stage 12 — PCC analysis (SRPBM vs FPKM)
# ============================================================

df <- shared_all %>%
  mutate(
    SRPBM = as.numeric(SRPBM),
    FPKM  = as.numeric(FPKM)
  )

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(c(r = NA, p = NA))
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(c(r = NA, p = NA))
  ct <- cor.test(x[ok], y[ok], method = "pearson")
  c(r = unname(ct$estimate), p = ct$p.value)
}

pcc_tbl <- df %>%
  group_by(Subtype, circ_id) %>%
  summarise(
    Gene.ID = paste(unique(na.omit(Gene.ID)), collapse = ";"),
    n_pairs = sum(is.finite(SRPBM) & is.finite(FPKM)),
    PCC     = safe_cor(SRPBM, FPKM)["r"],
    p_value = safe_cor(SRPBM, FPKM)["p"],
    .groups = "drop"
  ) %>%
  group_by(Subtype) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(
  pcc_tbl,
  file.path(stage8_dir, "PCC_per_circRNA_per_Subtype.csv"),
  row.names = FALSE
)

global_pcc <- df %>%
  group_by(Subtype) %>%
  summarise(
    n_pairs = sum(is.finite(SRPBM) & is.finite(FPKM)),
    PCC_global = cor(SRPBM, FPKM, use = "complete.obs", method = "pearson"),
    .groups = "drop"
  )

write.csv(
  global_pcc,
  file.path(stage8_dir, "PCC_global_per_Subtype.csv"),
  row.names = FALSE
)

# ============================================================
# Stage 13 — circRNA–miRNA–mRNA network (Cytoscape-ready)
# Purpose:
#   - Filter a global circRNA–miRNA–mRNA interaction reference for target circRNAs
#   - Summarize key network statistics (top miRNAs, top mRNAs, binding energy)
#   - Export Cytoscape nodes/edges for:
#       (A) Full network (supplementary figure)
#       (B) Curated network (main figure)
#
# Assumptions about input reference:
#   - The reference file contains columns: circRNA, miRNA, mRNA, Energy, Support Type
#   - circRNA ids match the format: chr:start-end
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# ============================================================
# Configuration
# ============================================================
base_dir  <- "F:/Me/Expresion/Articles/My_Data/Final_Words/R_Analysis/Redrawing_shapes/3.1"
stage11_dir <- file.path(base_dir, "Stage11_circRNA_miRNA_mRNA_network")
dir.create(stage11_dir, recursive = TRUE, showWarnings = FALSE)

ref_file <- file.path(
  "F:/Me/Expresion/Articles/My_Data/Final_Words/R_Analysis/Mutual_Analysis",
  "circRNA-miRNA",
  "circRNA-miRNA-mRNA",
  "all",
  "All_circRNA_miRNA_mRNA.TSV"
)

# Output names (written into stage11_dir)
out_triplets_15    <- file.path(stage11_dir, "circRNA_miRNA_mRNA_15circRNAs.tsv")
out_top4_miRNA     <- file.path(stage11_dir, "top4_miRNAs_freq.tsv")
out_top4_mRNA      <- file.path(stage11_dir, "top4_mRNAs_freq.tsv")
out_top5_circEner  <- file.path(stage11_dir, "top5_circRNAs_most_negative_avg_energy.tsv")

out_nodes_full <- file.path(stage11_dir, "Cytoscape_nodes_FullNetwork.tsv")
out_edges_full <- file.path(stage11_dir, "Cytoscape_edges_FullNetwork.tsv")

out_nodes_main <- file.path(stage11_dir, "Cytoscape_nodes_MainFigure_core5_15miR_30mRNA.tsv")
out_edges_main <- file.path(stage11_dir, "Cytoscape_edges_MainFigure_core5_15miR_30mRNA.tsv")

# ============================================================
# Target circRNAs (15 shared circRNAs)
# ============================================================
target_circs <- c(
  "chr10:37153572-37188019",
  "chr10:37165094-37185712",
  "chr11:96091892-96093517",
  "chr17:82563354-82568201",
  "chr19:39877663-39893527",
  "chr19:39879555-39894495",
  "chr19:39883305-39896143",
  "chr19:39885673-39899697",
  "chr19:39885991-39902124",
  "chr19:39889479-39902837",
  "chr8:51860845-51861246",
  "chr8:98706467-98707311",
  "chr15:73702772-73703997",
  "chr6:4891713-4892379",
  "chrX:131749306-131794466"
)

# ============================================================
# Helpers
# ============================================================
normalize_circ_id <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\\|", "-") %>%
    str_replace_all("\u2013|\u2014", "-") %>%
    str_remove_all("\\s+") %>%
    str_trim()
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3) return(list(r = NA_real_, p = NA_real_, n = n))
  if (sd(x) == 0 || sd(y) == 0) return(list(r = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  list(r = unname(ct$estimate), p = ct$p.value, n = n)
}

# ============================================================
# 1) Read reference interactions
# ============================================================
if (!file.exists(ref_file)) stop("Reference interaction file not found: ", ref_file)

all_int <- read.delim(
  ref_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Stage11 | Input columns:\n")
print(colnames(all_int))

required_cols <- c("circRNA", "miRNA", "mRNA")
missing_cols <- setdiff(required_cols, colnames(all_int))
if (length(missing_cols) > 0) {
  stop("Missing required columns in reference file: ", paste(missing_cols, collapse = ", "))
}

all_int <- all_int %>%
  mutate(
    circRNA = normalize_circ_id(circRNA),
    miRNA   = as.character(miRNA),
    mRNA    = as.character(mRNA)
  )

# ============================================================
# 2) Filter interactions for target circRNAs
# ============================================================
sub_int <- all_int %>%
  filter(circRNA %in% normalize_circ_id(target_circs))

n_circ <- n_distinct(sub_int$circRNA)
cat("\nStage11 | Filtered interactions\n")
cat("  - Triplets rows:", nrow(sub_int), "\n")
cat("  - Unique circRNAs:", n_circ, "\n")

if (nrow(sub_int) == 0) {
  stop(
    "No interactions found for the 15 target circRNAs in the reference file.\n",
    "Check circRNA ids and formatting."
  )
}

# Save filtered triplets
write.table(
  sub_int,
  file = out_triplets_15,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("Saved:", out_triplets_15, "\n")

# ============================================================
# 3) Basic network statistics
# ============================================================
unique_miRNAs <- unique(sub_int$miRNA)
unique_mRNAs  <- unique(sub_int$mRNA)

cat("\nStage11 | Unique nodes (from filtered triplets)\n")
cat("  - Unique miRNAs:", length(unique_miRNAs), "\n")
cat("  - Unique mRNAs :", length(unique_mRNAs), "\n")

miRNA_counts <- sort(table(sub_int$miRNA), decreasing = TRUE)
top4_miRNA <- head(miRNA_counts, 4)

mRNA_counts <- sort(table(sub_int$mRNA), decreasing = TRUE)
top4_mRNA <- head(mRNA_counts, 4)

cat("\nTop 4 miRNAs (frequency):\n")
print(top4_miRNA)

cat("\nTop 4 mRNAs (frequency):\n")
print(top4_mRNA)

# ============================================================
# 4) Top 5 circRNAs by most negative mean binding energy
# ============================================================
if (!"Energy" %in% colnames(sub_int)) {
  warning("Column 'Energy' not found in input. Skipping energy-based ranking.")
  energy_by_circ <- tibble(circRNA = character(), avg_Energy = numeric())
} else {
  sub_int <- sub_int %>% mutate(Energy_num = safe_numeric(Energy))
  
  if (all(is.na(sub_int$Energy_num))) {
    stop("Energy column cannot be converted to numeric (all NA). Check Energy format.")
  }
  
  energy_by_circ <- sub_int %>%
    group_by(circRNA) %>%
    summarise(avg_Energy = mean(Energy_num, na.rm = TRUE), .groups = "drop") %>%
    arrange(avg_Energy)
  
  top5_circ_neg <- energy_by_circ %>% slice_head(n = 5)
  
  cat("\nTop 5 circRNAs with most negative mean Energy:\n")
  print(top5_circ_neg)
  
  write.table(top5_circ_neg, out_top5_circEner, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved:", out_top5_circEner, "\n")
}

# Save top4 tables
write.table(top4_miRNA, out_top4_miRNA, sep = "\t", quote = FALSE, col.names = NA)
write.table(top4_mRNA,  out_top4_mRNA,  sep = "\t", quote = FALSE, col.names = NA)
cat("Saved:", out_top4_miRNA, "\n")
cat("Saved:", out_top4_mRNA, "\n")

# ============================================================
# 5) Cytoscape export — Full network (Supplementary)
# ============================================================
# Harmonize support column name (optional)
if ("Support Type" %in% colnames(sub_int)) {
  sub_int <- sub_int %>% rename(Support_Type = `Support Type`)
} else if (!("Support_Type" %in% colnames(sub_int))) {
  sub_int <- sub_int %>% mutate(Support_Type = NA_character_)
}

triplets_full <- sub_int %>%
  mutate(
    circRNA = as.character(circRNA),
    miRNA   = as.character(miRNA),
    mRNA    = as.character(mRNA)
  )

nodes_circ_full <- triplets_full %>%
  distinct(id = circRNA) %>%
  mutate(type = "circRNA", color = "#2ECC71", shape = "V")

nodes_miR_full <- triplets_full %>%
  distinct(id = miRNA) %>%
  mutate(type = "miRNA", color = "#3498DB", shape = "ellipse")

nodes_mRNA_full <- triplets_full %>%
  distinct(id = mRNA) %>%
  mutate(type = "mRNA", color = "#FFC000", shape = "diamond")

nodes_full <- bind_rows(nodes_circ_full, nodes_miR_full, nodes_mRNA_full) %>%
  distinct(id, .keep_all = TRUE)

edges_circ_miR_full <- triplets_full %>%
  transmute(
    source      = circRNA,
    target      = miRNA,
    interaction = "circRNA-miRNA",
    Energy      = if ("Energy" %in% colnames(triplets_full)) as.character(Energy) else NA_character_,
    Support     = Support_Type
  ) %>%
  distinct(source, target, .keep_all = TRUE)

edges_miR_mRNA_full <- triplets_full %>%
  transmute(
    source      = miRNA,
    target      = mRNA,
    interaction = "miRNA-mRNA",
    Support     = Support_Type
  ) %>%
  distinct(source, target, .keep_all = TRUE)

edges_full <- bind_rows(edges_circ_miR_full, edges_miR_mRNA_full)

write_tsv(nodes_full, out_nodes_full)
write_tsv(edges_full, out_edges_full)

cat("\nCytoscape (Full network) files saved:\n")
cat("  - ", out_nodes_full, "\n", sep = "")
cat("  - ", out_edges_full, "\n", sep = "")
cat("Full network | nodes:", nrow(nodes_full), " edges:", nrow(edges_full), "\n")

# ============================================================
# 6) Cytoscape export — Curated network (Main figure)
# Rules:
#   - core5 circRNAs: 5 circRNAs with most negative mean Energy
#   - top15 miRNAs: most frequent within core5
#   - top30 mRNAs : most frequent within (core5 + top15 miRNAs)
# ============================================================
if (nrow(energy_by_circ) == 0) {
  warning("Energy-based circRNA ranking not available. Skipping curated network export.")
} else {
  circ_core5 <- energy_by_circ %>% slice_head(n = 5) %>% pull(circRNA)
  
  triplets_core <- triplets_full %>%
    filter(circRNA %in% circ_core5)
  
  miR_top15_tbl <- triplets_core %>%
    count(miRNA, sort = TRUE) %>%
    slice_head(n = 15)
  
  miR_top15 <- miR_top15_tbl$miRNA
  
  triplets_core_miR <- triplets_core %>%
    filter(miRNA %in% miR_top15)
  
  mRNA_top30_tbl <- triplets_core_miR %>%
    count(mRNA, sort = TRUE) %>%
    slice_head(n = 30)
  
  mRNA_top30 <- mRNA_top30_tbl$mRNA
  
  triplets_main <- triplets_core_miR %>%
    filter(mRNA %in% mRNA_top30)
  
  cat("\nCurated network summary:\n")
  cat("  - circRNAs:", n_distinct(triplets_main$circRNA), "\n")
  cat("  - miRNAs :", n_distinct(triplets_main$miRNA), "\n")
  cat("  - mRNAs  :", n_distinct(triplets_main$mRNA), "\n")
  cat("  - triplets:", nrow(triplets_main), "\n")
  
  nodes_circ_main <- triplets_main %>%
    distinct(id = circRNA) %>%
    mutate(type = "circRNA", color = "#2ECC71", shape = "V", is_core5 = TRUE)
  
  nodes_miR_main <- triplets_main %>%
    distinct(id = miRNA) %>%
    mutate(type = "miRNA", color = "#3498DB", shape = "ellipse", is_core5 = TRUE)
  
  nodes_mRNA_main <- triplets_main %>%
    distinct(id = mRNA) %>%
    mutate(type = "mRNA", color = "#FFC000", shape = "diamond", is_core5 = FALSE)
  
  nodes_main <- bind_rows(nodes_circ_main, nodes_miR_main, nodes_mRNA_main) %>%
    distinct(id, .keep_all = TRUE)
  
  edges_circ_miR_main <- triplets_main %>%
    transmute(
      source      = circRNA,
      target      = miRNA,
      interaction = "circRNA-miRNA",
      Energy      = if ("Energy" %in% colnames(triplets_main)) as.character(Energy) else NA_character_
    ) %>%
    distinct(source, target, .keep_all = TRUE)
  
  edges_miR_mRNA_main <- triplets_main %>%
    transmute(
      source      = miRNA,
      target      = mRNA,
      interaction = "miRNA-mRNA"
    ) %>%
    distinct(source, target, .keep_all = TRUE)
  
  edges_main <- bind_rows(edges_circ_miR_main, edges_miR_mRNA_main)
  
  write_tsv(nodes_main, out_nodes_main)
  write_tsv(edges_main, out_edges_main)
  
  cat("\nCytoscape (Curated network) files saved:\n")
  cat("  - ", out_nodes_main, "\n", sep = "")
  cat("  - ", out_edges_main, "\n", sep = "")
  cat("Curated network | nodes:", nrow(nodes_main), " edges:", nrow(edges_main), "\n")
}

# ============================================================
# Stage 12 — GEO / DAVID enrichment plotting
# Purpose:
#   - Export unique mRNAs for GEO/DAVID
#   - Read DAVID outputs (xlsx) and generate barplots for BP/MF/CC/KEGG
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

stage12_dir <- file.path(base_dir, "Stage12_GEO_DAVID")
dir.create(stage12_dir, recursive = TRUE, showWarnings = FALSE)

# Input: filtered triplets (15 circRNAs) exported above
net_file <- out_triplets_15
if (!file.exists(net_file)) stop("Missing network file: ", net_file)

net_df <- read.delim(net_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!("mRNA" %in% colnames(net_df))) stop("Column 'mRNA' not found in: ", net_file)

mrna_unique <- sort(unique(as.character(net_df$mRNA)))
cat("\nStage12 | Unique mRNAs:", length(mrna_unique), "\n")

out_mrna_list <- file.path(stage12_dir, "Unique_mRNAs_for_GEO.txt")
write.table(mrna_unique, file = out_mrna_list, quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("Saved:", out_mrna_list, "\n")

# ============================================================
# DAVID helpers
# ============================================================
prepare_david <- function(path) {
  if (!file.exists(path)) stop("DAVID file not found: ", path)
  
  df <- readxl::read_xlsx(path)
  
  if (!all(c("Term", "Benjamini") %in% colnames(df))) {
    stop("DAVID file missing required columns 'Term' and/or 'Benjamini': ", path)
  }
  
  # Try to detect fold enrichment column name
  fe_col <- c("Fold Enrichment", "Fold.Enrichment", "Fold_Enrichment")
  fe_col <- fe_col[fe_col %in% colnames(df)][1]
  if (is.na(fe_col)) stop("Fold enrichment column not found in: ", path)
  
  df %>%
    mutate(
      Term_clean = as.character(Term),
      Benjamini  = safe_numeric(Benjamini),
      Fold.Enrichment = safe_numeric(.data[[fe_col]]),
      negLog10FDR = -log10(Benjamini)
    ) %>%
    filter(is.finite(Benjamini), is.finite(Fold.Enrichment), is.finite(negLog10FDR)) %>%
    arrange(Benjamini)
}

make_david_barplot <- function(df,
                               category_title,
                               out_pdf,
                               top_n = 10,
                               x_max_fixed = 100,
                               wrap_width = 40,
                               y_lineheight = 1.0) {
  
  df_top <- df %>%
    slice_head(n = top_n) %>%
    mutate(
      Term_plot   = stringr::str_wrap(Term_clean, width = wrap_width),
      Term_factor = factor(Term_plot, levels = rev(Term_plot))
    )
  
  p <- ggplot(df_top, aes(x = negLog10FDR, y = Term_factor, fill = Fold.Enrichment)) +
    geom_col() +
    scale_fill_gradient(name = "Fold\nEnrichment", low = "blue", high = "red", na.value = "grey") +
    scale_x_continuous(
      limits = c(0, x_max_fixed),
      breaks = seq(0, x_max_fixed, by = 20),
      expand = c(0, 0)
    ) +
    labs(title = category_title, x = "-log10(FDR)", y = "Term") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", size = 16, hjust = 0),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.text.y  = element_text(size = 9, lineheight = y_lineheight),
      axis.text.x  = element_text(size = 10),
      legend.position   = "right",
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.border      = element_blank(),
      axis.line.x       = element_line(colour = "black"),
      axis.ticks.x      = element_line(colour = "black"),
      axis.line.y       = element_blank(),
      axis.ticks.y      = element_blank()
    )
  
  ggsave(
    filename = out_pdf,
    plot     = p,
    device   = cairo_pdf,
    width    = 7,
    height   = 5,
    units    = "in"
  )
  
  p
}

# ============================================================
# DAVID inputs (place these in stage12_dir or change paths)
# ============================================================
bp_xlsx   <- file.path(stage12_dir, "BP.xlsx")
mf_xlsx   <- file.path(stage12_dir, "MF.xlsx")
cc_xlsx   <- file.path(stage12_dir, "CC.xlsx")
kegg_xlsx <- file.path(stage12_dir, "KEGG.xlsx")

bp_df   <- prepare_david(bp_xlsx)
mf_df   <- prepare_david(mf_xlsx)
cc_df   <- prepare_david(cc_xlsx)
kegg_df <- prepare_david(kegg_xlsx)

global_x_max_raw <- max(
  bp_df$negLog10FDR,
  mf_df$negLog10FDR,
  cc_df$negLog10FDR,
  kegg_df$negLog10FDR,
  na.rm = TRUE
)
global_x_max <- ceiling(global_x_max_raw / 10) * 10

make_david_barplot(
  bp_df,
  category_title = "a. Biological Processes",
  out_pdf = file.path(stage12_dir, "Fig_DAVID_BP.pdf"),
  top_n = 10,
  x_max_fixed = global_x_max,
  wrap_width = 40,
  y_lineheight = 1.0
)

make_david_barplot(
  mf_df,
  category_title = "b. Molecular Functions",
  out_pdf = file.path(stage12_dir, "Fig_DAVID_MF.pdf"),
  top_n = 10,
  x_max_fixed = global_x_max,
  wrap_width = 40,
  y_lineheight = 1.3
)

make_david_barplot(
  cc_df,
  category_title = "c. Cellular Components",
  out_pdf = file.path(stage12_dir, "Fig_DAVID_CC.pdf"),
  top_n = 10,
  x_max_fixed = global_x_max,
  wrap_width = 40,
  y_lineheight = 1.0
)

make_david_barplot(
  kegg_df,
  category_title = "d. KEGG Pathways",
  out_pdf = file.path(stage12_dir, "Fig_DAVID_KEGG.pdf"),
  top_n = 10,
  x_max_fixed = global_x_max,
  wrap_width = 40,
  y_lineheight = 1.0
)

cat("\nStage11–Stage12 completed.\n")


# ==============================
# END OF SECTION 3.1 (Stages 1–12)
# ==============================

###############################################################################
# Section 3.2 — Differential circRNA expression and DE-refined regulatory network
#
# GitHub-ready pipeline (Stages 1–7)
#   - Loads raw CIRIquant circRNA counts for three breast cancer subtypes
#   - Runs DESeq2 differential expression for two contrasts
#   - Generates volcano plots (Top10 labels + no-label version)
#   - Builds a DE-refined circRNA–miRNA–mRNA network from a reference triplet file
#   - Summarizes network connectivity and IBC-enriched loci
#   - Exports Cytoscape node/edge tables (DE network + curated figure network)
#   - Extracts unique mRNA targets for enrichment and visualizes DAVID/KEGG outputs
#   - Builds circRNA × KEGG pathway target-count heatmap (top 10 pathways)
#
# Reproducibility notes
#   - This script assumes all input paths exist and are consistent with your layout.
#   - Output files are written into stage-specific folders under `section_dir`.
#   - circRNA IDs are harmonized: CIRIquant uses "chr:start|end" while triplets use
#     "chr:start-end". We convert "|" -> "-" before matching.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(readxl)
  library(tools)
  library(scales)
})

# =============================================================================
# Stage 1 — Configuration and raw count loading
# =============================================================================

base_dir    <- "F:/Me/Expresion/Articles/My_Data/Final_Words/R_Analysis"
section_dir <- file.path(base_dir, "Redrawing_shapes", "3.3")

dirs <- list(
  stage1 = file.path(section_dir, "Stage1_LoadCounts"),
  stage2 = file.path(section_dir, "Stage2_Metadata"),
  stage3 = file.path(section_dir, "Stage3_DESeq2"),
  stage4 = file.path(section_dir, "Stage4_Volcano"),
  stage5 = file.path(section_dir, "Stage5_DE_Refined_Network"),
  stage6 = file.path(section_dir, "Stage6_mRNA_for_Enrichment"),
  stage7 = file.path(section_dir, "Stage7_Enrichment_and_Heatmap")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# ---- Input: CIRIquant expressionsfinal.txt (comma-separated) ----------------
# Format:
#   - Row names: circ_id
#   - Columns: samples (numeric counts; 0, 1.0, 2.0, ...)
read_ciri_counts <- function(path) {
  message("Reading counts: ", path)
  df <- read.csv(path, header = TRUE, row.names = 1, check.names = FALSE)
  df[] <- lapply(df, as.numeric)
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  message("Loaded (rows × cols): ", paste(dim(mat), collapse = " × "))
  mat
}

# ---- Define raw count inputs ------------------------------------------------
ciri_inputs <- list(
  IBC  = file.path(base_dir, "GSE207248_Completed_hg38", "CIRIquant_GSE207248", "expressionsfinal.txt"),
  TNBC = file.path(base_dir, "GSE142731_Completed_hg38", "CIRIquant_Final", "expressionsfinal.txt"),
  HER2 = file.path(base_dir, "GSE191230_Completed_hg38", "prerCIRIquant_GSE191230", "expressionsfinal.txt")
)

setwd(dirs$stage1)

IBC_all  <- read_ciri_counts(ciri_inputs$IBC)
TNBC_all <- read_ciri_counts(ciri_inputs$TNBC)
HER2_all <- read_ciri_counts(ciri_inputs$HER2)

# ---- Select tumor samples per subtype --------------------------------------
# IBC: all samples retained (Pre_ and Post_ considered tumor samples for subtype comparison)
IBC <- IBC_all

# TNBC: all samples are tumor samples
TNBC <- TNBC_all

# HER2+: keep columns containing "Tumor" (e.g., Metastatic_Tumor_* and Primary_Tumor_*)
HER2 <- HER2_all[, grepl("Tumor", colnames(HER2_all)), drop = FALSE]

message("Final matrices (rows × cols):")
message("  TNBC: ", paste(dim(TNBC), collapse = " × "))
message("  IBC : ", paste(dim(IBC),  collapse = " × "))
message("  HER2: ", paste(dim(HER2), collapse = " × "))

# Save a lightweight snapshot of sample names
write_lines(colnames(TNBC), file.path(dirs$stage1, "samples_TNBC.txt"))
write_lines(colnames(IBC),  file.path(dirs$stage1, "samples_IBC.txt"))
write_lines(colnames(HER2), file.path(dirs$stage1, "samples_HER2.txt"))

# =============================================================================
# Stage 2 — Metadata construction (subtype + batch)
# =============================================================================

setwd(dirs$stage2)

make_meta <- function(mat, subtype_label, batch_label) {
  data.frame(
    sample  = colnames(mat),
    subtype = subtype_label,
    batch   = batch_label,
    row.names = colnames(mat),
    stringsAsFactors = FALSE
  )
}

meta_TNBC <- make_meta(TNBC, subtype_label = "TNBC", batch_label = "GSE142731")
meta_IBC  <- make_meta(IBC,  subtype_label = "IBC",  batch_label = "GSE207248")
meta_HER2 <- make_meta(HER2, subtype_label = "HER2", batch_label = "GSE191230")

write.csv(meta_TNBC, "metadata_TNBC.csv", row.names = FALSE)
write.csv(meta_IBC,  "metadata_IBC.csv",  row.names = FALSE)
write.csv(meta_HER2, "metadata_HER2.csv", row.names = FALSE)

# =============================================================================
# Stage 3 — Differential expression (DESeq2) with robust filtering
# =============================================================================

setwd(dirs$stage3)

# Filtering and significance cutoffs
padj_cutoff <- 0.05
lfc_cutoff  <- 1

# Harmonize circRNA IDs consistently across downstream steps
normalize_circ_id <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\\|", "-") %>%
    str_replace_all("\u2013|\u2014", "-") %>%
    str_remove_all("\\s+") %>%
    str_trim()
}

run_deseq_pair <- function(mat_group1, mat_group2,
                           meta_group1, meta_group2,
                           group1_name, group2_name,
                           contrast_prefix,
                           padj_cutoff = 0.05,
                           lfc_cutoff  = 1) {
  
  message("\n====================================================")
  message("DESeq2 contrast: ", contrast_prefix)
  message("====================================================")
  
  # 3.1 Common circRNAs
  common_circs <- intersect(rownames(mat_group1), rownames(mat_group2))
  if (length(common_circs) == 0) stop("No shared circRNAs found between groups.")
  
  mat1 <- mat_group1[common_circs, , drop = FALSE]
  mat2 <- mat_group2[common_circs, , drop = FALSE]
  expr_mat <- cbind(mat1, mat2)
  
  meta_all <- rbind(meta_group1, meta_group2)
  meta_all <- meta_all[colnames(expr_mat), , drop = FALSE]
  
  # 3.2 Check batch×subtype confounding
  tab_sub_batch <- table(meta_all$subtype, meta_all$batch)
  message("Subtype × batch table:")
  print(tab_sub_batch)
  
  row_nonzero_counts <- apply(tab_sub_batch > 0, 1, sum)
  if (all(row_nonzero_counts == 1)) {
    message("Batch is fully confounded with subtype; using design = ~ subtype")
    design_formula <- ~ subtype
  } else {
    message("Batch can be modeled; using design = ~ batch + subtype")
    design_formula <- ~ batch + subtype
  }
  
  # 3.3 Strong low-expression filtering
  #   - total counts across all samples >= 10
  #   - count >= 2 in at least 3 samples
  total_reads <- rowSums(expr_mat, na.rm = TRUE)
  samples_with_2_or_more <- rowSums(expr_mat >= 2, na.rm = TRUE)
  keep_rows <- total_reads >= 10 & samples_with_2_or_more >= 3
  
  message("circRNAs before filter: ", nrow(expr_mat))
  message("circRNAs retained:      ", sum(keep_rows))
  
  expr_filtered <- expr_mat[keep_rows, , drop = FALSE]
  if (nrow(expr_filtered) == 0) stop("No circRNAs remained after filtering.")
  
  # 3.4 Drop zero-sum samples
  col_sums <- colSums(expr_filtered)
  if (any(col_sums == 0)) {
    zero_samples <- names(col_sums[col_sums == 0])
    message("Dropping zero-sum samples: ", paste(zero_samples, collapse = ", "))
    expr_filtered <- expr_filtered[, col_sums > 0, drop = FALSE]
    meta_all <- meta_all[colnames(expr_filtered), , drop = FALSE]
  }
  
  # 3.5 Integer counts for DESeq2
  expr_filtered <- round(expr_filtered)
  storage.mode(expr_filtered) <- "integer"
  
  meta_all$subtype <- factor(meta_all$subtype)
  meta_all$batch   <- factor(meta_all$batch)
  
  dds <- DESeqDataSetFromMatrix(
    countData = expr_filtered,
    colData   = meta_all,
    design    = design_formula
  )
  
  dds <- DESeq(dds, sfType = "poscounts")
  
  # Ensure reference level is group2 (so log2FC is group1 vs group2)
  dds$subtype <- relevel(dds$subtype, ref = group2_name)
  res <- results(dds, contrast = c("subtype", group1_name, group2_name))
  res <- res[order(res$pvalue), ]
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("circ_id") %>%
    mutate(
      circ_id = as.character(circ_id),
      status = case_when(
        !is.na(padj) & padj < padj_cutoff & log2FoldChange >=  lfc_cutoff ~ paste0("Up in ", group1_name),
        !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ paste0("Up in ", group2_name),
        TRUE ~ "NS"
      )
    )
  
  out_file <- paste0("DESeq2_results_", contrast_prefix, ".csv")
  write.csv(res_df, out_file, row.names = FALSE)
  message("Saved: ", file.path(getwd(), out_file))
  
  res_df
}

# ---- Run the two primary contrasts -----------------------------------------
res_TNBC_vs_IBC <- run_deseq_pair(
  mat_group1      = TNBC,
  mat_group2      = IBC,
  meta_group1     = meta_TNBC,
  meta_group2     = meta_IBC,
  group1_name     = "TNBC",
  group2_name     = "IBC",
  contrast_prefix = "TNBC_vs_IBC",
  padj_cutoff     = padj_cutoff,
  lfc_cutoff      = lfc_cutoff
)

res_IBC_vs_HER2 <- run_deseq_pair(
  mat_group1      = IBC,
  mat_group2      = HER2,
  meta_group1     = meta_IBC,
  meta_group2     = meta_HER2,
  group1_name     = "IBC",
  group2_name     = "HER2",
  contrast_prefix = "IBC_vs_HER2",
  padj_cutoff     = padj_cutoff,
  lfc_cutoff      = lfc_cutoff
)

message("DE circRNAs (TNBC_vs_IBC): ", sum(res_TNBC_vs_IBC$status != "NS"))
message("DE circRNAs (IBC_vs_HER2): ", sum(res_IBC_vs_HER2$status != "NS"))

# =============================================================================
# Stage 4 — Volcano plots with Top10 labeling
# =============================================================================

setwd(dirs$stage4)

load_deseq <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(all(c("circ_id", "pvalue", "padj", "log2FoldChange") %in% colnames(df)))
  df
}

res_TNBC_vs_IBC <- load_deseq(file.path(dirs$stage3, "DESeq2_results_TNBC_vs_IBC.csv"))
res_IBC_vs_HER2 <- load_deseq(file.path(dirs$stage3, "DESeq2_results_IBC_vs_HER2.csv"))

make_volcano_top10 <- function(res_df,
                               contrast_name,
                               group1_label,
                               group2_label,
                               padj_cutoff = 0.05,
                               lfc_cutoff  = 1,
                               out_prefix  = NULL) {
  
  if (is.null(out_prefix)) out_prefix <- paste0("Volcano_", contrast_name)
  
  df <- res_df %>%
    filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
    mutate(
      status = case_when(
        !is.na(padj) & padj < padj_cutoff & log2FoldChange >=  lfc_cutoff ~ paste0("Up in ", group1_label),
        !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ paste0("Up in ", group2_label),
        TRUE ~ "NS"
      ),
      minusLog10P = -log10(pvalue + 1e-300)
    )
  
  sig_df <- df %>% filter(status != "NS")
  top10 <- sig_df %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_head(n = 10)
  
  df <- df %>%
    mutate(is_top10 = circ_id %in% top10$circ_id)
  
  col_map <- c(
    "NS" = "black",
    paste0("Up in ", group1_label) = "#0072B2",
    paste0("Up in ", group2_label) = "#FF7F00"
  )
  
  base_theme <- theme_minimal(base_size = 13) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x      = element_line(colour = "black", linewidth = 0.5),
      axis.line.y      = element_line(colour = "black", linewidth = 0.5),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      legend.position  = "top",
      plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title       = element_text(size = 13, face = "bold")
    )
  
  p_core <- ggplot(df, aes(x = log2FoldChange, y = minusLog10P)) +
    geom_point(data = df %>% filter(status == "NS"),
               aes(color = status), size = 2.5, alpha = 0.7) +
    geom_point(data = df %>% filter(status != "NS"),
               aes(color = status), size = 4, alpha = 0.85) +
    geom_point(data = df %>% filter(is_top10),
               shape = 21, size = 5, stroke = 1.0, color = "black", fill = NA) +
    scale_color_manual(values = col_map, name = "Status") +
    labs(
      title = paste("Volcano plot —", contrast_name),
      x = "log2(Fold Change)",
      y = "-log10(p-value)"
    ) +
    base_theme
  
  p_labels <- p_core +
    geom_text_repel(
      data = df %>% filter(is_top10),
      aes(label = circ_id),
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3
    )
  
  # Save both versions
  ggsave(
    filename    = paste0(out_prefix, "_Top10_labels.tiff"),
    plot        = p_labels,
    dpi         = 1200,
    width       = 8,
    height      = 6,
    units       = "in",
    compression = "lzw"
  )
  ggsave(
    filename    = paste0(out_prefix, "_Top10_no_labels.tiff"),
    plot        = p_core,
    dpi         = 1200,
    width       = 8,
    height      = 6,
    units       = "in",
    compression = "lzw"
  )
  
  list(df = df, sig = sig_df, top10 = top10)
}

TNBC_IBC_vol <- make_volcano_top10(
  res_df        = res_TNBC_vs_IBC,
  contrast_name = "TNBC_vs_IBC",
  group1_label  = "TNBC",
  group2_label  = "IBC",
  padj_cutoff   = padj_cutoff,
  lfc_cutoff    = lfc_cutoff
)

IBC_HER2_vol <- make_volcano_top10(
  res_df        = res_IBC_vs_HER2,
  contrast_name = "IBC_vs_HER2",
  group1_label  = "IBC",
  group2_label  = "HER2",
  padj_cutoff   = padj_cutoff,
  lfc_cutoff    = lfc_cutoff
)

# =============================================================================
# Stage 5 — DE-refined circRNA–miRNA–mRNA network + Cytoscape exports
# =============================================================================

setwd(dirs$stage5)

triplets_path <- file.path(
  base_dir,
  "Mutual_Analysis", "circRNA–miRNA", "circRNA–miRNA–mRNA", "all",
  "All_circRNA_miRNA_mRNA.tsv"
)

if (!file.exists(triplets_path)) {
  stop("Missing triplets reference file: ", triplets_path)
}

triplets_all <- read_tsv(triplets_path, show_col_types = FALSE)

# Standardize mRNA column name
if ("mRNA" %in% colnames(triplets_all)) {
  triplets_all <- triplets_all %>% rename(mRNA_target = mRNA)
}
stopifnot(all(c("circRNA", "miRNA", "mRNA_target") %in% colnames(triplets_all)))

# 5.1 Identify DE circRNAs across both contrasts (harmonized IDs)
sig_TNBC_IBC <- res_TNBC_vs_IBC %>%
  filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)

sig_IBC_HER2 <- res_IBC_vs_HER2 %>%
  filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)

all_DE_circs <- unique(c(sig_TNBC_IBC$circ_id, sig_IBC_HER2$circ_id)) %>%
  normalize_circ_id()

message("Total DE circRNAs (unique; harmonized): ", length(all_DE_circs))

if (length(all_DE_circs) == 0) stop("No DE circRNAs found with current cutoffs.")

# 5.2 Restrict triplets to DE circRNAs
triplets_DE <- triplets_all %>%
  mutate(circRNA = normalize_circ_id(circRNA)) %>%
  filter(circRNA %in% all_DE_circs) %>%
  distinct(circRNA, miRNA, mRNA_target, .keep_all = TRUE)

if (nrow(triplets_DE) == 0) {
  stop("DE-based network is empty. Check circRNA ID formats and overlap.")
}

write_tsv(triplets_DE, "Network_triplets_DE_circRNAs.tsv")

message("DE network summary:")
message("  circRNAs: ", n_distinct(triplets_DE$circRNA))
message("  miRNAs  : ", n_distinct(triplets_DE$miRNA))
message("  mRNAs   : ", n_distinct(triplets_DE$mRNA_target))
message("  triplets: ", nrow(triplets_DE))

# 5.3 Connectivity summary per circRNA (supplementary)
circ_connectivity <- triplets_DE %>%
  group_by(circRNA) %>%
  summarise(
    n_miRNA   = n_distinct(miRNA),
    n_mRNA    = n_distinct(mRNA_target),
    n_triplet = n(),
    .groups   = "drop"
  ) %>%
  arrange(desc(n_mRNA), desc(n_miRNA))

write.csv(circ_connectivity,
          "Supp_Table_DE_circRNA_connectivity_network.csv",
          row.names = FALSE)

# 5.4 IBC-enriched flag from DE directions across contrasts
# TNBC_vs_IBC: "Up in IBC" => log2FC < -lfc_cutoff
# IBC_vs_HER2: "Up in IBC" => log2FC >  lfc_cutoff  (because group1=IBC, group2=HER2)
IBC_status <- bind_rows(
  res_TNBC_vs_IBC %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    transmute(circRNA = normalize_circ_id(circ_id), IBC_enriched = TRUE),
  
  res_IBC_vs_HER2 %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange >  lfc_cutoff) %>%
    transmute(circRNA = normalize_circ_id(circ_id), IBC_enriched = TRUE)
) %>% distinct()

circ_connectivity_IBC <- circ_connectivity %>%
  left_join(IBC_status, by = "circRNA") %>%
  mutate(IBC_enriched = ifelse(is.na(IBC_enriched), FALSE, IBC_enriched))

write.csv(
  circ_connectivity_IBC,
  "Supp_Table_DE_circRNA_connectivity_with_IBC_flag.csv",
  row.names = FALSE
)

# 5.5 Table: circRNAs enriched in IBC in BOTH contrasts (if any)
TNBC_IBC_IBCenr <- res_TNBC_vs_IBC %>%
  filter(!is.na(padj), padj < padj_cutoff, log2FoldChange <= -lfc_cutoff) %>%
  transmute(circRNA = normalize_circ_id(circ_id), log2FC_TNBC_IBC = log2FoldChange)

IBC_HER2_IBCenr <- res_IBC_vs_HER2 %>%
  filter(!is.na(padj), padj < padj_cutoff, log2FoldChange >=  lfc_cutoff) %>%
  transmute(circRNA = normalize_circ_id(circ_id), log2FC_IBC_HER2 = log2FoldChange)

IBC_both <- inner_join(TNBC_IBC_IBCenr, IBC_HER2_IBCenr, by = "circRNA")

network_summary <- triplets_DE %>%
  group_by(circRNA) %>%
  summarise(
    miRNA = n_distinct(miRNA),
    mRNA  = n_distinct(mRNA_target),
    .groups = "drop"
  )

table_IBC_both <- IBC_both %>%
  left_join(network_summary, by = "circRNA") %>%
  mutate(
    Comparisons = "TNBC vs IBC; IBC vs HER2",
    log2FC = paste0(
      "TNBCvsIBC: ", round(log2FC_TNBC_IBC, 2),
      "; IBCvsHER2: ", round(log2FC_IBC_HER2, 2)
    )
  ) %>%
  select(circRNA, Comparisons, log2FC, miRNA, mRNA) %>%
  arrange(desc(mRNA), desc(miRNA))

write.csv(table_IBC_both, "Table_IBC_enriched_in_both_contrasts.csv", row.names = FALSE)

# 5.6 Cytoscape exports for a curated "main-figure" network
# Strategy:
#   - pick top 5 circRNAs by mRNA targets (optionally restricted to IBC-enriched)
#   - within those, pick top 10 shared/recurring miRNAs
#   - then pick top 20 multi-regulated mRNAs (n_miRNA >= 2)
circ_core_main <- circ_connectivity_IBC %>%
  filter(IBC_enriched) %>%
  arrange(desc(n_mRNA), desc(n_miRNA)) %>%
  slice_head(n = 5) %>%
  pull(circRNA)

if (length(circ_core_main) == 0) {
  message("No IBC-enriched circRNAs in network; falling back to top 5 by mRNA degree.")
  circ_core_main <- circ_connectivity %>%
    arrange(desc(n_mRNA), desc(n_miRNA)) %>%
    slice_head(n = 5) %>%
    pull(circRNA)
}

miRNA_core <- triplets_DE %>%
  filter(circRNA %in% circ_core_main) %>%
  count(miRNA, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(miRNA)

mRNA_core <- triplets_DE %>%
  filter(circRNA %in% circ_core_main, miRNA %in% miRNA_core) %>%
  group_by(mRNA_target) %>%
  summarise(
    n_miRNA  = n_distinct(miRNA),
    n_edges  = n(),
    .groups  = "drop"
  ) %>%
  filter(n_miRNA >= 2) %>%
  arrange(desc(n_miRNA), desc(n_edges)) %>%
  slice_head(n = 20) %>%
  pull(mRNA_target)

triplets_fig <- triplets_DE %>%
  filter(circRNA %in% circ_core_main,
         miRNA %in% miRNA_core,
         mRNA_target %in% mRNA_core)

write_tsv(triplets_fig, "Network_triplets_Figure_main.tsv")

# Nodes
circ_degree_fig <- triplets_fig %>% count(circRNA, name = "degree")
miRNA_degree_fig <- triplets_fig %>% count(miRNA, name = "degree")
mRNA_degree_fig <- triplets_fig %>% count(mRNA_target, name = "degree")

nodes_circ <- circ_degree_fig %>%
  transmute(
    id = circRNA,
    type = "circRNA",
    shape = "V",
    color = "#2ECC71",
    size = scales::rescale(degree, to = c(60, 90))
  )

nodes_miR <- miRNA_degree_fig %>%
  transmute(
    id = miRNA,
    type = "miRNA",
    shape = "ellipse",
    color = "#3498DB",
    size = scales::rescale(degree, to = c(30, 55))
  )

nodes_mRNA <- mRNA_degree_fig %>%
  transmute(
    id = mRNA_target,
    type = "mRNA",
    shape = "diamond",
    color = "#F1C40F",
    size = scales::rescale(degree, to = c(20, 40))
  )

nodes_fig <- bind_rows(nodes_circ, nodes_miR, nodes_mRNA) %>%
  distinct(id, .keep_all = TRUE) %>%
  arrange(type, desc(size))

write_tsv(nodes_fig, "Cytoscape_nodes_Figure_main.tsv")

# Edges
edges_circ_miR <- triplets_fig %>%
  transmute(source = circRNA, target = miRNA, interaction = "circRNA-miRNA") %>%
  distinct()

edges_miR_mRNA <- triplets_fig %>%
  transmute(source = miRNA, target = mRNA_target, interaction = "miRNA-mRNA") %>%
  distinct()

edges_fig <- bind_rows(edges_circ_miR, edges_miR_mRNA) %>%
  distinct(source, target, interaction)

write_tsv(edges_fig, "Cytoscape_edges_Figure_main.tsv")

# =============================================================================
# Stage 6 — Extract unique mRNAs from the DE-refined network
# =============================================================================

setwd(dirs$stage6)

triplets_DE <- read_tsv(file.path(dirs$stage5, "Network_triplets_DE_circRNAs.tsv"),
                        show_col_types = FALSE)

mRNA_unique <- triplets_DE %>%
  distinct(mRNA_target) %>%
  arrange(mRNA_target)

write_csv(mRNA_unique, "unique_mRNAs_from_DE_network.csv")
write_lines(mRNA_unique$mRNA_target, "unique_mRNAs_from_DE_network.txt")

message("Unique mRNAs exported: ", nrow(mRNA_unique))

# =============================================================================
# Stage 7 — Enrichment plotting (DAVID outputs) + circRNA×KEGG heatmap
# =============================================================================

setwd(dirs$stage7)

# ---- Stage 7.1: Read DAVID chart files (xlsx or txt) -----------------------
read_david_chart <- function(path) {
  ext <- tools::file_ext(path)
  
  if (ext %in% c("xls", "xlsx")) {
    df <- readxl::read_xlsx(path)
  } else {
    df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # Choose the best available multiple-testing column
  if ("FDR" %in% names(df)) {
    df$FDR_use <- as.numeric(df$FDR)
  } else if ("Benjamini" %in% names(df)) {
    df$FDR_use <- as.numeric(df$Benjamini)
  } else if ("PValue" %in% names(df)) {
    df$FDR_use <- as.numeric(df$PValue)
  } else if ("P-Value" %in% names(df)) {
    df$FDR_use <- as.numeric(df$`P-Value`)
  } else {
    stop("No FDR/Benjamini/PValue column found in: ", path)
  }
  
  df$FDR_safe <- pmax(df$FDR_use, 1e-300)
  
  # Fold enrichment column name may vary
  if (!"Fold.Enrichment" %in% names(df)) {
    fe_col <- grep("Fold", names(df), value = TRUE)[1]
    df$Fold.Enrichment <- if (!is.na(fe_col)) as.numeric(df[[fe_col]]) else NA_real_
  }
  
  if (!"Term" %in% names(df)) stop("No 'Term' column found in: ", path)
  
  df %>%
    filter(!is.na(FDR_safe)) %>%
    mutate(
      Term_clean  = sub("^[^~]*~", "", Term),
      negLog10FDR = -log10(FDR_safe)
    ) %>%
    arrange(FDR_use)
}

make_david_barplot <- function(df,
                               category_title,
                               out_pdf,
                               top_n = 10,
                               x_max_fixed = 100,
                               wrap_width = 40,
                               y_lineheight = 1.0) {
  
  df_top <- df %>%
    slice_head(n = top_n) %>%
    mutate(Term_plot = str_wrap(Term_clean, width = wrap_width)) %>%
    distinct(Term_plot, .keep_all = TRUE) %>%
    mutate(Term_factor = factor(Term_plot, levels = rev(Term_plot)))
  
  p <- ggplot(df_top, aes(x = negLog10FDR, y = Term_factor, fill = Fold.Enrichment)) +
    geom_col() +
    scale_fill_gradient(name = "Fold\nEnrichment",
                        low = "#b3cde3", high = "#081d58", na.value = "grey80") +
    scale_x_continuous(limits = c(0, x_max_fixed),
                       breaks = seq(0, x_max_fixed, by = 20),
                       expand = c(0, 0)) +
    labs(title = category_title, x = "-log10(FDR)", y = "Term") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      axis.text.y = element_text(size = 9, lineheight = y_lineheight),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.ticks.x = element_line(colour = "black"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  ggsave(out_pdf, plot = p, device = cairo_pdf, width = 7, height = 5, units = "in")
  p
}

# ---- Provide DAVID output file names (adjust to your repository) ------------
bp_file   <- "BIOLOGICAL_PROCESS.xlsx"
mf_file   <- "MOLECULAR_FUNCTION.xlsx"
cc_file   <- "CELLULAR_COMPONENT.xlsx"
kegg_file <- "KEGG.xlsx"

bp_geo   <- read_david_chart(bp_file)
mf_geo   <- read_david_chart(mf_file)
cc_geo   <- read_david_chart(cc_file)
kegg_geo <- read_david_chart(kegg_file)

write.csv(bp_geo   %>% slice_head(n = 30), "Supp_BP_top30.csv",   row.names = FALSE)
write.csv(mf_geo   %>% slice_head(n = 30), "Supp_MF_top30.csv",   row.names = FALSE)
write.csv(cc_geo   %>% slice_head(n = 30), "Supp_CC_top30.csv",   row.names = FALSE)
write.csv(kegg_geo %>% slice_head(n = 30), "Supp_KEGG_top30.csv", row.names = FALSE)

global_x_max <- max(
  bp_geo$negLog10FDR, mf_geo$negLog10FDR, cc_geo$negLog10FDR, kegg_geo$negLog10FDR,
  na.rm = TRUE
)
global_x_max <- ceiling(global_x_max / 20) * 20

p_bp   <- make_david_barplot(bp_geo,   "A. Biological Process",   "Fig_BP.pdf",   top_n = 10, x_max_fixed = global_x_max)
p_mf   <- make_david_barplot(mf_geo,   "B. Molecular Function",   "Fig_MF.pdf",   top_n = 10, x_max_fixed = global_x_max, y_lineheight = 1.3)
p_cc   <- make_david_barplot(cc_geo,   "C. Cellular Component",   "Fig_CC.pdf",   top_n = 10, x_max_fixed = global_x_max)
p_kegg <- make_david_barplot(kegg_geo, "D. KEGG Pathways",        "Fig_KEGG.pdf", top_n = 10, x_max_fixed = global_x_max)

# ---- Stage 7.2: circRNA × KEGG heatmap (target counts) ----------------------
triplets_DE <- read_tsv(file.path(dirs$stage5, "Network_triplets_DE_circRNAs.tsv"),
                        show_col_types = FALSE)

circ_mrna <- triplets_DE %>%
  distinct(circRNA, mRNA_target) %>%
  filter(!is.na(mRNA_target), mRNA_target != "")

# Provide the 10 KEGG pathways and corresponding Excel files (one per pathway)
kegg_top10_names <- c(
  "Pathways in cancer",
  "Proteoglycans in cancer",
  "Breast cancer",
  "Gastric cancer",
  "MAPK signaling pathway",
  "PI3K-Akt signaling pathway",
  "Signaling pathways regulating pluripotency of stem cells",
  "Hepatitis B",
  "AGE-RAGE signaling pathway in diabetic complications",
  "Lipid and atherosclerosis"
)

kegg_files <- c(
  "Pathways in cancer.xlsx",
  "Proteoglycans in cancer.xlsx",
  "Breast cancer.xlsx",
  "Gastric cancer.xlsx",
  "MAPK signaling pathway.xlsx",
  "PI3K-Akt signaling pathway.xlsx",
  "Signaling pathways regulating pluripotency of stem cells.xlsx",
  "Hepatitis B.xlsx",
  "AGE-RAGE signaling pathway in diabetic complications.xlsx",
  "Lipid and atherosclerosis.xlsx"
)

read_kegg_gene_list <- function(file, pathway_name) {
  df <- readxl::read_xlsx(file)
  
  gene_col <- dplyr::case_when(
    "Gene Symbol" %in% names(df) ~ "Gene Symbol",
    "Gene"        %in% names(df) ~ "Gene",
    "Symbol"      %in% names(df) ~ "Symbol",
    TRUE                          ~ names(df)[1]
  )
  
  tibble(
    Pathway = pathway_name,
    Gene    = df[[gene_col]] %>% as.character() %>% str_trim()
  ) %>%
    filter(!is.na(Gene), Gene != "")
}

kegg_gene_map <- bind_rows(purrr::map2(kegg_files, kegg_top10_names, read_kegg_gene_list))

circ_path <- circ_mrna %>%
  inner_join(kegg_gene_map, by = c("mRNA_target" = "Gene"))

if (nrow(circ_path) == 0) stop("No overlap between network mRNA targets and KEGG gene lists.")

circ_path_counts <- circ_path %>%
  count(circRNA, Pathway, name = "n_targets")

heat_mat <- circ_path_counts %>%
  pivot_wider(names_from = Pathway, values_from = n_targets, values_fill = 0) %>%
  as.data.frame()

for (p in kegg_top10_names) if (!p %in% colnames(heat_mat)) heat_mat[[p]] <- 0
heat_mat <- heat_mat[, c("circRNA", kegg_top10_names)]

rownames(heat_mat) <- heat_mat$circRNA
heat_mat$circRNA <- NULL

write.csv(circ_path_counts, "circRNA_KEGG_target_counts_long.csv", row.names = FALSE)
write.csv(heat_mat, "circRNA_KEGG_target_counts_matrix.csv", row.names = TRUE)

heat_mat_log <- log10(as.matrix(heat_mat) + 1)
my_cols <- colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100)

pheatmap(
  heat_mat_log,
  color        = my_cols,
  cluster_rows = nrow(heat_mat_log) >= 2,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize     = 8,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col    = 45,
  main         = "DE circRNAs vs Top 10 KEGG pathways (target counts)",
  filename     = "Fig_circRNA_KEGG_heatmap.pdf",
  width        = 7,
  height       = 5
)

# ==============================
# END OF SECTION 3.2
# ==============================
