## ============================================================
## One-click pipeline:
## MetaPhlAn4 taxonomy.tsv + metadata (High / Low calcium)
## Species-level (s__) LEfSe using microeco
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(microeco)
})

## ----------------------------
## Step 0: file paths
## ----------------------------
meta_file <- "metadata.txt"     # columns: sample_id, group
tax_file  <- "taxonomy.tsv"     # MetaPhlAn4 output

out_tax_mat   <- "taxonomy_species_matrix.tsv"
out_lefse_csv <- "lefse_species_high_vs_low.csv"

## ----------------------------
## Step 1: read metadata (High / Low only)
## ----------------------------
GROUP <- read.table(
  meta_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(all(c("sample_id", "group") %in% colnames(GROUP)))
rownames(GROUP) <- GROUP$sample_id

# keep High / Low only (case-insensitive)
grp_lower <- tolower(trimws(GROUP$group))
keep_idx  <- grp_lower %in% c("high", "low")

GROUP <- GROUP[keep_idx, , drop = FALSE]
GROUP$group <- factor(GROUP$group, levels = c("Low", "High"))

cat(sprintf(
  "[Step1] metadata loaded: %d samples (High/Low)\n",
  nrow(GROUP)
))

## ----------------------------
## Step 2: read taxonomy.tsv (MetaPhlAn4-safe)
## ----------------------------
tx_lines <- read_lines(tax_file)

# remove MetaPhlAn header line if present
if (length(tx_lines) >= 2 && grepl("^#metaphlan", tx_lines[2], ignore.case = TRUE)) {
  tx_lines <- tx_lines[-2]
}

TAX <- read_tsv(I(tx_lines), show_col_types = FALSE)

# first column = feature ID
colnames(TAX)[1] <- "Feature"
TAX <- as.data.frame(TAX)
rownames(TAX) <- TAX$Feature
TAX <- TAX[, -1, drop = FALSE]

# drop UNMAPPED / UNINTEGRATED
rm_rows <- grepl("UNMAPPED|UNINTEGRATED", rownames(TAX), ignore.case = TRUE)
if (any(rm_rows)) {
  TAX <- TAX[!rm_rows, , drop = FALSE]
  cat(sprintf("[Step2] removed %d UNMAPPED/UNINTEGRATED rows\n", sum(rm_rows)))
}

# clean column names
colnames(TAX) <- str_trim(colnames(TAX))

# force numeric, NA -> 0
TAX[] <- lapply(TAX, function(x) {
  suppressWarnings(xn <- as.numeric(x))
  xn[is.na(xn)] <- 0
  xn
})

## ----------------------------
## Step 3: keep species level only (s__)
## ----------------------------
is_species <- grepl("\\|s__|^s__", rownames(TAX))

if (!any(is_species)) {
  stop("No s__ features detected. Check MetaPhlAn output.")
}

TAX_s <- TAX[is_species, , drop = FALSE]

cat(sprintf(
  "[Step3] kept species-level features (s__): %d\n",
  nrow(TAX_s)
))

## ----------------------------
## Step 4: align samples strictly
## ----------------------------
common_samples <- intersect(colnames(TAX_s), rownames(GROUP))

if (length(common_samples) < 2) {
  stop("Too few overlapping samples between taxonomy and metadata.")
}

common_samples <- sort(common_samples)

TAX_s   <- TAX_s[, common_samples, drop = FALSE]
GROUP_s <- GROUP[common_samples, , drop = FALSE]

stopifnot(identical(colnames(TAX_s), rownames(GROUP_s)))

# save species abundance matrix (feature × sample)
write.table(
  TAX_s,
  file = out_tax_mat,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat(sprintf(
  "[Step4] species matrix written: %s (%d samples)\n",
  out_tax_mat,
  ncol(TAX_s)
))

## ----------------------------
## Step 5: build microeco object
## ----------------------------
sample_table <- data.frame(
  group = GROUP_s$group,
  row.names = rownames(GROUP_s),
  check.names = FALSE
)

# minimal tax_table: feature name as Taxon
tax_table_min <- data.frame(
  Taxon = rownames(TAX_s),
  row.names = rownames(TAX_s),
  check.names = FALSE
)

micro_dat <- microtable$new(
  sample_table = sample_table,
  otu_table    = TAX_s,
  tax_table    = tax_table_min,
  auto_tidy    = FALSE
)

## ----------------------------
## Step 6: LEfSe (High vs Low)
## ----------------------------
lefse_res <- trans_diff$new(
  dataset         = micro_dat,
  method          = "lefse",
  group           = "group",
  alpha           = 0.05,
  lda_cutoff      = 2.0,
  p_adjust_method = "none"
)

write.csv(
  lefse_res$res_diff,
  out_lefse_csv,
  row.names = FALSE
)

cat(sprintf(
  "[Step6] LEfSe finished. Results saved to: %s\n",
  out_lefse_csv
))

cat("\nPipeline completed successfully (species-level LEfSe, High vs Low calcium)\n")
