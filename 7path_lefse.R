suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(microeco)
})

## ---------- Step 0: paths ----------
meta_file <- "metadata.txt"                      # sample_id, group
path_file <- "pathabundance_relab_stratified.tsv"
out_csv   <- "lefse_pathway_high_vs_low.csv"
out_mat   <- "pathway_unstratified_matrix.tsv"   # 可选：保存聚合后的矩阵方便复现

diff_groups <- c("High", "Low")

## ---------- helper functions ----------
extract_sample_id <- function(x) {
  # 从 "F1945_combined" / "F1945" / 其他混杂形式里抽出 F数字
  str_extract(x, "F\\d+")
}

extract_pathway <- function(x) {
  # stratified: "PWY-2941: xxx|g__...s__..."
  # unstratified: "PWY-2941: xxx"
  str_split_fixed(x, "\\|", 2)[, 1]
}

safe_lefse <- function(dataset, group_col = "group") {
  tryCatch(
    trans_diff$new(
      dataset = dataset,
      method  = "lefse",
      group   = group_col,
      alpha   = 0.05,
      lda_cutoff = 2.0,          # 你要改阈值就改这行
      p_adjust_method = "none"
    ),
    error = function(e) {
      message("LEfSe 出错：", conditionMessage(e))
      return(NULL)
    }
  )
}

## ---------- Step 1: read metadata (High/Low only) ----------
GROUP <- read_tsv(meta_file, show_col_types = FALSE) %>%
  rename(sample_id = 1, group = 2)

GROUP <- GROUP %>%
  mutate(group = str_trim(group)) %>%
  filter(tolower(group) %in% tolower(diff_groups)) %>%
  mutate(group = factor(group, levels = c("Low", "High")))

stopifnot(all(c("sample_id","group") %in% colnames(GROUP)))
GROUP <- as.data.frame(GROUP)
rownames(GROUP) <- GROUP$sample_id

cat(sprintf("[Step1] metadata: %d samples kept (High/Low)\n", nrow(GROUP)))

## ---------- Step 2: read pathway table ----------
# comment.char="" 保留以 # 开头的真实列名行（你的第一列是 "# Pathway"）
PATH_raw <- read.table(
  path_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

# 第一列是 pathway ID（可能叫 "# Pathway" 或 "Pathway"）
first_col <- colnames(PATH_raw)[1]
rownames(PATH_raw) <- PATH_raw[[1]]
PATH_raw <- PATH_raw[, -1, drop = FALSE]

# 列名抽 sample_id（Fxxxx）
sid <- extract_sample_id(colnames(PATH_raw))
keep_c <- !is.na(sid)
PATH_raw <- PATH_raw[, keep_c, drop = FALSE]
colnames(PATH_raw) <- sid[keep_c]

# 若抽完后出现重复列名（极少见，但保险），合并同名列（取和）
if (any(duplicated(colnames(PATH_raw)))) {
  PATH_raw <- as.data.frame(PATH_raw)
  PATH_raw$.__row__ <- rownames(PATH_raw)
  PATH_raw <- PATH_raw %>%
    tidyr::pivot_longer(-.__row__, names_to = "sample_id", values_to = "val") %>%
    group_by(.__row__, sample_id) %>%
    summarise(val = sum(as.numeric(val), na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample_id, values_from = val)
  rn <- PATH_raw$.__row__
  PATH_raw <- PATH_raw %>% select(-.__row__)
  PATH_raw <- as.data.frame(PATH_raw)
  rownames(PATH_raw) <- rn
}

# 去掉 UNMAPPED/UNINTEGRATED 行
rm_rows <- grepl("UNMAPPED|UNINTEGRATED", rownames(PATH_raw), ignore.case = TRUE)
if (any(rm_rows)) {
  PATH_raw <- PATH_raw[!rm_rows, , drop = FALSE]
  cat(sprintf("[Step2] removed %d UNMAPPED/UNINTEGRATED rows\n", sum(rm_rows)))
}

# 全部转数值，NA -> 0
PATH_raw[] <- lapply(PATH_raw, function(x) {
  xn <- suppressWarnings(as.numeric(x))
  xn[is.na(xn)] <- 0
  xn
})

## ---------- Step 3: stratified -> unstratified (collapse by pathway) ----------
path_names <- extract_pathway(rownames(PATH_raw))

PATH_un <- PATH_raw %>%
  as.data.frame() %>%
  mutate(Pathway = path_names) %>%
  group_by(Pathway) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

rownames(PATH_un) <- PATH_un$Pathway
PATH_un <- PATH_un[, setdiff(colnames(PATH_un), "Pathway"), drop = FALSE]

cat(sprintf("[Step3] collapsed to unstratified pathways: %d pathways\n", nrow(PATH_un)))

## ---------- Step 4: align samples with metadata ----------
common <- intersect(colnames(PATH_un), rownames(GROUP))
if (length(common) < 2) stop("交集样本数不足，无法做差异分析。")

common <- sort(common)
PATH_un <- PATH_un[, common, drop = FALSE]
GROUP2  <- GROUP[common, , drop = FALSE]

stopifnot(identical(colnames(PATH_un), rownames(GROUP2)))
cat(sprintf("[Step4] aligned samples: %d\n", length(common)))

# 可选：保存聚合后的矩阵
write.table(PATH_un, out_mat, sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf("[Step4] wrote unstratified pathway matrix: %s\n", out_mat))

## ---------- Step 5: microeco LEfSe ----------
sample_table <- data.frame(
  group = GROUP2$group,
  row.names = rownames(GROUP2),
  check.names = FALSE
)

df_path <- microtable$new(
  sample_table = sample_table,
  otu_table    = PATH_un,
  auto_tidy    = FALSE
)
df_path$taxa_abund <- list(OTU = PATH_un)

lefse_res <- safe_lefse(df_path, group_col = "group")

if (!is.null(lefse_res)) {
  write.csv(lefse_res$res_diff, out_csv, row.names = FALSE)
  message("[Step5] LEfSe done. Output: ", out_csv)
} else {
  message("[Step5] LEfSe failed or no results.")
}
