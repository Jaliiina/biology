suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

# =============================
# INPUTS
# =============================
tax_file  <- "taxonomy.tsv"
pwy_file  <- "pathabundance_relab_stratified.tsv"  # HUMAnN3 stratified table
out_csv   <- "species_pathway_spearman_ponly.csv"

# =============================
# YOUR SIGNIFICANT SETS (去重)
# =============================
sig_species <- unique(c(
  "Lachnospiraceae_bacterium",
  "Holdemanella_porci",
  "Brucella_pituitosa",
  "Bacteroides_nordii",
  "Actinomyces_sp_ph3"
))

# 这里保持用 MetaCyc pathway ID 来筛（冒号前那段）
sig_pathways <- c(
  "PWY-6305",
  "PWY0-41",
  "PWY-7874",
  "PWY-6992",
  "P105-PWY",
  "PWY-7094",
  "PWY0-301",
  "P108-PWY"
)

# =============================
# HELPERS
# =============================
# 样本列名规范化：优先提取 F####，提取不到就原样
norm_sample <- function(x) {
  f <- str_extract(x, "F\\d+")
  ifelse(!is.na(f), f, x)
}

# 从 taxonomy feature 里抽 species 名：
# - 若含 s__Xxx_yyy 取 yyy
# - 否则直接用原字符串
extract_species <- function(x) {
  s <- str_extract(x, "s__[^|;\\s]+")
  s <- str_replace(s, "^s__", "")
  ifelse(!is.na(s), s, x)
}

# 保留到 | 之前的整段（含 "PWY-7237: xxx"）
extract_pathway_left <- function(x) {
  str_trim(str_replace(x, "\\|.*$", ""))
}

# 从 "PWY-7237: xxx" 提取 PWY-7237（冒号前）
extract_pathway_id <- function(x) {
  str_trim(str_extract(x, "^[^:]+"))
}

# 安全转数值（NA -> 0）
to_num0 <- function(v) {
  vv <- suppressWarnings(as.numeric(v))
  vv[is.na(vv)] <- 0
  vv
}

# =============================
# 1) READ taxonomy.tsv -> species x sample matrix
# =============================
cat("Reading taxonomy.tsv ...\n")
tax <- read_tsv(tax_file, show_col_types = FALSE)
names(tax)[1] <- "feature_raw"
names(tax)[-1] <- norm_sample(names(tax)[-1])

# 转数值
tax[-1] <- lapply(tax[-1], to_num0)

tax2 <- tax %>%
  mutate(species = extract_species(feature_raw)) %>%
  filter(species %in% sig_species) %>%
  select(species, where(is.numeric))

if (!("species" %in% names(tax2)) || nrow(tax2) == 0) {
  stop("taxonomy.tsv 里没有匹配到你给的显著菌名。请检查：taxonomy 第一列是否包含 s__ 物种名、或者物种拼写是否一致。")
}

# 如果同一物种出现多行，按样本求和
tax_mat <- tax2 %>%
  group_by(species) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

# =============================
# 2) READ HUMAnN3 stratified pathway table -> unstratified pathway x sample matrix
#    用 stratified 求和得到总通路丰度（对每个 pathway_id 把不同物种行加总）
#    同时保留 pathway_left（到 | 之前的完整描述）
# =============================
cat("Reading HUMAnN3 stratified pathways ...\n")
pwy <- read_tsv(pwy_file, comment = "", show_col_types = FALSE)
names(pwy)[1] <- "pathway_str"
names(pwy)[-1] <- norm_sample(names(pwy)[-1])
pwy[-1] <- lapply(pwy[-1], to_num0)

pwy_un <- pwy %>%
  filter(str_detect(pathway_str, "\\|")) %>%  # stratified 行
  mutate(
    pathway_left = extract_pathway_left(pathway_str),   # 你要保留的整段
    pathway_id   = extract_pathway_id(pathway_left)     # 用于筛选的 ID
  ) %>%
  filter(pathway_id %in% sig_pathways) %>%
  select(pathway_id, pathway_left, where(is.numeric)) %>%
  group_by(pathway_id, pathway_left) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

if (nrow(pwy_un) == 0) {
  stop("在 stratified 表里没找到你给的通路ID。请确认通路ID（如 PWY0-41）是否完全一致。")
}

# 生成 pathway 映射（用于输出保留完整描述）
pwy_map <- pwy_un %>% distinct(pathway_id, pathway_left)

# =============================
# 3) ALIGN SAMPLES
# =============================
sp_samples <- setdiff(names(tax_mat), "species")
pw_samples <- setdiff(names(pwy_un), c("pathway_id", "pathway_left"))
common_samples <- intersect(sp_samples, pw_samples)

if (length(common_samples) < 5) {
  stop(paste0("共同样本太少：", length(common_samples),
              "。请检查 taxonomy.tsv 与 pathway 表的样本列名是否能对齐（F#### 是否存在/一致）。"))
}

# 转成 sample 行，feature 列的宽表
sp_wide <- tax_mat %>%
  select(species, all_of(common_samples)) %>%
  pivot_longer(-species, names_to = "sample", values_to = "abund") %>%
  pivot_wider(names_from = species, values_from = abund) %>%
  arrange(sample)

pw_wide <- pwy_un %>%
  select(pathway_id, all_of(common_samples)) %>%
  pivot_longer(-pathway_id, names_to = "sample", values_to = "abund") %>%
  pivot_wider(names_from = pathway_id, values_from = abund) %>%
  arrange(sample)

dat <- sp_wide %>% inner_join(pw_wide, by = "sample")

# 可选：log(x + pseudo)（建议保留）
pseudo <- 1e-10
dat_t <- dat %>%
  mutate(across(all_of(c(sig_species, sig_pathways)), ~ log(.x + pseudo)))

# =============================
# 4) SPEARMAN (p only, no FDR)
# =============================
pairs <- expand.grid(
  species = sig_species,
  pathway_id = sig_pathways,
  stringsAsFactors = FALSE
)

calc_one <- function(sp, pw) {
  x <- dat_t[[sp]]
  y <- dat_t[[pw]]
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]

  # 变成常数就没法算相关
  if (length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(tibble(
      species = sp,
      pathway_id = pw,
      rho = NA_real_,
      p_value = NA_real_,
      n = sum(ok),
      n_nonzero_species = sum(dat[[sp]] > 0, na.rm = TRUE),
      n_nonzero_pathway = sum(dat[[pw]] > 0, na.rm = TRUE)
    ))
  }

  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(
    species = sp,
    pathway_id = pw,
    rho = unname(ct$estimate),
    p_value = ct$p.value,
    n = sum(ok),
    n_nonzero_species = sum(dat[[sp]] > 0, na.rm = TRUE),
    n_nonzero_pathway = sum(dat[[pw]] > 0, na.rm = TRUE)
  )
}

res <- purrr::pmap_dfr(pairs, ~ calc_one(..1, ..2)) %>%
  left_join(pwy_map, by = "pathway_id") %>%
  relocate(pathway_left, .after = pathway_id) %>%
  arrange(p_value, desc(abs(rho)))

write.csv(res, out_csv, row.names = FALSE)

cat("✅ 输出完成：", out_csv, "\n",
    "共同样本数：", length(common_samples), "\n",
    "菌数：", length(sig_species), "；通路数：", length(sig_pathways), "\n",
    "已保留通路完整描述（到 | 之前）：pathway_left\n",
    sep = "")
