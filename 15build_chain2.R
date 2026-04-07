suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

# =============================
# INPUTS
# =============================
f_sp_pw  <- "species_pathway_spearman_ponly.csv"
f_sp_sev <- "regression_species_multi_severity_sig_p005.csv"
f_pw_sev <- "regression_pathway_multi_severity_sig_p005.csv"

# =============================
# OUTPUTS
# =============================
out_sp_pw <- "filtered_species_pathway_cor.csv"      # 菌-通路（只保留三表共同的菌/通路）
out_pw_sev <- "filtered_pathway_severity_reg.csv"    # 通路-症状（只保留三表共同的通路）

# =============================
# HELPERS
# =============================
# 从 "PWY-6305: xxx" 里提取 ID（PWY-6305 / P105-PWY 等）
get_pathway_id <- function(x) {
  x <- as.character(x)
  str_trim(str_extract(x, "^[^:]+"))
}

# =============================
# READ
# =============================
sp_pw <- read_csv(f_sp_pw, show_col_types = FALSE) %>%
  mutate(
    species = as.character(species),
    pathway_id = as.character(pathway_id)
  )

sp_sev <- read_csv(f_sp_sev, show_col_types = FALSE) %>%
  mutate(
    severity = as.character(severity),
    species  = as.character(species)
  )

pw_sev_raw <- read_csv(f_pw_sev, show_col_types = FALSE)

# pw~sev 这张表有时列名叫 pathway，有时叫 pathway_left/pathway_id
# 这里做一个兼容：优先找 pathway_id；否则从 pathway（全称）里提取
if ("pathway_id" %in% names(pw_sev_raw)) {
  pw_sev <- pw_sev_raw %>%
    mutate(
      severity = as.character(severity),
      pathway_id = as.character(pathway_id),
      pathway_left = if ("pathway_left" %in% names(.)) as.character(pathway_left) else pathway_id
    )
} else if ("pathway" %in% names(pw_sev_raw)) {
  pw_sev <- pw_sev_raw %>%
    mutate(
      severity = as.character(severity),
      pathway_left = as.character(pathway),
      pathway_id = get_pathway_id(pathway_left)
    )
} else if ("pathway_left" %in% names(pw_sev_raw)) {
  pw_sev <- pw_sev_raw %>%
    mutate(
      severity = as.character(severity),
      pathway_left = as.character(pathway_left),
      pathway_id = get_pathway_id(pathway_left)
    )
} else {
  stop("regression_pathway_multi_severity_sig_p005.csv 中未找到 pathway / pathway_left / pathway_id 列。")
}

# =============================
# 1) 找“三表都存在”的菌 & 通路
# =============================
common_species <- intersect(unique(sp_sev$species), unique(sp_pw$species))
common_pathways <- intersect(unique(pw_sev$pathway_id), unique(sp_pw$pathway_id))

cat("Common species:", length(common_species), "\n")
cat("Common pathways:", length(common_pathways), "\n")

# =============================
# 2) 过滤并保存：菌-通路相关（只保留共同菌/共同通路）
# =============================
sp_pw_f <- sp_pw %>%
  filter(
    species %in% common_species,
    pathway_id %in% common_pathways
  ) %>%
  arrange(p_value, desc(abs(rho)))

write.csv(sp_pw_f, out_sp_pw, row.names = FALSE)

# =============================
# 3) 过滤并保存：通路-症状回归（只保留共同通路）
# =============================
pw_sev_f <- pw_sev %>%
  filter(pathway_id %in% common_pathways) %>%
  arrange(severity, p)

write.csv(pw_sev_f, out_pw_sev, row.names = FALSE)

cat("✅ Saved:\n",
    "  - ", out_sp_pw, " (rows=", nrow(sp_pw_f), ")\n",
    "  - ", out_pw_sev, " (rows=", nrow(pw_sev_f), ")\n",
    sep = "")
