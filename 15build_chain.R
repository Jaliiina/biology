suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

# =============================
# INPUTS
# =============================
f_sp_pw  <- "filtered_species_pathway_cor.csv"      # 菌-通路相关（含 rho, p_value）
f_pw_sev <- "filtered_pathway_severity_reg.csv"     # 通路-症状回归（含 beta, p）

# =============================
# OUTPUT
# =============================
out_chain <- "species_pathway_severity_merged.csv"

# =============================
# READ
# =============================
sp_pw <- read_csv(f_sp_pw, show_col_types = FALSE) %>%
  mutate(
    species = as.character(species),
    pathway_id = as.character(pathway_id),
    pathway_left = if ("pathway_left" %in% names(.)) as.character(pathway_left) else pathway_id,
    rho = as.numeric(rho),
    p_sp_pw = as.numeric(p_value),
    n_sp_pw = if ("n" %in% names(.)) as.integer(n) else NA_integer_
  ) %>%
  select(species, pathway_id, pathway_left, rho, p_sp_pw, n_sp_pw, everything())

pw_sev_raw <- read_csv(f_pw_sev, show_col_types = FALSE)

# 兼容不同列名：pathway_left / pathway / pathway_id
get_pathway_id <- function(x) str_trim(str_extract(as.character(x), "^[^:]+"))

pw_sev <- pw_sev_raw %>%
  mutate(severity = as.character(severity)) %>%
  {
    if ("pathway_id" %in% names(.)) {
      mutate(., pathway_id = as.character(pathway_id),
                pathway_left = if ("pathway_left" %in% names(.)) as.character(pathway_left) else
                              if ("pathway" %in% names(.)) as.character(pathway) else as.character(pathway_id))
    } else if ("pathway_left" %in% names(.)) {
      mutate(., pathway_left = as.character(pathway_left),
                pathway_id = get_pathway_id(pathway_left))
    } else if ("pathway" %in% names(.)) {
      mutate(., pathway_left = as.character(pathway),
                pathway_id = get_pathway_id(pathway_left))
    } else {
      stop("filtered_pathway_severity_reg.csv 中未找到 pathway / pathway_left / pathway_id 列。")
    }
  } %>%
  mutate(
    beta_pw_sev = as.numeric(beta),
    p_pw_sev    = as.numeric(p),
    n_pw_sev    = if ("n" %in% names(.)) as.integer(n) else NA_integer_
  )

# =============================
# MERGE: species-pathway  ×  pathway-severity  =>  species-pathway-severity
# 不管方向是否违背：不计算 coherence，直接做笛卡尔式扩展（同一通路可对应多个 severity）
# =============================
chain <- sp_pw %>%
  inner_join(
    pw_sev %>% select(severity, pathway_id, pathway_left, beta_pw_sev, p_pw_sev, n_pw_sev, everything()),
    by = c("pathway_id", "pathway_left")
  ) %>%
  # 输出更清爽的核心列（其余列你也可以保留）
  select(
    severity,
    species,
    pathway_id,
    pathway_left,
    rho,
    p_sp_pw,
    n_sp_pw,
    beta_pw_sev,
    p_pw_sev,
    n_pw_sev,
    everything()
  ) %>%
  arrange(severity, p_pw_sev, p_sp_pw, desc(abs(rho)))

write.csv(chain, out_chain, row.names = FALSE)

cat("✅ 已合并输出：", out_chain, "\n",
    "   行数：", nrow(chain), "\n",
    "   物种数：", n_distinct(chain$species), "\n",
    "   通路数：", n_distinct(chain$pathway_id), "\n",
    "   severity 数：", n_distinct(chain$severity), "\n",
    sep = "")
