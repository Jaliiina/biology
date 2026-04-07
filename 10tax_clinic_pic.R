## ====== Species–Severity regression Circos plot (p<0.05) ======
library(tidyverse)
library(circlize)

# 1) 读数据 ---------------------------------------------------------
in_file <- "regression_species_multi_severity_sig_p005.csv"
df_raw <- read.csv(in_file, stringsAsFactors = FALSE, check.names = FALSE)

# 小工具：自动匹配列名
pick_col <- function(candidates, cols) {
  hit <- intersect(candidates, cols)
  if (length(hit) == 0) stop("找不到列：", paste(candidates, collapse = " / "))
  hit[1]
}

cols <- colnames(df_raw)

col_trait <- pick_col(c("severity", "trait", "clinical", "phenotype"), cols)
col_sp    <- pick_col(c("species", "species_label", "taxon"), cols)
col_beta  <- pick_col(c("beta", "coef", "estimate"), cols)
col_p     <- pick_col(c("p", "p_value", "P"), cols)

# 你关心的 severity 顺序（按你之前的顺序）
clinical_vars <- c("ATEC-4","ATEC-3","ATEC-2","ATEC-1","ATEC-Total","ABC","CARS")

# 2) 筛选：severity 范围 + p + |beta| -------------------------------
P_CUT <- 0.05
B_CUT <- 0.00   # 如果你想去掉很小的效应，改成比如 0.1/0.2

df <- df_raw %>%
  transmute(
    species = .data[[col_sp]],
    trait   = .data[[col_trait]],
    beta    = as.numeric(.data[[col_beta]]),
    p       = as.numeric(.data[[col_p]])
  ) %>%
  filter(
    trait %in% clinical_vars,
    !is.na(beta), !is.na(p),
    p < P_CUT,
    abs(beta) >= B_CUT
  )

if (nrow(df) == 0) {
  stop("在 p <", P_CUT, " 且 |beta| >=", B_CUT, " 条件下，没有任何显著的菌种–严重程度关联。")
}

# 3) 构建边数据：from = species, to = severity -----------------------
edges <- df %>%
  transmute(
    from     = species,
    to       = trait,
    beta     = beta,
    abs_beta = abs(beta)
  )

# 4) 节点集合 --------------------------------------------------------
species_nodes <- unique(edges$from)

# severity 节点按 clinical_vars 的顺序，但只保留真实出现过的
symptom_nodes <- clinical_vars[clinical_vars %in% unique(edges$to)]

# 扇区顺序：物种在前，症状在后
all_sectors <- c(species_nodes, symptom_nodes)

# 锁定顺序
edges <- edges %>%
  mutate(
    from = factor(from, levels = all_sectors),
    to   = factor(to,   levels = all_sectors)
  )

# 节点颜色：菌种一色、症状一色（沿用你参考的配色逻辑）
grid.col <- c(
  setNames(rep("#4D9221", length(species_nodes)), species_nodes),  # species: 绿
  setNames(rep("#525252", length(symptom_nodes)), symptom_nodes)   # symptom: 深灰
)

# 边颜色：按 beta 渐变（蓝负红正）
# 为了适配你的数据范围，用对称的最大绝对值做映射
max_abs_beta <- max(abs(edges$beta), na.rm = TRUE)
if (!is.finite(max_abs_beta) || max_abs_beta == 0) max_abs_beta <- 1

col_fun <- colorRamp2(
  c(-max_abs_beta, 0, max_abs_beta),
  c("#4575b4", "grey90", "#d73027")
)
edge_col <- col_fun(edges$beta)

# 5) 画 Circos --------------------------------------------------------
circos.clear()

gap_after <- c(
  rep(2, length(species_nodes) - 1), 10,      # 物种之间 2 度，物种组后面 10 度
  rep(2, length(symptom_nodes) - 1), 10       # 症状之间 2 度，最后再留一点空
)

circos.par(
  start.degree = 90,
  gap.after = gap_after,
  track.margin = c(0.01, 0.01)
)

out_png <- "species_severity_regression_circos_p005.png"
png(out_png, width = 2400, height = 2400, res = 300)

circos.clear()
circos.par(
  start.degree = 90,
  gap.after = gap_after,
  track.margin = c(0.01, 0.01)
)

# chordDiagram：权重用 abs_beta
chordDiagram(
  x = edges[, c("from", "to", "abs_beta")],
  grid.col = grid.col,
  col = edge_col,
  transparency = 0.25,
  directional = 0,
  order = all_sectors,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1)
)

# 外圈标签
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    circos.text(
      x = mean(xlim),
      y = 0,
      labels = sector.name,
      facing = "bending.out",
      niceFacing = TRUE,
      cex = 1.0,
      font = 1.0
    )
  }
)

title("Associations Between Microbial Species and ASD Severity")

dev.off()

cat("✅ Circos 图已保存：", out_png, "\n", sep = "")
