suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# =============================
# INPUT
# =============================
in_csv  <- "species_pathway_spearman_ponly.csv"
out_png <- "species_pathway_heatmap_strong.png"

# =============================
# 1) 读入 & 基础处理
# =============================
df <- read_csv(in_csv, show_col_types = FALSE) %>%
  mutate(
    species = as.character(species),
    pathway_left = as.character(pathway_left),
    rho = as.numeric(rho),
    p_value = as.numeric(p_value),
    sig = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# =============================
# 2) 排序（让图更好看）
#    - species：按最大 |rho|
#    - pathway：按最小 p 值
# =============================
species_order <- df %>%
  group_by(species) %>%
  summarise(score = max(abs(rho), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(score)) %>%
  pull(species)

pathway_order <- df %>%
  group_by(pathway_left) %>%
  summarise(score = min(p_value, na.rm = TRUE), .groups = "drop") %>%
  arrange(score) %>%
  pull(pathway_left)

df <- df %>%
  mutate(
    species = factor(species, levels = species_order),
    pathway_left = factor(pathway_left, levels = pathway_order)
  )

# =============================
# 3) 画热图（高对比配色）
# =============================
p <- ggplot(df, aes(x = species, y = pathway_left, fill = rho)) +
  geom_tile(color = "grey85", linewidth = 0.4) +
  geom_text(
    aes(label = sig),
    size = 5,
    color = "black",
    na.rm = TRUE
  ) +
  scale_fill_gradient2(
    low = "#4575b4",     # 深蓝（强负相关）
    mid = "white",
    high = "#d73027",    # 深红（强正相关）
    midpoint = 0,
    limits = c(-0.35, 0.35),   # 🔥 收紧范围，颜色立刻变重
    name = "Spearman ρ"
  ) +
  labs(
    title = "Associations between symptom-related gut microbial species and metabolic pathways",
    x = "Species",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(
    hjust = 0.5,    # 居中
    face = "bold",  # 加粗
    size = 14       # 字号（可调）
  )
  )
print(p)

ggsave(out_png, p, width = 12, height = 6.5, dpi = 300)
cat("✅ 热图已保存：", out_png, "\n", sep = "")
