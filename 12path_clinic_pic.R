library(tidyverse)

# 1) 读入“通路 ~ 严重程度”的回归显著结果 -----------------------------
# 你这一步的输出文件（p<0.05 的那份）
df <- read.csv("regression_pathway_multi_severity_sig_p005.csv",
               stringsAsFactors = FALSE,
               check.names = FALSE)

# 你关心的症状变量顺序（横轴）
clinical_vars <- c("ATEC-Total", "ATEC-1", "ATEC-2",
                   "ATEC-3", "ATEC-4", "ABC", "CARS")

# p 阈值、最多展示多少条通路
p_cut <- 0.05
N_TOP <- 10  # 你只有 8 条的话可以设 10，最终也只会画 8 条；想固定 8 条就改 8

# 2) 清理：通路名 + 过滤症状 + 去掉不想要的通路 ------------------------
# 你的回归结果里 pathway 一般已经是 base pathway（无 |g__...），但我仍做一次保险处理
df2 <- df %>%
  mutate(
    pathway = sub("\\|.*$", "", pathway),
    trait   = severity
  ) %>%
  filter(trait %in% clinical_vars) %>%
  # 如果你也想剔除某些通路前缀，按需要开关（和你原代码一致风格）
  filter(!grepl("^COA-PWY:", pathway)) %>%
  filter(!grepl("UNINTEGRATED", pathway, ignore.case = TRUE))

# 只看 p < 阈值（你要求的）
df_sig <- df2 %>% filter(p < p_cut)

if (nrow(df_sig) == 0) {
  stop("在 p <", p_cut, " 下没有任何显著的通路–严重程度关联。")
}

# 3) 用 max |beta| 选出 Top N 通路（逻辑对齐你的 max |rho|） -----------
top_pathways <- df_sig %>%
  group_by(pathway) %>%
  summarise(max_abs_beta = max(abs(beta), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_beta)) %>%
  slice_head(n = N_TOP) %>%
  pull(pathway)

df_plot <- df_sig %>% filter(pathway %in% top_pathways)

# y 轴通路顺序：按 max_abs_beta 从小到大（越强越靠上）
pathway_order <- df_plot %>%
  group_by(pathway) %>%
  summarise(max_abs_beta = max(abs(beta), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_abs_beta) %>%
  pull(pathway)

df_plot <- df_plot %>%
  mutate(
    pathway = factor(pathway, levels = pathway_order),
    trait   = factor(trait, levels = clinical_vars),
    neg_log10_p = -log10(p)
  )

# 4) 画 dot plot：颜色 = beta，大小 = -log10(p) ------------------------
p <- ggplot(df_plot, aes(x = pathway, y = trait)) +
  geom_point(
    aes(size = neg_log10_p, color = beta),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low  = "#4575b4",  # 负相关（蓝）
    mid  = "white",
    high = "#d73027",  # 正相关（红）
    limits = c(-0.5, 0.5),  
    oob = scales::squish,    
    midpoint = 0,
    name = "Regression beta"
  ) +
  scale_size_continuous(
    name = "-log10(p)",
    range = c(1.5, 6)
  ) +
  xlab("Metabolic pathways") +
  ylab("Clinical severity scores") +
  labs(title = "Associations Between Metabolic Pathway and ASD Severity") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 7),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 80),
    legend.position = "right"
  )

# 5) 保存 + 输出 -------------------------------------------------------
ggsave("pathway_clinical_dotplot_topN.png",
       p, width = 8.5, height = 6, dpi = 600)

print(p)
