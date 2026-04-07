# ## 安装 & 加载依赖包
# pkgs <- c("dplyr", "ggplot2", "forcats", "readr")
# to_install <- setdiff(pkgs, installed.packages()[,1])
# if (length(to_install)) install.packages(to_install)

# library(dplyr)
# library(ggplot2)
# library(forcats)
# library(readr)

# ## 1. 读入结果文件（就是你刚给的这个）
# res_file <- "lefse_pathway_stratified_diff.csv"
# lefse_df <- read_csv(res_file, show_col_types = FALSE)

# ## 2. 去掉指定通路：COA-PWY: coenzyme A biosynthesis I (prokaryotic)
# ## 数据里是 Taxa 列，后面大概率会带一个 |unclassified，所以用 grepl 匹配前半段
# lefse_df <- lefse_df %>%
#   filter(
#     !grepl(
#       "COA-PWY: coenzyme A biosynthesis I \\(prokaryotic\\)",
#       Taxa
#     )
#   )

# ## 3. 准备 FDR 和 -log10(FDR)
# ## 如果你就用表里的 P.adj，当 FDR 用：
# lefse_df <- lefse_df %>%
#   mutate(
#     FDR = P.adj
#   )

# ## 如果你想自己做 BH 校正，可以改成：
# # lefse_df <- lefse_df %>%
# #   mutate(
# #     FDR = p.adjust(P.unadj, method = "BH")
# #   )

# lefse_df <- lefse_df %>%
#   mutate(
#     neg_log10_FDR = -log10(FDR)
#   )

# ## 4. 取 top30，按 |LDA| 最大
# lefse_top30 <- lefse_df %>%
#   arrange(desc(abs(LDA))) %>%
#   slice_head(n = 30) %>%
#   ## y 轴用通路名；如果你嫌后面的 |unclassified 短得难看，可以只保留前半段：
#   mutate(
#     pathway = sub("\\|.*$", "", Taxa),     # 只保留 "PWY-XXXX: ..." 这部分
#     pathway = fct_reorder(pathway, LDA)    # 按 LDA 排序，方便展示
#   )

# ## 5. 画气泡图
# p <- ggplot(lefse_top30,
#             aes(x = LDA, y = pathway)) +
#   geom_point(
#     aes(
#       size  = neg_log10_FDR,
#       color = Group
#     ),
#     alpha = 0.8
#   ) +
#   scale_color_manual(
#     values = c(
#       "Low"  = "#FFA500",
#       "High" = "#1F77B4"
#     ),
#     name = "Group"
#   ) +
#   scale_x_continuous(
#     name   = "LDA score",
#     breaks = 0:3,
#     limits = c(0, 3)
#   ) +
#   scale_size_continuous(
#     name = "-log10(FDR)"
#   ) +
#   ylab("Pathway") +
#   labs(title = "High vs. Low Calcium Groups — Metabolic Pathway Differences") +
#   theme_bw() +
#   theme(
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_line(),
#     panel.grid.major.y = element_blank(),
#     axis.text.y = element_text(size = 8),
#     axis.text.x = element_text(size = 10),
#     axis.title  = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
#     legend.position = "right"
#   )

# print(p)


# print(p)

# ## 6. 保存成图片
# ggsave("pathways_dotplot_LEfSe_top30_bubble.png",
#        p, width = 9, height = 6, dpi = 300)

## 安装 & 加载依赖包
pkgs <- c("dplyr", "ggplot2", "forcats", "readr")
to_install <- setdiff(pkgs, installed.packages()[,1])
if (length(to_install)) install.packages(to_install)

library(dplyr)
library(ggplot2)
library(forcats)
library(readr)

## 1. 读入结果文件
res_file <- "lefse_pathway_high_vs_low.csv"
lefse_df <- read_csv(res_file, show_col_types = FALSE)

# ## 2. 去掉指定通路：COA-PWY: coenzyme A biosynthesis I (prokaryotic)
# lefse_df <- lefse_df %>%
#   filter(
#     !grepl(
#       "COA-PWY: coenzyme A biosynthesis I \\(prokaryotic\\)",
#       Taxa
#     )
#   )

## 3. 准备 FDR 和 -log10(FDR)
lefse_df <- lefse_df %>%
  mutate(
    FDR = P.adj,
    neg_log10_FDR = -log10(FDR)
  )

## 4. 取 top30，按 |LDA| 最大，并做坐标变换 -------------------------
lefse_top30 <- lefse_df %>%
  arrange(desc(abs(LDA))) %>%
  slice_head(n = 10) %>%
  mutate(
    pathway = sub("\\|.*$", "", Taxa),   # 只保留 "PWY-XXXX: ..." 这部分
    # 把 LDA 加上方向：High 为正，Low 为负
    signed_lda = ifelse(tolower(Group) == "high", LDA, -LDA),
    pathway = fct_reorder(pathway, signed_lda)  # 按带符号 LDA 排序
  )

## 5. 画泡泡棒棒糖图（lollipop + bubble） --------------------------
p <- ggplot(lefse_top30,
            aes(x = signed_lda, y = pathway)) +
  # 棒棒糖的“棒”：从 0 连到 LDA
  geom_segment(
    aes(x = 0, xend = signed_lda, y = pathway, yend = pathway,
        color = Group),
    linewidth = 0.6
  ) +
  # 棒棒糖的“糖”：末端的泡泡，大小代表 -log10(FDR)
  geom_point(
    aes(
      size  = neg_log10_FDR,
      color = Group
    ),
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "Low"  = "#FFA500",  # 低钙组：橙色
      "High" = "#1F77B4"   # 高钙组：蓝色
    ),
    name = "Group"
  ) +
  scale_size_continuous(
    name = "-log10(FDR)"
  ) +
  scale_x_continuous(
    name   = "LDA score",
    # 让坐标左右稍微留点空，不强行限定到 0~3
    breaks = pretty(range(lefse_top30$signed_lda))
  ) +
  ylab("Pathway") +
  labs(title = "Metabolic Pathway Differences between High and Low Calcium Groups") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    axis.title  = element_text(size = 10),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "right"
  )

print(p)

## 6. 保存成图片
ggsave("pathways_lollipop_LEfSe_top30.png",
       p, width = 9, height = 6, dpi = 600)
