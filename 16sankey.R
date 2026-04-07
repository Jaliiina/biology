library(tidyverse)
library(ggalluvial)
library(stringr)

## ============================
## 0. 参数
## ============================
TOP_N   <- 40
IN_FILE <- "species_pathway_severity_merged.csv"
OUT_PNG <- "sankey_species_pathway_severity_hd.png"

## 1. 读数据 ----------------------------------------------------
chains_raw <- readr::read_csv(IN_FILE, show_col_types = FALSE)

## 2. 清理/统一字段（基于我们自己的数据集） ----------------------
# 目标字段：species, pathway_feature(可含描述), trait(severity), weight
chains_clean <- chains_raw %>%
  mutate(
    species_label = as.character(species),
    trait = as.character(severity),

    # 通路字段优先用 pathway_left（带描述），否则 fallback
    pathway_feature = case_when(
      "pathway_left" %in% names(.) ~ as.character(pathway_left),
      "pathway"      %in% names(.) ~ as.character(pathway),
      "pathway_id"   %in% names(.) ~ as.character(pathway_id),
      TRUE ~ NA_character_
    ),

    # 去掉 | 后面的（保险）
    pathway_clean = str_remove(pathway_feature, "\\s*\\|.*$"),

    # 取我们链路强度（你也可以换成 abs(rho)）
    rho = as.numeric(rho),
    beta_pw_sev = case_when(
      "beta_pw_sev" %in% names(.) ~ as.numeric(beta_pw_sev),
      "beta"        %in% names(.) ~ as.numeric(beta),
      TRUE ~ NA_real_
    ),

    # weight：链条强度（建议 = |rho|*|beta|）
    weight = abs(rho) * abs(beta_pw_sev)
  ) %>%
  filter(!is.na(species_label), !is.na(pathway_clean), !is.na(trait)) %>%
  filter(is.finite(weight), weight > 0)

## 3. 在 (species, pathway, trait) 级别重新聚合 -----------------
chains_trip <- chains_clean %>%
  group_by(species_label, pathway_clean, trait) %>%
  summarise(
    species  = first(species_label),
    pathway  = first(pathway_clean),
    clinical = first(trait),

    # 聚合权重
    weight = mean(weight, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    chain_id = row_number()
  )

## 3.5 选择 Top N 链 --------------------------------------------
chains_trip <- chains_trip %>%
  slice_max(order_by = weight, n = TOP_N)

## 4. 转成长表 -------------------------------------------------
alluvial_long <- chains_trip %>%
  select(chain_id, weight, species, pathway, clinical) %>%
  pivot_longer(
    cols      = c(species, pathway, clinical),
    names_to  = "axis",
    values_to = "stratum"
  ) %>%
  mutate(
    axis = factor(
      axis,
      levels = c("species", "pathway", "clinical"),
      labels = c("Gut Microbiota", "Pathways", "Clinical scores")
    )
  )

## 5. 分配颜色（完全照抄你给的美工：按 stratum 上色） ----------
n_stratum <- n_distinct(alluvial_long$stratum)

stratum_palette <- colorRampPalette(
  c("#618ebd", "#8bcd83", "#F6D55C", "#f9ba7a",
    "#d7a1c8", "#bbb9b9")
)(n_stratum)

stratum_cols <- alluvial_long %>%
  distinct(stratum) %>%
  arrange(stratum) %>%
  mutate(col = stratum_palette)

alluvial_long <- alluvial_long %>%
  left_join(stratum_cols, by = "stratum")

## 6. 画桑基图 --------------------------------------------------
g <- ggplot(
  alluvial_long,
  aes(
    x        = axis,
    stratum  = stratum,
    alluvium = chain_id,
    y        = weight
  )
) +
  geom_flow(
    aes(fill = col, colour = col),
    alpha       = 0.6,
    knot.pos    = 0.4,
    show.legend = FALSE
  ) +
  geom_stratum(
    aes(fill = col),
    colour = "black",
    width  = 0.4
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3,
    colour = "black"
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "Microbiota–Pathway–Clinical Associations",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.x      = element_blank(),
    axis.title.x     = element_blank(),
    axis.title.y     = element_blank(),
    axis.text.x.top  = element_text(size = 12, face = "bold", vjust = 3),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  scale_x_discrete(position = "top")

g

## 7. 保存 -----------------------------------------------------
ggsave(
  OUT_PNG,
  g, width = 11, height = 6, dpi = 600
)

cat("✅ Saved: ", OUT_PNG, "\n", sep = "")
