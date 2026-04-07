library(tidyverse)
library(ggalluvial)

TOP_N   <- 40
IN_FILE <- "species_pathway_severity_merged.csv"
OUT_PNG <- "sankey_species_pathway_severity.png"

## 1) 读数据 + 清洗（只用这5列）
chains_clean <- readr::read_csv(
  IN_FILE,
  show_col_types = FALSE,
  col_types = cols(
    severity    = col_character(),
    species     = col_character(),
    pathway     = col_character(),
    rho_sp_pw   = col_double(),
    beta_pw_sev = col_double()
  )
) %>%
  transmute(
    species  = str_trim(species),
    pathway  = str_trim(pathway),
    clinical = str_trim(severity),
    rho_sp_pw   = as.numeric(rho_sp_pw),
    beta_pw_sev = as.numeric(beta_pw_sev),
    weight = abs(rho_sp_pw) * abs(beta_pw_sev)
  ) %>%
  filter(
    species  != "", pathway != "", clinical != "",
    is.finite(rho_sp_pw), is.finite(beta_pw_sev),
    is.finite(weight), weight > 0
  )

## 2) 聚合到 (species, pathway, clinical)
##    保留 rho/beta 的聚合值，用来判正负
chains_trip <- chains_clean %>%
  group_by(species, pathway, clinical) %>%
  summarise(
    weight = mean(weight, na.rm = TRUE),
    rho_sp_pw   = mean(rho_sp_pw, na.rm = TRUE),
    beta_pw_sev = mean(beta_pw_sev, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(chain_id = row_number()) %>%
  slice_max(order_by = weight, n = TOP_N)

## 3) 做长表（柱子用）
axis_levels <- c("species", "pathway", "clinical")
axis_labels <- c("Gut Microbiota", "Pathways", "Clinical scores")

alluvial_long <- chains_trip %>%
  select(chain_id, weight, species, pathway, clinical) %>%
  pivot_longer(
    cols = c(species, pathway, clinical),
    names_to = "axis",
    values_to = "stratum"
  ) %>%
  mutate(axis = factor(axis, levels = axis_levels, labels = axis_labels))

## 4) 柱子颜色（按所有 stratum 去重后调色）
n_stratum <- dplyr::n_distinct(alluvial_long$stratum)

stratum_palette <- colorRampPalette(
  c("#85add8", "#8bcd83", "#F6D55C", "#f9ba7a", "#d7a1c8", "#e7aaaa")
)(n_stratum)

stratum_cols <- alluvial_long %>%
  distinct(stratum) %>%
  arrange(stratum) %>%
  mutate(col = stratum_palette[dplyr::row_number()])

alluvial_long <- alluvial_long %>%
  left_join(stratum_cols, by = "stratum")

## 5) 两段流带数据分别准备 + 各用各的符号列上色
flow_cols <- chains_trip %>%
  transmute(
    chain_id,
    flow_col_sp_pw  = if_else(rho_sp_pw > 0, "#ce5e7e", "#5B8CFF"),   # rho_sp_pw: 正红负蓝
    flow_col_pw_sev = if_else(beta_pw_sev > 0, "#ce5e7e", "#5B8CFF")  # beta_pw_sev: 正红负蓝
  )

# A 段：species -> pathway
flow_sp_pw <- alluvial_long %>%
  filter(axis %in% axis_labels[1:2]) %>%   # 只保留前两列
  left_join(flow_cols, by = "chain_id")

# B 段：pathway -> clinical
flow_pw_sev <- alluvial_long %>%
  filter(axis %in% axis_labels[2:3]) %>%   # 只保留后两列
  left_join(flow_cols, by = "chain_id")

## 6) 画图：先画两段流带，再画柱子
g <- ggplot() +
  # species -> pathway 段：rho_sp_pw 决定红/蓝
  geom_flow(
    data = flow_sp_pw,
    aes(
      x = axis, stratum = stratum, alluvium = chain_id, y = weight,
      fill = flow_col_sp_pw, colour = flow_col_sp_pw
    ),
    alpha = 0.6, knot.pos = 0.4, show.legend = FALSE
  ) +
  # pathway -> clinical 段：beta_pw_sev 决定红/蓝
  geom_flow(
    data = flow_pw_sev,
    aes(
      x = axis, stratum = stratum, alluvium = chain_id, y = weight,
      fill = flow_col_pw_sev, colour = flow_col_pw_sev
    ),
    alpha = 0.6, knot.pos = 0.4, show.legend = FALSE
  ) +
  # 柱子：按 stratum palette（不受 flow 影响）
  geom_stratum(
    data = alluvial_long,
    aes(x = axis, stratum = stratum, y = weight, fill = col),
    colour = "black",
    width = 0.4
  ) +
  geom_text(
    data = alluvial_long,
    stat = "stratum",
    aes(x = axis, stratum = stratum, y = weight, label = after_stat(stratum)),
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

print(g)

ggsave(OUT_PNG, g, width = 11, height = 6, dpi = 600)
cat("✅ Saved: ", OUT_PNG, "\n", sep = "")
