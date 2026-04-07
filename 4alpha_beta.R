## ----------------------------
## 0) packages
## ----------------------------
pkgs <- c("tidyverse", "readr", "vegan", "patchwork")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

## ----------------------------
## 1) input paths
## ----------------------------
meta_file <- "metadata.txt"
tax_file  <- "taxonomy.tsv"

## ----------------------------
## 2) read metadata
## ----------------------------
meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  rename(sample_id = 1, group = 2) %>%
  mutate(group = factor(group, levels = c("Low", "High")))

## ----------------------------
## 3) read taxonomy table
##    - first column: taxonomy ID string
##    - remaining columns: samples
## ----------------------------
tax_raw <- read_tsv(tax_file, show_col_types = FALSE)

# ensure first column is named ID
if (!("ID" %in% names(tax_raw))) names(tax_raw)[1] <- "ID"

## ----------------------------
## 4) collapse to species level (s__)
## ----------------------------
extract_species <- function(x) {
  # capture s__... until next | (ignore trailing t__ etc.)
  stringr::str_extract(x, "s__[^|]+")
}

tax_sp <- tax_raw %>%
  mutate(species = extract_species(ID)) %>%
  filter(!is.na(species)) %>%
  select(-ID)

# convert sample columns to numeric
sample_cols <- setdiff(colnames(tax_sp), "species")
tax_sp[sample_cols] <- lapply(tax_sp[sample_cols], function(x) as.numeric(as.character(x)))

# sum rows that map to the same species (e.g., multiple t__ under same s__)
tax_sp_sum <- tax_sp %>%
  group_by(species) %>%
  summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

## ----------------------------
## 5) make matrix: samples x species
## ----------------------------
otu_mat <- tax_sp_sum %>%
  column_to_rownames("species") %>%
  as.matrix()

# transpose: species x samples -> samples x species
otu_mat <- t(otu_mat)

# align samples with metadata
common_ids <- intersect(rownames(otu_mat), meta$sample_id)
otu_mat <- otu_mat[common_ids, , drop = FALSE]
meta    <- meta %>% filter(sample_id %in% common_ids)

# row-normalize to relative abundance (even if already %)
otu_rel <- sweep(otu_mat, 1, rowSums(otu_mat), FUN = "/")
otu_rel[is.na(otu_rel)] <- 0

## ----------------------------
## 6) alpha diversity
## ----------------------------
alpha_df <- tibble(
  sample_id = rownames(otu_rel),
  Observed  = vegan::specnumber(otu_rel > 0),                 # richness: count of present species
  Shannon   = vegan::diversity(otu_rel, index = "shannon")    # Shannon
) %>%
  left_join(meta, by = "sample_id") %>%
  mutate(group = factor(group, levels = c("Low", "High")))

# Wilcoxon tests (to match your aesthetics)
wil_shannon  <- wilcox.test(Shannon ~ group, data = alpha_df)
wil_observed <- wilcox.test(Observed ~ group, data = alpha_df)

## ----------------------------
## 7) beta diversity: Bray-Curtis + PCoA + PERMANOVA
## ----------------------------
bray <- vegan::vegdist(otu_rel, method = "bray")

# PCoA (classical MDS)
pcoa <- cmdscale(bray, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$sample_id <- rownames(pcoa_df)

# variance explained (positive eigenvalues)
var_expl <- pcoa$eig[pcoa$eig > 0] / sum(pcoa$eig[pcoa$eig > 0])
pc1_lab <- paste0("PC1 (", round(var_expl[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_expl[2] * 100, 1), "%)")

pcoa_df <- pcoa_df %>%
  left_join(meta, by = "sample_id") %>%
  mutate(group = factor(group, levels = c("Low", "High")))

# PERMANOVA
set.seed(123)
adon <- vegan::adonis2(bray ~ group, data = meta, permutations = 999)
pval <- adon$`Pr(>F)`[1]
r2   <- adon$R2[1]

## ----------------------------
## 8) plotting (match your reference aesthetics)
## ----------------------------

# A: Shannon
p_alpha_shannon <- ggplot(alpha_df, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  labs(
    title = "Shannon diversity index",
    x = "Calcium level",
    y = "Shannon index"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(alpha_df$Shannon, na.rm = TRUE) * 1.03,
    label = paste0("Wilcoxon p = ", signif(wil_shannon$p.value, 3))
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# B: Observed richness
p_alpha_observed <- ggplot(alpha_df, aes(x = group, y = Observed, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  labs(
    title = "Observed species richness",
    x = "Calcium level",
    y = "Number of species"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(alpha_df$Observed, na.rm = TRUE) * 1.03,
    label = paste0("Wilcoxon p = ", signif(wil_observed$p.value, 3))
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# C: PCoA
p_beta <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.9) +
  labs(
    title = "PCoA (Bray–Curtis)",
    x = pc1_lab,
    y = pc2_lab,
    caption = paste0(
      "PERMANOVA: p = ",
      signif(pval, 3),
      ", R² = ",
      signif(r2, 3)
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

## ----------------------------
## 9) combine into one figure with tags A/B/C
## ----------------------------
design <- "
AB
CC
"

final_plot <- p_alpha_shannon + p_alpha_observed + p_beta +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A")

print(final_plot)

ggsave(
  "alpha_beta_diversity_calcium.png",
  final_plot,
  width = 10,
  height = 8,
  dpi = 300
)
