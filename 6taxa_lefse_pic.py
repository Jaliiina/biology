import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import re

# -----------------------------
# Config
# -----------------------------
in_csv = "lefse_species_high_vs_low.csv"
out_png = "lefse_species_LDA_gt1.5_diverging_clean.png"

HIGH_COLOR = "#2F4B7C"  # High (blue)
LOW_COLOR  = "#C44E52"  # Low (red)
LDA_CUTOFF = 2.0

# -----------------------------
# Load
# -----------------------------
df = pd.read_csv(in_csv)

df["Taxa"] = df["Taxa"].astype(str)
df["Group"] = df["Group"].astype(str).str.strip().str.title()
df["LDA"] = pd.to_numeric(df["LDA"], errors="coerce")

# keep High / Low only
df = df[df["Group"].str.lower().isin(["high", "low"])].copy()
df = df[df["LDA"].notna()].copy()

# -----------------------------
# Map to species level (keep t__, collapse to s__)
# -----------------------------
df["SpeciesTag"] = df["Taxa"].str.extract(r"(s__[^|]+)", expand=False)
df = df[df["SpeciesTag"].notna()].copy()

df["Species"] = df["SpeciesTag"].str.replace("s__", "", regex=False)

# -----------------------------
# Remove placeholder species (GGB / SGB)
# -----------------------------
# 规则：species 名以 GGB 或 SGB 开头，或包含 _SGBxxxx → 删除
df = df[
    ~df["Species"].str.match(r"^(GGB|SGB)", case=False) &
    ~df["Species"].str.contains(r"_SGB\d+", regex=True)
].copy()

# -----------------------------
# LDA cutoff
# -----------------------------
df = df[df["LDA"] > LDA_CUTOFF].copy()

# -----------------------------
# Collapse duplicated species
# keep row with max LDA per species
# -----------------------------
df = (
    df.sort_values(["Species", "LDA"], ascending=[True, False])
      .drop_duplicates(subset=["Species"], keep="first")
)

# -----------------------------
# Diverging values + group-wise sorting
# -----------------------------
high = df[df["Group"] == "High"].sort_values("LDA", ascending=False).copy()
low  = df[df["Group"] == "Low"].sort_values("LDA", ascending=False).copy()

high["LDA_signed"] = high["LDA"]
low["LDA_signed"]  = -low["LDA"]

plot_df = pd.concat([high, low], axis=0).reset_index(drop=True)

if plot_df.empty:
    raise ValueError("No species left after filtering (LDA > 1.5, remove GGB/SGB).")

labels = plot_df["Species"].tolist()
vals   = plot_df["LDA_signed"].tolist()
colors = [HIGH_COLOR if g == "High" else LOW_COLOR for g in plot_df["Group"]]

# -----------------------------
# Plot
# -----------------------------
fig_h = 0.40 * len(plot_df) + 2.0
fig, ax = plt.subplots(figsize=(8, fig_h))

y = list(range(len(plot_df)))
ax.barh(y, vals, color=colors, edgecolor="white", linewidth=0.6, height=0.82)

ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9)
ax.invert_yaxis()

# center line
ax.axvline(0, color="gray", linestyle="--", linewidth=1)

ax.set_title("Species Differences between High and Low Calcium Groups", fontsize=11)
ax.set_xlabel("LDA Score", fontsize=10)
ax.set_ylabel("Feature", fontsize=10)

# symmetric limits
max_abs = max(abs(min(vals)), abs(max(vals)))
ax.set_xlim(-max_abs * 1.08, max_abs * 1.08)

ax.grid(False)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

legend_handles = [
    Patch(facecolor=HIGH_COLOR, label="High"),
    Patch(facecolor=LOW_COLOR, label="Low"),
]
ax.legend(handles=legend_handles, frameon=False, loc="lower right")

plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()

print(f"Saved: {out_png}")
print(f"High: {len(high)} species, Low: {len(low)} species, Total: {len(plot_df)} species")
