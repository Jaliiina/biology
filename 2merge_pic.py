import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

# 读取 partial Spearman 结果
df_sp = pd.read_csv("Ca_vs_scales_partialSpearman_ageAdjusted.csv").set_index("量表")

# 4 个显著量表
targets = ["CARS", "ATEC-Total", "ATEC-1", "ATEC-3"]
# titles  = ["CARS", "ATEC-Total", "ATEC-1", "ATEC-3"]
labels  = ["A", "B", "C", "D"]

titles = [
    "Association between Calcium and CARS",
    "Association between Calcium and ATEC-Total",
    "Association between Calcium and ATEC-1",
    "Association between Calcium and ATEC-3"
]


# 读取数据
df = pd.read_excel("data0116.xlsx")

# --- 画图开始 ---
fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=300)
axes = axes.flatten()

for i, scale in enumerate(targets):
    ax = axes[i]

    # 散点
    ax.scatter(df["Calcium"], df[scale], s=35, alpha=0.7, color="#4C72B0")

    # LOWESS
    lowess = sm.nonparametric.lowess
    smoothed = lowess(df[scale], df["Calcium"], frac=0.3)
    ax.plot(smoothed[:, 0], smoothed[:, 1], color="#4C72B0", linewidth=2)

    # Spearman（使用 partial 结果但不写 partial）
    rho = df_sp.loc[scale, "Partial_Spearman_rho"]
    p   = df_sp.loc[scale, "p值"]
    ax.set_title(f"{titles[i]}\nSpearman ρ={rho:.3f}, p={p:.3f}", fontsize=13)

    ax.set_xlabel("Calcium")
    ax.set_ylabel(scale)

    # 横坐标刻度
    ax.set_xticks([100,200,300,400,500,600,700,800,900])
    ax.set_xlim(80, 920)

# ===========================
# ★ 在 subplot 外层左上角写 ABCD
# ===========================
# 使用 figure 坐标系（0–1）
panel_positions = [
    (0.02, 0.97),  # A 左上
    (0.52, 0.97),  # B 右上
    (0.02, 0.48),  # C 左下
    (0.52, 0.48),  # D 右下
]

for (x, y), lab in zip(panel_positions, labels):
    fig.text(x, y, lab, fontsize=18, fontweight="bold")

plt.tight_layout()
plt.savefig("Ca_scatter_4panel_ABCD.png", dpi=300)
plt.show()

print("已生成：Ca_scatter_4panel_ABCD.png")
