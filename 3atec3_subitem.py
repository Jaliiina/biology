import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import statsmodels.api as sm

# ========== 1) Load data ==========
df = pd.read_excel("data0116.xlsx")

# Target & age
target = "Calcium"
age_col = "年龄"

# English translation of behavioral items (ATEC-3 subitems)
items_en = {
    "知道别人叫自己的名字": "Responds to own name",
    "会对称赞做出回应": "Responds to praise",
    "会看人和动物": "Looks at people/animals",
    "会看图片（和电视节目）": "Looks at pictures/TV",
    "会画画和涂鸦": "Draws or scribbles",
    "会正确使用玩具": "Uses toys correctly",
    "面部表情正常": "Normal facial expressions",
    "能理解电视剧情或者卡通片剧情": "Understands TV/cartoon plots",
    "能够听懂说明，例如，如何打开一个盒子": "Understands instructions",
    "能够发现周围环境的改变，例如妈妈换发型了": "Notices environmental changes",
    "危险感知": "Sense of danger",
    "具备想象力": "Imagination",
    "积极主动": "Initiative",
    "自己穿衣服": "Dresses self",
    "好奇、感兴趣": "Curiosity/interest",
    "敢于冒险和探索": "Adventurous/exploratory",
    "会自我调整，不浮躁": "Self-regulation",
    "会去看别人正在看的东西": "Joint attention"
}

items_cn = list(items_en.keys())     # 中文列名
items = list(items_en.values())      # 英文列名

# ========== 2) 数值化 ==========
df[target] = pd.to_numeric(df[target], errors="coerce")
df[age_col] = pd.to_numeric(df[age_col], errors="coerce")
for c in items_cn:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# 回归取残差函数（y ~ age）
def get_residual(y, age):
    X = sm.add_constant(age)
    model = sm.OLS(y, X, missing="drop").fit()
    return model.resid

# ========== 3) 计算“校正年龄后的 Spearman” ==========
records = []

for cn, en in items_en.items():
    sub = df[[target, age_col, cn]].dropna()
    if len(sub) < 5:
        rho, p = np.nan, np.nan
    else:
        # 去掉年龄影响
        y_res  = get_residual(sub[cn],      sub[age_col])
        ca_res = get_residual(sub[target],  sub[age_col])
        rho, p = spearmanr(ca_res, y_res)

    records.append((en, rho, p))

corr_df = pd.DataFrame(records, columns=["Variable", "rho", "p"])

# ========== 4) Significance stars & severity ==========
def p_to_stars(p):
    if pd.isna(p):
        return ""
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

corr_df["Significance"] = corr_df["p"].apply(p_to_stars)

def p_to_severity(p):
    if pd.isna(p) or p <= 0:
        return np.nan
    return -np.log10(p)

corr_df["Severity"] = corr_df["p"].apply(p_to_severity)

corr_df_sorted = corr_df.sort_values("Severity", ascending=False)
corr_df_sorted.to_csv("Ca_behavior_partialSpearman_ageAdjusted_results.csv", index=False)
print("Saved: Ca_behavior_partialSpearman_ageAdjusted_results.csv")

# ========== 5) Heatmap (English, only significant items) ==========
sig_df = corr_df_sorted[corr_df_sorted["p"] < 0.05].copy()

if sig_df.empty:
    print("No significant items (p < 0.05). No heatmap generated.")
else:
    plt.figure(figsize=(9, len(sig_df) * 0.45))

    pivot_df = sig_df.set_index("Variable")[["Severity"]]

    # Annotation text: rho + stars
    annot_text = sig_df.apply(
        lambda row: f"{row['rho']:.2f}\n{row['Significance']}", axis=1
    ).values.reshape(-1, 1)

    sns.heatmap(
        pivot_df,
        annot=annot_text,
        cmap="Reds",
        vmin=0,
        vmax=np.nanmax(pivot_df.values),
        linewidths=0.5,
        cbar_kws={"label": "-log10(p) (Significance severity)"},
        fmt="",
        annot_kws={"fontsize": 10, "fontweight": "bold", "color": "black"}
    )

    plt.title("Significant associations between calcium levels and ATEC-3 subitems",
              fontsize=14, pad=10)
    plt.xlabel("")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig("Ca_behavior_severity_heatmap_sig_only_ageAdjusted.png", dpi=300)
    plt.show()

    print("Saved: Ca_behavior_severity_heatmap_sig_only_ageAdjusted.png")
