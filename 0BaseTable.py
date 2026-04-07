import pandas as pd
import numpy as np
import statsmodels.formula.api as smf

# =========================
# 0) 参数：你需要改的地方
# =========================
file_path = "data0116.xlsx"   # <- 改成你的Excel路径
sheet_name = 0                 # <- 默认第一个sheet；也可以写 "Sheet1"

# 你要控制的协变量（按你的表头写）
# 常见：年龄、性别、胃肠道症状、挑食、睡眠障碍、过敏、情绪...（看你论文设计）
# ⚠️ baseline表一般只控制年龄/性别；但你要求“需要做控制变量”，这里给你完全可配。
COVARS_RAW = ["年龄", "性别"]  # <- 想加就加，比如 ["年龄","性别","挑食","睡眠障碍"]

# 症状评分0-5：把它变成“有/无” 或 “中重度”
# 方案A：>0 视为有症状
# 方案B：>=3 视为中重度（0-5里更常用）
SYMPTOM_RULE = ">=3"  # <- 你可改成 ">0"

# 输出文件
out_path = "baseline_table_Ca_Q1_vs_Q4_adjusted.xlsx"


# =========================
# 1) 读取与清洗
# =========================
df = pd.read_excel(file_path, sheet_name=sheet_name)
df.columns = [str(c).strip() for c in df.columns]

# 需要数值化的列（按你给的字段）
num_cols = ["Calcium", "年龄", "ABC", "CARS", "ATEC-1", "ATEC-2", "ATEC-3", "ATEC-4", "ATEC-Total",
            "胃肠道症状", "挑食", "睡眠障碍", "兴奋", "过敏", "情绪"]
for c in num_cols:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

# 性别标准化
if "性别" in df.columns:
    df["性别"] = df["性别"].astype(str).str.strip()

# =========================
# 2) 四分位分组：Q1 vs Q4
# =========================
q1 = df["Calcium"].quantile(0.25)
q3 = df["Calcium"].quantile(0.75)

df["Ca_group"] = np.nan
df.loc[df["Calcium"] <= q1, "Ca_group"] = 0   # 低Ca
df.loc[df["Calcium"] >= q3, "Ca_group"] = 1   # 高Ca

sub = df[df["Ca_group"].isin([0, 1])].copy()
sub["Ca_group"] = sub["Ca_group"].astype(int)

n_low = int((sub["Ca_group"] == 0).sum())
n_high = int((sub["Ca_group"] == 1).sum())
print(f"Calcium Q1={q1:.3f}, Q3={q3:.3f} | low n={n_low}, high n={n_high}")

# =========================
# 3) 处理协变量：把性别做成0/1，避免公式里中文类别报错
# =========================
# 把 COVARS_RAW 映射为可用于模型的列名
sub = sub.copy()
covars_model = []

for c in COVARS_RAW:
    if c not in sub.columns:
        print(f"[WARN] 协变量列不存在：{c}（已跳过）")
        continue

    if c == "性别":
        # 你数据是 男/女；如果不是，改这里
        sub["sex_bin"] = (sub["性别"] == "男").astype(float)
        covars_model.append("sex_bin")
    else:
        covars_model.append(c)

# =========================
# 4) 工具函数：显示、二值化、检验
# =========================
def mean_sd(x):
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) == 0:
        return np.nan, np.nan
    return x.mean(), x.std(ddof=1)

def fmt_mean_sd(x):
    m, s = mean_sd(x)
    if pd.isna(m):
        return ""
    return f"{m:.2f}±{s:.2f}"

def star(p):
    try:
        return "*" if float(p) < 0.05 else ""
    except:
        return ""

def symptom_binary(series, rule=">0"):
    x = pd.to_numeric(series, errors="coerce")
    if rule == ">0":
        return (x > 0).astype(float)
    elif rule == ">=3":
        return (x >= 3).astype(float)
    elif rule == ">=2":
        return (x >= 2).astype(float)
    else:
        raise ValueError("SYMPTOM_RULE 只能用 '>0' 或 '>=2' 或 '>=3'")

def proportion_fmt(bin_series):
    s = pd.to_numeric(bin_series, errors="coerce")
    s = s.dropna()
    if len(s) == 0:
        return ""
    k = int(s.sum())
    n = int(len(s))
    pct = 100.0 * k / n if n > 0 else np.nan
    return f"{pct:.1f}% ({k}/{n})"

def ancova_pvalue(df_, y, covars):
    cols = [y, "Ca_group"] + covars
    d = df_[cols].dropna()
    if d.shape[0] < 5:
        return np.nan

    formula = f"{y} ~ Ca_group"
    for c in covars:
        formula += f" + {c}"

    try:
        model = smf.ols(formula, data=d).fit()
        return float(model.pvalues.get("Ca_group", np.nan))
    except Exception as e:
        print(f"[WARN] ANCOVA失败: {y} | {e}")
        return np.nan

def logistic_pvalue(df_, y_bin, covars):
    cols = [y_bin, "Ca_group"] + covars
    d = df_[cols].dropna()
    if d.shape[0] < 8:
        return np.nan

    # 如果全是0或全是1，logit会崩
    if d[y_bin].nunique() < 2:
        return np.nan

    formula = f"{y_bin} ~ Ca_group"
    for c in covars:
        formula += f" + {c}"

    try:
        model = smf.logit(formula, data=d).fit(disp=False)
        return float(model.pvalues.get("Ca_group", np.nan))
    except Exception as e:
        # 分离问题常见（complete separation）
        print(f"[WARN] Logistic失败: {y_bin} | {e}")
        return np.nan

# =========================
# 5) 你要输出的 baseline 行定义
# kind:
#   - continuous: mean±sd + ANCOVA(调整协变量)
#   - male_pct  : 男性比例 + Logistic(调整协变量)
#   - symptom_pct: 症状二值比例 + Logistic(调整协变量)
# =========================
rows = [
    ("年龄", "continuous", "年龄"),
    ("性别（男，%）", "male_pct", "性别"),
    ("胃肠道症状（%）", "symptom_pct", "胃肠道症状"),
    ("挑食（%）", "symptom_pct", "挑食"),
    ("兴奋（%）", "symptom_pct", "兴奋"),
    ("ABC", "continuous", "ABC"),
    ("CARS", "continuous", "CARS"),
    ("ATEC-Total", "continuous", "ATEC-Total"),
    ("Calcium", "continuous", "Calcium"),
]

# =========================
# 6) 生成表
# =========================
low = sub[sub["Ca_group"] == 0].copy()
high = sub[sub["Ca_group"] == 1].copy()

out = []
for label, kind, col in rows:
    if col not in sub.columns and kind != "male_pct":
        out.append([label, "", "", "", ""])
        continue

    if kind == "continuous":
        low_val = fmt_mean_sd(low[col])
        high_val = fmt_mean_sd(high[col])
        p = ancova_pvalue(sub, col, covars_model)
        out.append([label, low_val, high_val, p, star(p)])

    elif kind == "male_pct":
        # 显示：直接按原始男/女统计
        low_val = proportion_fmt((low["性别"] == "男").astype(float))
        high_val = proportion_fmt((high["性别"] == "男").astype(float))

        # 模型：用 sex_bin 做因变量（注意：此处不把 sex_bin 放入协变量，避免自己控制自己）
        if "sex_bin" not in sub.columns:
            sub["sex_bin"] = (sub["性别"] == "男").astype(float)

        covars_for_sex = [c for c in covars_model if c != "sex_bin"]
        p = logistic_pvalue(sub, "sex_bin", covars_for_sex)
        out.append([label, low_val, high_val, p, star(p)])

    elif kind == "symptom_pct":
        bin_col = f"{col}_bin"
        sub[bin_col] = symptom_binary(sub[col], rule=SYMPTOM_RULE)
        low_bin = symptom_binary(low[col], rule=SYMPTOM_RULE)
        high_bin = symptom_binary(high[col], rule=SYMPTOM_RULE)

        low_val = proportion_fmt(low_bin)
        high_val = proportion_fmt(high_bin)

        p = logistic_pvalue(sub, bin_col, covars_model)
        out.append([label, low_val, high_val, p, star(p)])

    else:
        raise ValueError(f"未知kind: {kind}")

table = pd.DataFrame(out, columns=[
    "变量",
    f"低Ca组（n={n_low}）",
    f"高Ca组（n={n_high}）",
    "P值（调整后）",
    "显著性"
])

# p值格式化
def fmt_p(x):
    try:
        x = float(x)
        return f"{x:.4f}"
    except:
        return ""
table["P值（调整后）"] = table["P值（调整后）"].apply(fmt_p)

print(table)

# =========================
# 7) 导出Excel
# =========================
with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
    table.to_excel(writer, index=False, sheet_name="Baseline_Q1Q4")

print(f"Saved: {out_path}")
