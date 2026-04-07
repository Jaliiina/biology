# -*- coding: utf-8 -*-
"""
全样本：通路 ~ 严重程度（不做 FDR，只看 p < 0.05）
候选通路来自 LEfSe（High vs Low Ca）

模型：
  severity ~ pathway + Calcium + 年龄 + 性别

输出：
  1) regression_pathway_multi_severity_all.csv
  2) regression_pathway_multi_severity_sig_p005.csv
"""

import re
import numpy as np
import pandas as pd
import statsmodels.api as sm

# =============================
# 0) 文件路径
# =============================
LEFSE_FILE = "lefse_pathway_high_vs_low.csv"
PATH_FILE  = "pathabundance_relab_stratified.tsv"
PHENO_FILE = "data0116.xlsx"

OUT_ALL = "regression_pathway_multi_severity_all.csv"
OUT_SIG = "regression_pathway_multi_severity_sig_p005.csv"

# =============================
# 1) 参数
# =============================
SEVERITY_COLS = [
    "ATEC-1", "ATEC-2", "ATEC-3", "ATEC-4",
    "ATEC-Total", "ABC", "CARS"
]

CALCIUM_COL = "Calcium"
AGE_COL     = "年龄"
SEX_COL     = "性别"

PSEUDOCOUNT = 1e-6
USE_ROBUST_SE = True   # HC3 稳健标准误

# =============================
# 2) 小工具函数
# =============================
def clr(df, pc=1e-6):
    x = df.astype(float) + pc
    logx = np.log(x)
    return logx.sub(logx.mean(axis=1), axis=0)

def base_pathway_name(s):
    return str(s).split("|")[0].strip()

def extract_Fid(x):
    m = re.search(r"(F\d+)", str(x))
    return m.group(1) if m else None

# =============================
# 3) 读取 phenotype（全样本）
# =============================
pheno = pd.read_excel(PHENO_FILE)
pheno["sample_id"] = pheno["F号"].astype(str).str.strip()

need_cols = ["sample_id", CALCIUM_COL, AGE_COL, SEX_COL] + SEVERITY_COLS
ph = pheno[need_cols].copy()

for c in [CALCIUM_COL, AGE_COL] + SEVERITY_COLS:
    ph[c] = pd.to_numeric(ph[c], errors="coerce")

sex_map = {
    "女": 1, "女性": 1, "F": 1, "female": 1,
    "男": 0, "男性": 0, "M": 0, "male": 0
}
ph["sex_bin"] = ph[SEX_COL].astype(str).str.strip().map(sex_map)

ph = ph.dropna(subset=[CALCIUM_COL, AGE_COL, "sex_bin"])
ph = ph.set_index("sample_id")

# =============================
# 4) 读取 LEfSe 候选通路
# =============================
lefse = pd.read_csv(LEFSE_FILE)

if "Taxa" not in lefse.columns:
    raise ValueError("LEfSe 文件中找不到列 'Taxa'")

lefse["pathway"] = lefse["Taxa"].map(base_pathway_name)
candidate_paths = sorted(lefse["pathway"].dropna().unique().tolist())

if len(candidate_paths) == 0:
    raise ValueError("LEfSe 候选通路为空")

# =============================
# 5) 读取通路丰度表并汇总到通路层
# =============================
path_raw = pd.read_csv(PATH_FILE, sep="\t", dtype=str)
first_col = path_raw.columns[0]
path_raw = path_raw.rename(columns={first_col: "Pathway"})

mask_rm = path_raw["Pathway"].str.contains("UNMAPPED|UNINTEGRATED", case=False, na=False)
path_raw = path_raw.loc[~mask_rm].copy()

sample_cols = [c for c in path_raw.columns if c != "Pathway"]
col_map = {c: extract_Fid(c) for c in sample_cols}
keep_cols = [c for c in sample_cols if col_map[c] is not None]

path_raw = path_raw[["Pathway"] + keep_cols]
path_raw = path_raw.rename(columns={c: col_map[c] for c in keep_cols})

for c in path_raw.columns:
    if c != "Pathway":
        path_raw[c] = pd.to_numeric(path_raw[c], errors="coerce").fillna(0)

path_raw["pathway_base"] = path_raw["Pathway"].map(base_pathway_name)
path_sum = path_raw.groupby("pathway_base")[path_raw.columns[1:-1]].sum()

path_abund = path_sum.T
path_abund.index = path_abund.index.astype(str).str.strip()

# 样本对齐
common = sorted(set(path_abund.index) & set(ph.index))
path_abund = path_abund.loc[common]
ph = ph.loc[common]

# 只保留 LEfSe 候选通路
candidate_paths = [p for p in candidate_paths if p in path_abund.columns]
path_abund = path_abund[candidate_paths]

# 相对丰度 + CLR
path_rel = path_abund.div(path_abund.sum(axis=1), axis=0).fillna(0)
path_clr = clr(path_rel, PSEUDOCOUNT)

# =============================
# 6) 回归分析（无 FDR）
# =============================
X_base = pd.DataFrame({
    "Calcium": ph[CALCIUM_COL],
    "Age": ph[AGE_COL],
    "Sex": ph["sex_bin"]
}, index=ph.index)

rows = []

for sev in SEVERITY_COLS:
    y = ph[sev].dropna()

    for pw in candidate_paths:
        dat = pd.concat(
            [y, X_base, path_clr[pw].rename("Pathway")],
            axis=1
        ).dropna()

        if dat.shape[0] < 10:
            continue

        y_m = dat[sev]
        X_m = sm.add_constant(dat[["Calcium", "Age", "Sex", "Pathway"]])

        model = sm.OLS(y_m, X_m)
        fit = model.fit(cov_type="HC3") if USE_ROBUST_SE else model.fit()

        beta = fit.params["Pathway"]
        se   = fit.bse["Pathway"]
        pval = fit.pvalues["Pathway"]

        rows.append({
            "severity": sev,
            "pathway": pw,
            "n": dat.shape[0],
            "beta": beta,
            "ci95_low": beta - 1.96 * se,
            "ci95_up":  beta + 1.96 * se,
            "p": pval
        })

res = pd.DataFrame(rows).sort_values(["severity", "p"])

res.to_csv(OUT_ALL, index=False, encoding="utf-8-sig")
sig = res[res["p"] < 0.05]
sig.to_csv(OUT_SIG, index=False, encoding="utf-8-sig")

print("全部结果：", OUT_ALL)
print("p < 0.05 显著结果：", OUT_SIG)
print("显著条目数：", sig.shape[0])
