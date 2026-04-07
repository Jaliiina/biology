# -*- coding: utf-8 -*-
"""
全样本稳健关联（多 severity，不做 FDR）
severity ∈ {ATEC-1, ATEC-2, ATEC-3, ATEC-4, ATEC-Total, ABC, CARS}

模型：
  severity ~ microbe (CLR) + Calcium + 年龄 + 性别
显著性标准：p < 0.05
"""

import re
import numpy as np
import pandas as pd
import statsmodels.api as sm

# =============================
# 0) 文件路径
# =============================
TAX_FILE   = "taxonomy.tsv"
LEFSE_FILE = "lefse_species_high_vs_low.csv"
PHENO_FILE = "data0116.xlsx"

OUT_ALL = "regression_species_multi_severity_all_p005.csv"
OUT_SIG = "regression_species_multi_severity_sig_p005.csv"

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
USE_ROBUST_SE = True
FILTER_PLACEHOLDER = True   # 去掉 GGB/SGB

# =============================
# 2) 工具函数
# =============================
def extract_species_tag(x):
    m = re.search(r"(s__[^|]+)", str(x))
    return m.group(1) if m else None

def clean_species(tag):
    return tag.replace("s__", "") if tag else None

def is_placeholder(name):
    if name is None:
        return True
    if re.match(r"^(GGB|SGB)", name, flags=re.IGNORECASE):
        return True
    if re.search(r"_SGB\d+", name):
        return True
    return False

def clr(df, pc=1e-6):
    x = df.astype(float) + pc
    logx = np.log(x)
    return logx.sub(logx.mean(axis=1), axis=0)

# =============================
# 3) 读取 phenotype
# =============================
pheno = pd.read_excel(PHENO_FILE)
pheno["sample_id"] = pheno["F号"].astype(str).str.strip()

need = ["sample_id", CALCIUM_COL, AGE_COL, SEX_COL] + SEVERITY_COLS
ph = pheno[need].copy()

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
# 4) LEfSe 候选菌
# =============================
lefse = pd.read_csv(LEFSE_FILE)
lefse["species_tag"] = lefse["Taxa"].astype(str).apply(extract_species_tag)
lefse = lefse.dropna(subset=["species_tag"])

candidate_tags = lefse["species_tag"].unique().tolist()

# =============================
# 5) taxonomy → species abundance
# =============================
tax = pd.read_csv(TAX_FILE, sep="\t", dtype=str)
tax = tax.rename(columns={tax.columns[0]: "ID"})
tax = tax[~tax["ID"].str.contains("UNMAPPED|UNINTEGRATED", case=False, na=False)]

tax["species_tag"] = tax["ID"].apply(extract_species_tag)
tax = tax.dropna(subset=["species_tag"])

sample_cols = [c for c in tax.columns if c not in ["ID", "species_tag"]]
tax[sample_cols] = tax[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

tax_sp = tax.groupby("species_tag")[sample_cols].sum()
abund = tax_sp.T
abund.index = abund.index.astype(str).str.strip()

# 样本对齐
common = sorted(set(abund.index) & set(ph.index))
abund = abund.loc[common]
ph = ph.loc[common]

# 只保留候选菌
candidate_tags = [t for t in candidate_tags if t in abund.columns]
abund = abund[candidate_tags]

# 相对丰度 + CLR
abund_rel = abund.div(abund.sum(axis=1), axis=0).fillna(0)
abund_clr = clr(abund_rel, PSEUDOCOUNT)

# =============================
# 6) 回归（多 severity）
# =============================
rows = []

X_base = pd.DataFrame({
    "Calcium": ph[CALCIUM_COL],
    "Age": ph[AGE_COL],
    "Sex": ph["sex_bin"]
}, index=ph.index)

for sev in SEVERITY_COLS:
    y = ph[sev].dropna()

    for tag in candidate_tags:
        name = clean_species(tag)
        if FILTER_PLACEHOLDER and is_placeholder(name):
            continue

        x_micro = abund_clr[tag]
        dat = pd.concat([y, X_base, x_micro], axis=1).dropna()
        if dat.shape[0] < 10:
            continue

        y_m = dat[sev]
        X_m = dat[["Calcium", "Age", "Sex", tag]].rename(columns={tag: "Microbe"})
        X_m = sm.add_constant(X_m)

        model = sm.OLS(y_m, X_m)
        fit = model.fit(cov_type="HC3") if USE_ROBUST_SE else model.fit()

        beta = fit.params["Microbe"]
        se   = fit.bse["Microbe"]
        p    = fit.pvalues["Microbe"]

        rows.append({
            "severity": sev,
            "species": name,
            "n": dat.shape[0],
            "beta": beta,
            "ci95_low": beta - 1.96 * se,
            "ci95_up":  beta + 1.96 * se,
            "p": p
        })

res = pd.DataFrame(rows).sort_values(["severity", "p"])
res.to_csv(OUT_ALL, index=False, encoding="utf-8-sig")

sig = res[res["p"] < 0.05]
sig.to_csv(OUT_SIG, index=False, encoding="utf-8-sig")

print("全部结果：", OUT_ALL)
print("p < 0.05 显著结果：", OUT_SIG)
