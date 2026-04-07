import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import statsmodels.api as sm
import matplotlib.pyplot as plt

# ======================
# 1) 读入数据
# ======================
df = pd.read_excel("data0116.xlsx")

plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False

# ======================
# 2) 变量
# ======================
target = "Calcium"
age_col = "年龄"
score_cols = ["ATEC-Total","ABC","CARS","ATEC-1","ATEC-2","ATEC-3","ATEC-4"]

# 转数值
df[target] = pd.to_numeric(df[target], errors="coerce")
df[age_col] = pd.to_numeric(df[age_col], errors="coerce")
for c in score_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# 回归残差函数
def get_residual(y, x):
    x = sm.add_constant(x)
    model = sm.OLS(y, x, missing="drop").fit()
    return model.resid

# ======================
# 3) 计算“控制年龄后的偏 Spearman”
# ======================
records = []

for col in score_cols:
    sub = df[[target, age_col, col]].dropna()
    
    # 计算残差（去掉年龄贡献）：
    y_res = get_residual(sub[col], sub[age_col])      # 量表残差
    ca_res = get_residual(sub[target], sub[age_col])  # 钙残差
    
    # 偏 Spearman
    rho, p = spearmanr(ca_res, y_res)

    # 方向解释
    if np.isnan(rho):
        direction = ""
    elif rho > 0:
        direction = "校正年龄后：Ca 越高，该量表越高（同向）"
    elif rho < 0:
        direction = "校正年龄后：Ca 越高，该量表越低（反向）"
    else:
        direction = "无相关"

    records.append([col, len(sub), rho, p, direction])

result_df = pd.DataFrame(
    records,
    columns=["量表", "有效样本数", "Partial_Spearman_rho", "p值", "方向解释"]
)

print("=== 控制年龄后的 Partial Spearman 相关 ===")
print(result_df)

result_df.to_csv("Ca_vs_scales_partialSpearman_ageAdjusted.csv", index=False)
