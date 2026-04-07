import pandas as pd
import numpy as np
import plotly.graph_objects as go

# ========= 1) 读入数据 =========
# 你的 Excel 里至少要有这几列：
# severity, species, pathway, rho_sp_pw, beta_pw_sev
excel_path = "species_pathway_severity_merged.csv"   # 改成你的路径
               # 或者 "Sheet1"

df = pd.read_csv(excel_path)

needed = ["severity", "species", "pathway", "rho_sp_pw", "beta_pw_sev"]
missing = [c for c in needed if c not in df.columns]
if missing:
    raise ValueError(f"缺少列: {missing}")

df = df.dropna(subset=["severity", "species", "pathway"]).copy()

# 如果 severity 是数字（比如 0/1/2/3），也可以先转成更好看的标签
# df["severity"] = df["severity"].map({0:"mild",1:"moderate",2:"severe"}).fillna(df["severity"])

# ========= 2) 颜色映射：红正蓝负，强度按 |v| =========
def signed_red_blue_rgba(v, vmax, min_alpha=0.15, max_alpha=0.90):
    """
    v: signed value
    vmax: max(abs(v)) 用来归一化
    return rgba string
    """
    if pd.isna(v):
        return f"rgba(180,180,180,{min_alpha})"
    vmax = float(vmax) if vmax and vmax > 0 else 1.0
    a = min_alpha + (max_alpha - min_alpha) * min(abs(v) / vmax, 1.0)

    if v >= 0:
        # 红
        return f"rgba(220,20,60,{a})"
    else:
        # 蓝
        return f"rgba(30,144,255,{a})"

# ========= 3) 聚合成两段边：species->pathway，pathway->severity =========
# 物种->通路：同一对 (species, pathway) 可能出现多行，聚合一下
sp_pw = (
    df.groupby(["species", "pathway"], as_index=False)
      .agg(rho_sp_pw=("rho_sp_pw", "mean"))
)
# 通路->severity：同一对 (pathway, severity) 聚合
pw_sev = (
    df.groupby(["pathway", "severity"], as_index=False)
      .agg(beta_pw_sev=("beta_pw_sev", "mean"))
)

# 条带粗细：默认用 abs(系数)，你可以替换成别的权重列
sp_pw["value"] = sp_pw["rho_sp_pw"].abs()
pw_sev["value"] = pw_sev["beta_pw_sev"].abs()

# 可选：避免全是 0 导致图画不出来
eps = 1e-6
sp_pw.loc[sp_pw["value"] <= 0, "value"] = eps
pw_sev.loc[pw_sev["value"] <= 0, "value"] = eps

# ========= 4) 节点（node）编码 =========
species_list = sp_pw["species"].astype(str).unique().tolist()
pathway_list  = pd.Index(sp_pw["pathway"].astype(str).unique()).union(
               pw_sev["pathway"].astype(str).unique()).tolist()
severity_list = pw_sev["severity"].astype(str).unique().tolist()

# 为了避免不同层同名冲突，给每层加前缀做“内部ID”，展示时再去前缀
sp_nodes  = [f"sp::{x}" for x in species_list]
pw_nodes  = [f"pw::{x}" for x in pathway_list]
sev_nodes = [f"sev::{x}" for x in severity_list]

all_nodes = sp_nodes + pw_nodes + sev_nodes

node_index = {name: i for i, name in enumerate(all_nodes)}

def disp_label(node_id: str) -> str:
    # 去掉前缀用于展示
    return node_id.split("::", 1)[1]

labels = [disp_label(n) for n in all_nodes]

# ========= 5) 构造 links =========
# 颜色归一化尺度（分别按各自最大绝对值）
rho_max  = float(np.nanmax(np.abs(sp_pw["rho_sp_pw"].values))) if len(sp_pw) else 1.0
beta_max = float(np.nanmax(np.abs(pw_sev["beta_pw_sev"].values))) if len(pw_sev) else 1.0

sources = []
targets = []
values  = []
colors  = []

# species -> pathway (rho_sp_pw 着色)
for _, r in sp_pw.iterrows():
    s = f"sp::{str(r['species'])}"
    t = f"pw::{str(r['pathway'])}"
    sources.append(node_index[s])
    targets.append(node_index[t])
    values.append(float(r["value"]))
    colors.append(signed_red_blue_rgba(r["rho_sp_pw"], rho_max))

# pathway -> severity (beta_pw_sev 着色)
for _, r in pw_sev.iterrows():
    s = f"pw::{str(r['pathway'])}"
    t = f"sev::{str(r['severity'])}"
    sources.append(node_index[s])
    targets.append(node_index[t])
    values.append(float(r["value"]))
    colors.append(signed_red_blue_rgba(r["beta_pw_sev"], beta_max))

# 节点颜色：给个中性灰（你也可以按层分色）
node_color = ["rgba(160,160,160,0.75)"] * len(all_nodes)

# ========= 6) 画图 =========
fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=15,
        thickness=18,
        line=dict(color="rgba(80,80,80,0.35)", width=0.5),
        label=labels,
        color=node_color,
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=colors,
    )
)])

fig.update_layout(
    title="Species → Pathway → Severity Sankey (rho / beta colored)",
    font=dict(size=12),
    width=1200,
    height=700,
)

fig.show()

# 如果你想导出 HTML：
# fig.write_html("sankey_species_pathway_severity.html")
