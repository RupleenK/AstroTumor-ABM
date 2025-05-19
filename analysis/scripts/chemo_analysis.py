import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, shapiro, levene, kruskal
import statsmodels.formula.api as smf
import statsmodels.api as sm
import scikit_posthocs as sp
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from sklearn.cluster import KMeans
import pingouin as pg

# ---------------- Global Style Settings ----------------
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "axes.titlesize": 20,
    "legend.fontsize": 18
})

# ---------------- Load Data and Filter to Final Time Step ----------------
df_chemo = pd.read_csv("chemo_targeted.csv")
last_time = df_chemo["TimeStep"].max()
df_final = df_chemo[df_chemo["TimeStep"] == last_time].copy()

# ---------------- Group by Condition ----------------
ctrl = df_final[df_final["GridIndex"] == 0]["TumorCellCount"].dropna()
switch = df_final[df_final["GridIndex"] == 1]["TumorCellCount"].dropna()

# ---------------- Statistics ----------------
mean_ctrl, std_ctrl = ctrl.mean(), ctrl.std()
mean_switch, std_switch = switch.mean(), switch.std()

print("=== Final Tumor Burden Comparison ===")
print(f"No Switch: mean = {mean_ctrl:.0f}, std = {std_ctrl:.0f}")
print(f"Switch:    mean = {mean_switch:.0f}, std = {std_switch:.0f}")

# ---------------- Mann–Whitney U Test ----------------
u_stat, p_val = mannwhitneyu(ctrl, switch, alternative='two-sided')
print(f"Mann–Whitney U test p = {p_val:.3g}")

# ================= Figure 8: Chemo Response Time Series Analysis =================
df_chemo = pd.read_csv("chemo_targeted.csv")
# Group data by condition and time step, then compute mean, standard deviation and count.
grouped = df_chemo.groupby(["GridIndex", "TimeStep"])["TumorCellCount"].agg(["mean", "std", "count"]).reset_index()
# Compute standard error of the mean
grouped["sem"] = grouped["std"] / np.sqrt(grouped["count"])

fig, ax = plt.subplots(figsize=(8,6))
colors = {0: "#377eb8", 1: "#984ea3"}
for grid in [0, 1]:
    sub = grouped[grouped["GridIndex"] == grid].sort_values("TimeStep")
    label = "No Switch" if grid == 0 else "Switch"
    ax.plot(sub["TimeStep"], sub["mean"], color=colors[grid], label=label, linewidth=3)
    ax.fill_between(sub["TimeStep"],
                    sub["mean"]-sub["sem"],
                    sub["mean"]+sub["sem"],
                    color=colors[grid], alpha=0.3)
ax.set_xlabel("Time Step", fontsize=20)
ax.set_ylabel("Tumor Cell Count\n(Mean ± SEM)", fontsize=20)
ax.set_title("Chemotherapy Response", fontsize=20)
ax.legend(fontsize=20)
plt.tight_layout()
plt.savefig("Chemo_Targeted.png", dpi=300, bbox_inches="tight")
plt.show()
