import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy.stats as stats
import scikit_posthocs as sp
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import mannwhitneyu, kruskal
from mpl_toolkits.mplot3d import Axes3D  

# -------------------- Global Plotting Style --------------------
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "axes.titlesize": 20,
    "legend.fontsize": 18
})

# -------------------- Figure S4 -------------------------------

df = pd.read_csv("base_targeted.csv")
df["GridIndex"] = df["GridIndex"].astype(int)
metrics = ["FractalDimension", "Eccentricity", "Lacunarity"]
metric_titles = {
    "FractalDimension": "Fractal Dimension",
    "Eccentricity": "Eccentricity",
    "Lacunarity": "Lacunarity"
}
label_map = {0: "No Switch", 1: "Switch"}
colors = {0: "#3366cc", 1: "#762a83"}
palette = {"No Switch": "#3366cc", "Switch": "#762a83"}

for metric in metrics:
    df[f"{metric}_resid"] = np.nan
    for grid_idx in [0, 1]:
        sub = df[df["GridIndex"] == grid_idx].copy()
        model = smf.ols(f"{metric} ~ TumorCellCount", data=sub).fit()
        df.loc[sub.index, f"{metric}_resid"] = model.resid

final_time = df["TimeStep"].max()
df_final = df[df["TimeStep"] == final_time].copy()
df_final["ConditionLabel"] = df_final["GridIndex"].map(label_map)

fig, axes = plt.subplots(3, 3, figsize=(20, 18))
axes = axes.flatten()

for i, metric in enumerate(metrics):
    ax = axes[i]
    summary = df.groupby(["TimeStep", "GridIndex"])[metric].agg(["mean", "std"]).reset_index()
    for grid_idx in [0, 1]:
        group = summary[summary["GridIndex"] == grid_idx]
        ax.plot(group["TimeStep"], group["mean"], '-o', color=colors[grid_idx], label=label_map[grid_idx])
        ax.fill_between(group["TimeStep"], group["mean"] - group["std"], group["mean"] + group["std"], color=colors[grid_idx], alpha=0.3)
    ax.set_title(f"{metric_titles[metric]}: Mean ± STD")
    ax.set_xlabel("Time Step")
    ax.set_ylabel(metric_titles[metric])
    ax.legend()

for i, metric in enumerate(metrics):
    ax = axes[i + 3]
    resid_col = f"{metric}_resid"
    summary = df.groupby(["TimeStep", "GridIndex"])[resid_col].agg(["mean", "std"]).reset_index()
    for grid_idx in [0, 1]:
        group = summary[summary["GridIndex"] == grid_idx]
        ax.plot(group["TimeStep"], group["mean"], '-o', color=colors[grid_idx], label=label_map[grid_idx])
        ax.fill_between(group["TimeStep"], group["mean"] - group["std"], group["mean"] + group["std"], color=colors[grid_idx], alpha=0.3)
    ax.set_title(f"{metric_titles[metric]}: Residuals Over Time")
    ax.set_xlabel("Time Step")
    ax.set_ylabel("Residual")
    ax.legend()

for i, metric in enumerate(metrics):
    ax = axes[i + 6]
    resid_col = f"{metric}_resid"
    sns.boxplot(x="ConditionLabel", y=resid_col, data=df_final, palette=palette, showfliers=False, ax=ax)
    sns.stripplot(x="ConditionLabel", y=resid_col, data=df_final, color='black', alpha=0.5, ax=ax)
    ax.set_title(f"{metric_titles[metric]}: Residuals at Final Time Step")
    ax.set_xlabel("")
    ax.set_ylabel("Residual")

plt.tight_layout()
plt.savefig("S2_Morphology_TemporalResiduals.png", dpi=600, bbox_inches="tight")
plt.show()

# -------------------- Figure S3 -------------------------------

df = pd.read_csv("base_psa.csv", on_bad_lines='skip')
print("Successfully loaded dataframe shape:", df.shape)

df = df.reset_index()
df = df.rename(columns={"index": "ParamSetIndex"})

desired_columns = [
    "ParamSetIndex", "Replication", "TimeStep", "GridIndex", 
    "effectAntiMet", "effectProMet", "conversionThreshold", 
    "S2", "S4", "effectPerTumorCell", "TumorCellCount", 
    "FractalDimension", "Lacunarity", "Eccentricity"
]
df = df.iloc[:, :len(desired_columns)]
df.columns = desired_columns

df = df.rename(columns={
    "S2": "switchSensitivity",
    "S4": "divisionSensitivity"
})

def plot_3d_scatter(df, x, y, z, title, filename):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")

    ax.scatter(df[x], df[y], df[z], alpha=0.5)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_zlabel(z)
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.show()

plot_3d_scatter(
    df,
    x="switchSensitivity",
    y="divisionSensitivity",
    z="effectPerTumorCell",
    title="3D Projection of Sensitivity vs. Tumor Effect",
    filename="S3_param_sampling_sensitivity_effect.png"
)

plot_3d_scatter(
    df,
    x="effectAntiMet",
    y="effectProMet",
    z="conversionThreshold",
    title="3D Projection of Metastasis Effects",
    filename="S3_param_sampling_effects.png"
)

# -------------------- Figure S5 -------------------------------

# Filter to the final time step only
final_time = df["TimeStep"].max()
df_final = df[df["TimeStep"] == final_time].copy()

# Aggregate by parameter set
df_avg = df_final.groupby("ParamSetIndex", as_index=False)[param_list + output_list].mean()

param_list = [
    "effectAntiMet", "effectProMet", "conversionThreshold",
    "switchSensitivity", "divisionSensitivity", "effectPerTumorCell"
]
output_list = [
    "TumorCellCount", "FractalDimension", "Lacunarity", "Eccentricity"
]
discrete_threshold = 10

fig, axes = plt.subplots(nrows=len(output_list), ncols=len(param_list), figsize=(18, 12), sharey=False)

for row_idx, output_name in enumerate(output_list):
    for col_idx, param_name in enumerate(param_list):
        ax = axes[row_idx, col_idx]
        x = df_avg[param_name].values
        y = df_avg[output_name].values
        unique_vals = np.unique(x)

        if len(unique_vals) < discrete_threshold:
            temp_df = pd.DataFrame({param_name: x, output_name: y})
            order = np.sort(unique_vals)
            sns.boxplot(x=param_name, y=output_name, data=temp_df, order=order, ax=ax, width=0.6)
            sns.stripplot(x=param_name, y=output_name, data=temp_df, order=order, ax=ax, color='black', alpha=0.5)
        else:
            ax.scatter(x, y, alpha=0.5, s=10, label="Data")
            lowess_results = sm.nonparametric.lowess(y, x, frac=0.3)
            ax.plot(lowess_results[:, 0], lowess_results[:, 1], color='purple', linewidth=2, label="LOWESS")

        if row_idx == len(output_list) - 1:
            ax.set_xlabel(param_name, fontsize=12)
        if col_idx == 0:
            ax.set_ylabel(output_name, fontsize=12)
        if row_idx == 0 and col_idx == 0:
            ax.legend(fontsize="small")

plt.tight_layout()
plt.savefig("SUPP_4_ParamVsOutputs_Grid.png", dpi=600, bbox_inches="tight")
plt.show()

# -------------------- Figure S7 -------------------------------

df0 = pd.read_csv("tumorMod_0.0.csv")
df05 = pd.read_csv("tumorMod_0.5.csv")
df1 = pd.read_csv("tumorMod_1.0.csv")

df0['resistance'] = 1.0 - df0['modFactor']
df05['resistance'] = 1.0 - df05['modFactor']
df1['resistance'] = 1.0 - df1['modFactor']

all_res = pd.concat([df0['resistance'], df05['resistance'], df1['resistance']])
vmin = all_res.min()
vmax = all_res.max()

fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
titles = [r"$G_f$ = 0", r"$G_f$ = 0.5", r"$G_f$ = 1.0"]
dfs = [df0, df05, df1]

for ax, df, title in zip(axes, dfs, titles):
    sc = ax.scatter(df['x'], df['y'], c=df['resistance'], cmap='viridis',
                    vmin=vmin, vmax=vmax, s=1, alpha=0.9)
    ax.set_title(title, fontsize=16)
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), orientation="vertical", fraction=0.025, pad=0.04)
cbar.set_label(r"Effective Drug Exposure: $1 - G_f \cdot \frac{N_{\mathrm{pro}}}{8}$", fontsize=12, labelpad=10)

plt.suptitle("S7: Tumor Cells Colored by Effective Drug Exposure Across $G_f$ Values", fontsize=16, y=0.98)
plt.tight_layout(rect=[1.7, 0, 0.9, 1.4])
plt.savefig("S7_ChemoResistance_Plot.png", dpi=600, bbox_inches="tight")
plt.show()
