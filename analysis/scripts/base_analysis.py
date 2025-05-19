import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu 
import pingouin as pg

# ---------------- Global Style Settings ---------------
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "axes.titlesize": 20,
    "legend.fontsize": 18,
})

# ---------------- Cohen's d Function ----------------
def cohens_d(x, y):
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx - 1)*np.var(x, ddof=1) + (ny - 1)*np.var(y, ddof=1)) / dof)

# ---------------- Data Loading and Preprocessing ----------------
df_base = pd.read_csv("base_targeted.csv")
last_time = df_base["TimeStep"].max()

# Filter data to only include the final time step
df_base_final = df_base[df_base["TimeStep"] == last_time].copy()

# Create separate dataframes for control (GridIndex=0) and reprogramming condition (GridIndex=1)
control = df_base_final[df_base_final["GridIndex"] == 0]
reprog  = df_base_final[df_base_final["GridIndex"] == 1]

# ---------------- Tumor Cell Count Metrics Analysis ----------------
mean_control = control["TumorCellCount"].mean()
std_control  = control["TumorCellCount"].std()
mean_reprog  = reprog["TumorCellCount"].mean()
std_reprog   = reprog["TumorCellCount"].std()

# Percent increase in tumor cell count from control to reprogramming
percent_increase = ((mean_reprog - mean_control) / mean_control) * 100

# Statistical tests: Mann–Whitney U and Cohen's d
u_stat, p_val = mannwhitneyu(control["TumorCellCount"], reprog["TumorCellCount"], alternative='two-sided')
d_value = cohens_d(control["TumorCellCount"], reprog["TumorCellCount"])

# Print analysis results
print("===== Tumor Cell Count Analysis =====")
print(f"Control: {mean_control:.2f} ± {std_control:.2f}")
print(f"Reprog: {mean_reprog:.2f} ± {std_reprog:.2f}")
print(f"Percent Increase: {percent_increase:.2f}%")
print(f"Mann-Whitney U: {u_stat}, p-value: {p_val:.4g}")
print(f"Cohen's d: {d_value:.3f}")
print("\n")

# ---------------- Morphology Metrics Analysis ----------------
metrics = ["FractalDimension", "Eccentricity", "Lacunarity"]
results = []  

# Loop over each metric
for metric in metrics:
    # Drop missing values in each group for this metric
    ctrl_vals = control[metric].dropna()
    reprog_vals = reprog[metric].dropna()
    
    # Calculate means and standard deviations for control and reprogramming groups
    mean_ctrl = ctrl_vals.mean()
    std_ctrl  = ctrl_vals.std()
    mean_reprog = reprog_vals.mean()
    std_reprog  = reprog_vals.std()
    
    # Run the Mann-Whitney U test
    u_stat, p_val = mannwhitneyu(ctrl_vals, reprog_vals, alternative='two-sided')
    
    # Append the computed statistics to the results list
    results.append({
        "Metric": metric,
        "No Switch Mean ± Std": f"{mean_ctrl:.3f} ± {std_ctrl:.3f}",
        "Switch Mean ± Std":  f"{mean_reprog:.3f} ± {std_reprog:.3f}",
        "U_stat": u_stat,
        "p_value": p_val
    })

# Print the morphology metrics analysis results
print("===== Morphology Metrics Analysis (No Switch vs. Switch) =====")
for res in results:
    print(f"Metric: {res['Metric']}")
    print(f"  No Switch: {res['No Switch - Mean ± Std']}")
    print(f"  Switch:  {res['Switch - Mean ± Std']}")
    print(f"  Mann-Whitney U: {res['U_stat']}, p-value: {res['p_value']:.3g}\n")

# --------------- Figure 4B: Final Tumor Cell Count Over Time ---------------
time_summary = df_base.groupby(["TimeStep", "GridIndex"])["TumorCellCount"].agg(["mean", "std", "count"]).reset_index()
time_summary["sem"] = time_summary["std"] / np.sqrt(time_summary["count"])
time_control = time_summary[time_summary["GridIndex"] == 0]
time_reprog  = time_summary[time_summary["GridIndex"] == 1]

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(time_control["TimeStep"], time_control["mean"], '-o', label="No Switch", color="#3366cc", zorder=2)
ax.fill_between(time_control["TimeStep"],
                time_control["mean"] - time_control["std"],
                time_control["mean"] + time_control["std"],
                color="#3366cc", alpha=0.3, zorder=1)
ax.plot(time_reprog["TimeStep"], time_reprog["mean"], '-o', label="Switch", color="#762a83", zorder=2)
ax.fill_between(time_reprog["TimeStep"],
                time_reprog["mean"] - time_reprog["std"],
                time_reprog["mean"] + time_reprog["std"],
                color="#762a83", alpha=0.3, zorder=1)
ax.set_xlabel("Time Step")
ax.set_ylabel("Tumor Cell Count (mean ± Std)")
ax.set_title("Final Tumor Cell Count") 
ax.legend()
plt.tight_layout()
plt.savefig("Figure_4B.png", dpi=300, bbox_inches="tight")
plt.show()

# --------------- Figure 4C: Final Tumor Cell Count Boxplot ---------------
fig, ax = plt.subplots(figsize=(6, 6))
sns.boxplot(x="GridIndex", y="TumorCellCount", data=df_base_final,
            palette=["#3366cc", "#762a83"], showfliers=False, ax=ax)
sns.stripplot(x="GridIndex", y="TumorCellCount", data=df_base_final,
              color="black", alpha=0.5, ax=ax)
ax.set_xticklabels(["No Switch", "Switch"])
ax.set_xlabel("")
ax.set_ylabel("Tumor Cell Count")
ax.set_title("Final Tumor Cell Count") 

plt.tight_layout()
plt.savefig("Figure_4C.png", dpi=300, bbox_inches="tight")
plt.show()

# --------------- Figure 4D/E/F: Morphology Boxplots ---------------
metrics_spatial = ["FractalDimension", "Eccentricity", "Lacunarity"]
label_map = {
    "FractalDimension": "Fractal Dimension",
    "Eccentricity": "Eccentricity",
    "Lacunarity": "Lacunarity"
}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
for i, metric in enumerate(metrics_spatial):
    df_metric = df_base_final[["GridIndex", metric]].dropna()
    sns.boxplot(x="GridIndex", y=metric, data=df_metric,
                palette=["#3366cc", "#762a83"], showfliers=False, ax=axes[i])
    sns.stripplot(x="GridIndex", y=metric, data=df_metric,
                  color="black", alpha=0.5, ax=axes[i])
    axes[i].set_xticklabels(["No Switch", "Switch"])
    axes[i].set_xlabel("")
    axes[i].set_ylabel(label_map[metric])
    axes[i].set_title(label_map[metric])

plt.tight_layout()
plt.savefig("Figure_4DEF.png", dpi=600, bbox_inches="tight")
plt.show()

# --------------- Figure 5: PRCC Analysis ---------------
df_psa = pd.read_csv("base_psa.csv", on_bad_lines='skip', index_col=False)

# Rename columns for clarity 
df_psa = df_psa.rename(columns={"S2": "switchSensitivity", "S4": "divisionSensitivity"})

# Filter the data to only include the final time step
last_time_psa = df_psa["TimeStep"].max()
df_psa = df_psa[df_psa["TimeStep"] == last_time_psa].copy()

# For output metrics, compute the mean - for input parameters, take the first value encountered
agg_dict = {
    "TumorCellCount": "mean",
    "FractalDimension": "mean",
    "Lacunarity": "mean",
    "Eccentricity": "mean",
    "effectAntiMet": "first",
    "effectProMet": "first",
    "conversionThreshold": "first",
    "switchSensitivity": "first",
    "divisionSensitivity": "first",
    "effectPerTumorCell": "first"
}

# Convert all columns used in the aggregation to numeric values (coercing errors to NaN)
for col in agg_dict.keys():
    df_psa[col] = pd.to_numeric(df_psa[col], errors='coerce')

# Group data by "ParamSetIndex" and aggregate
df_agg = df_psa.groupby("ParamSetIndex", as_index=False).agg(agg_dict)

# Inputs parameters for sensitivity analysis
input_cols = ["effectAntiMet", "effectProMet", "conversionThreshold",
              "switchSensitivity", "divisionSensitivity", "effectPerTumorCell"]

# Outputs metrics 
output_metrics = ["TumorCellCount", "FractalDimension", "Lacunarity", "Eccentricity"]

# Extract the input and output dataframes 
X = df_agg[input_cols].copy()
Y = df_agg[output_metrics].copy()

# Rank inputs and outputs 
X_ranked = X.rank()
Y_ranked = Y.rank()

# Create dictionary to store PRCC results for each output metric
prcc_results = {metric: {} for metric in output_metrics}

# Loop over each output metric and then each input parameter
for metric in output_metrics:
    # Skip the metric if all values are NaN in the ranked output
    if Y_ranked[metric].isna().all():
        continue
    for col in X_ranked.columns:
        # Define covariates: all input parameters except the current one
        covariates = [c for c in X_ranked.columns if c != col]
        # Copy the ranked input dataframe and add the current output metric column
        df_temp = X_ranked.copy()
        df_temp[metric] = Y_ranked[metric]
        # Compute partial correlation 
        result = pg.partial_corr(data=df_temp, x=col, y=metric, covar=covariates, method='spearman')
        # Save the computed correlation coefficient 
        prcc_results[metric][col] = result['r'].values[0]

# Convert the PRCC results into a dictionary of dataframes 
prcc_dfs = {
    metric: pd.DataFrame(list(res.items()), columns=["Parameter", "PRCC"])
    for metric, res in prcc_results.items()
}

display_names = {
    "TumorCellCount": "Tumor Cell Count",
    "FractalDimension": "Fractal Dimension",
    "Lacunarity": "Lacunarity",
    "Eccentricity": "Eccentricity"
}

fixed_order = [
    "effectPerTumorCell",
    "divisionSensitivity",
    "switchSensitivity",
    "conversionThreshold",
    "effectProMet",
    "effectAntiMet"
]

print("\n=== Partial Rank Correlation Coefficients (PRCC) ===")
for metric, df_metric in prcc_dfs.items():
    print(f"\n{display_names[metric]}:")
    # Reorder the dataframe using the fixed order and reset the index
    df_sorted = df_metric.set_index("Parameter").loc[fixed_order].reset_index()
    for _, row in df_sorted.iterrows():
        print(f"  {row['Parameter']}: {row['PRCC']:.3f}")

fig, axes = plt.subplots(1, len(output_metrics), figsize=(20, 5), sharex=True)

for i, metric in enumerate(output_metrics):
    if metric in prcc_dfs:
        df_plot = prcc_dfs[metric].set_index("Parameter").reindex(fixed_order).reset_index()
        axes[i].barh(df_plot["Parameter"], df_plot["PRCC"], color="skyblue", edgecolor="black")
        axes[i].set_title(display_names[metric], fontsize=14)
        axes[i].set_xlabel("PRCC", fontsize=12)
        axes[i].set_xlim(-1, 1)
        axes[i].axvline(0, linestyle="--", color="black")
        axes[i].tick_params(axis='y', labelsize=12)
        axes[i].tick_params(axis='x', labelsize=10)

plt.tight_layout()
plt.savefig("Figure_5_PRCC.png", dpi=300, bbox_inches="tight")
plt.show()
