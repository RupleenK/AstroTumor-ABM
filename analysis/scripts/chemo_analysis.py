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

# ---------------- Figure 8A: Chemo Response Time Series Analysis ----------------
df_chemo = pd.read_csv("chemo_targeted.csv")
# Group data by condition and time step, then compute mean, standard deviation and count.
grouped = df_chemo.groupby(["GridIndex", "TimeStep"])["TumorCellCount"].agg(["mean", "std", "count"]).reset_index()
# Compute standard error of the mean (SEM)
grouped["sem"] = grouped["std"] / np.sqrt(grouped["count"])

fig, ax = plt.subplots(figsize=(8, 6))
colors = {0: "#377eb8", 1: "#984ea3"}
for grid in [0, 1]:
    sub = grouped[grouped["GridIndex"] == grid].sort_values("TimeStep")
    label = "No Switch" if grid == 0 else "Switch"
    ax.plot(sub["TimeStep"], sub["mean"], color=colors[grid], label=label, linewidth=3)
    ax.fill_between(sub["TimeStep"],
                    sub["mean"] - sub["sem"],
                    sub["mean"] + sub["sem"],
                    color=colors[grid], alpha=0.3)
ax.set_xlabel("Time Step", fontsize=16)
ax.set_ylabel("Tumor Cell Count (Mean ± SEM)", fontsize=16)
ax.set_title("Chemotherapy Response", fontsize=16)
ax.legend(fontsize=16)
plt.tight_layout()
plt.savefig("Chemo_Targeted.png", dpi=300, bbox_inches="tight")
plt.show()

# ---------------- Time to Extinction and Other Metrics ----------------
# Ensure that the "Replication" column is treated as an integer.
df_chemo["Replicate"] = df_chemo["Replication"].astype(int)
final_time = df_chemo["TimeStep"].max()
conditions = [0, 1] 

results = {}

for cond in conditions:
    sub = df_chemo[df_chemo["GridIndex"] == cond].copy()
    total_runs = len(sub["Replicate"].unique())
    extinct_count = 0
    extinction_times = []     # Collect extinction times for replicates that went extinct
    final_sizes_survivors = [] # Store final tumor sizes for replicates that survive
    
    # Loop over each replicate
    for rep in sub["Replicate"].unique():
        rep_data = sub[sub["Replicate"] == rep].copy()
        rep_data.sort_values("TimeStep", inplace=True)
        
        # Use a cutoff time
        rep_data_after_cutoff = rep_data[rep_data["TimeStep"] >= 35]
        is_extinct = False
        
        if not rep_data_after_cutoff.empty:
            # Find the first time when TumorCellCount drops below 1000 (adjust threshold as needed)
            low_steps = rep_data_after_cutoff[rep_data_after_cutoff["TumorCellCount"] < 1000]["TimeStep"]
            if not low_steps.empty:
                first_low_time = low_steps.iloc[0]
                # If tumor count stays below 1000 for all subsequent timesteps, mark as extinct.
                subsequent_data = rep_data_after_cutoff[rep_data_after_cutoff["TimeStep"] >= first_low_time]
                if (subsequent_data["TumorCellCount"] < 1000).all():
                    extinct_count += 1
                    extinction_times.append(first_low_time)
                    is_extinct = True
        
        # If not extinct, record the final tumor count at final_time
        if not is_extinct:
            rep_final = rep_data[rep_data["TimeStep"] == final_time]
            if not rep_final.empty:
                final_tumor_count = rep_final["TumorCellCount"].iloc[0]
                final_sizes_survivors.append(final_tumor_count)
    
    # Compute fraction of replicates that went extinct
    extinct_fraction = (extinct_count / total_runs) * 100
    
    # Calculate mean ± std extinction time among extinct replicates
    if len(extinction_times) > 0:
        mean_ext_time = np.mean(extinction_times)
        std_ext_time = np.std(extinction_times)
    else:
        mean_ext_time = np.nan
        std_ext_time = np.nan
    
    # Calculate mean ± std of final tumor size among survivors
    if len(final_sizes_survivors) > 0:
        mean_final_size = np.mean(final_sizes_survivors)
        std_final_size = np.std(final_sizes_survivors)
    else:
        mean_final_size = np.nan
        std_final_size = np.nan
    
    results[cond] = {
        "Condition": "No Switch" if cond == 0 else "Switch",
        "TotalRuns": total_runs,
        "ExtinctCount": extinct_count,
        "ExtinctFraction": extinct_fraction,
        "MeanExtTime": mean_ext_time,
        "StdExtTime": std_ext_time,
        "MeanFinalSize": mean_final_size,
        "StdFinalSize": std_final_size
    }

# Print summary ---
for cond in conditions:
    r = results[cond]
    print(f"\nCondition = {r['Condition']} (GridIndex={cond})")
    print(f" - Total replicates: {r['TotalRuns']}")
    print(f" - Extinct replicates: {r['ExtinctCount']} ({r['ExtinctFraction']:.1f}%)")
    if not np.isnan(r["MeanExtTime"]):
        print(f" - Mean extinction time: {r['MeanExtTime']:.1f} ± {r['StdExtTime']:.1f}")
    if not np.isnan(r["MeanFinalSize"]):
        print(f" - Mean final size (survivors): {r['MeanFinalSize']:.1f} ± {r['StdFinalSize']:.1f}")
    else:
        print(" - No survivors or no data at final time.")
