import numpy as np
import pandas as pd
import scipy.stats as stats          
import statsmodels.api as sm          
import statsmodels.formula.api as smf     
import scikit_posthocs as sp            
from statsmodels.stats.multicomp import pairwise_tukeyhsd 
from scipy.stats import kruskal         
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.cluster import KMeans        
from matplotlib.colors import to_hex

# -------------------- Global Style Settings --------------------
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "axes.titlesize": 20,
    "legend.fontsize": 18
})

# -------------------- Targeted Analysis Metrics --------------------
df_spatial = pd.read_csv("spatial_targeted.csv", on_bad_lines='skip', index_col=False)
spatial_map = {
    0: "Uniform",
    1: "Random",
    2: "Clumped",
    3: "Radial",
    4: "Inverse Radial",
    5: "Gradient"
}
df_spatial["SpatialDistribution"] = df_spatial["GridIndex"].map(spatial_map)

# Filter to the final time step and drop rows missing key metrics
last_time_step = df_spatial["TimeStep"].max()
df_spatial = df_spatial[df_spatial["TimeStep"] == last_time_step].copy()
df_spatial = df_spatial.dropna(subset=["FractalDimension", "Eccentricity"])

# Define spatial metrics and names 
metrics = ["TumorCellCount", "FractalDimension", "Lacunarity", "Eccentricity"]
display_names = {
    "TumorCellCount": "Tumor Cell Count",
    "FractalDimension": "Fractal Dimension",
    "Lacunarity": "Lacunarity",
    "Eccentricity": "Eccentricity"
}

# Print summary statistics by spatial distribution 
print("=== Summary Statistics by Spatial Distribution ===")
for metric in metrics:
    print(f"\n{display_names[metric]}:")
    # Group by GridIndex and compute the mean and standard deviation.
    grouped = df_spatial.groupby("GridIndex")[metric].agg(["mean", "std"]).reset_index()
    for _, row in grouped.iterrows():
        distr = row["GridIndex"]
        mean_val = row["mean"]
        std_val = row["std"]
        print(f"  {spatial_map.get(distr, distr)}: {mean_val:.3f} ± {std_val:.3f}")

# -------------------- Statistical Analysis (Fig S7)--------------------
results = []
dunn_results_dict = {}  # to store Dunn's test DataFrames for later plotting

for metric in metrics:
    # Normality check for each spatial group (only if sample size >= 3)
    normality_p_values = []
    for distr in spatial_map.values():
        subset = df_spatial[df_spatial["SpatialDistribution"] == distr][metric]
        if len(subset) >= 3:
            _, p_val = stats.shapiro(subset)
            normality_p_values.append(p_val)
    
    # Check homogeneity of variance for groups with sufficient samples
    groups = [df_spatial[df_spatial["SpatialDistribution"] == distr][metric]
              for distr in spatial_map.values() if len(df_spatial[df_spatial["SpatialDistribution"] == distr]) >= 3]
    if groups:
        levene_stat, levene_p = stats.levene(*groups)
    else:
        levene_stat, levene_p = np.nan, np.nan

    # Choose test based on assumptions
    if normality_p_values and all(p > 0.05 for p in normality_p_values) and (not np.isnan(levene_p)) and levene_p > 0.05:
        test_used = "ANOVA"
        model = smf.ols(f"{metric} ~ C(SpatialDistribution)", data=df_spatial).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        overall_p = anova_table["PR(>F)"][0]
        posthoc = pairwise_tukeyhsd(df_spatial[metric], df_spatial["SpatialDistribution"]).summary()
    else:
        test_used = "Kruskal-Wallis"
        group_data = [df_spatial[df_spatial["SpatialDistribution"] == distr][metric] for distr in spatial_map.values()]
        valid_groups = [g for g in group_data if len(g) >= 3]
        if len(valid_groups) >= 2:
            kw_stat, kw_p = kruskal(*valid_groups)
            overall_p = kw_p
            try:
                # Perform Dunn's test with Bonferroni adjustment
                dunn_df = sp.posthoc_dunn(df_spatial, val_col=metric, group_col="SpatialDistribution", p_adjust="bonferroni")
            except ZeroDivisionError:
                dunn_df = pd.DataFrame({"Error": ["Insufficient samples for Dunn's test"]})
            posthoc = dunn_df
            dunn_results_dict[metric] = dunn_df  
        else:
            overall_p = np.nan
            posthoc = "Insufficient samples"
    
    results.append({
        "Metric": metric,
        "Test Used": test_used,
        "Overall p-value": overall_p,
        "Post-hoc Results": posthoc
    })

results_df = pd.DataFrame(results)
print("\n=== Overall Statistical Test Results ===")
print(results_df[["Metric", "Test Used", "Overall p-value"]].to_string(index=False))

dunn_keys = list(dunn_results_dict.keys())
if len(dunn_keys) > 0:
    n_plots = len(dunn_keys)
    nrows = (n_plots + 1) // 2  
    fig_dunn, axes_dunn = plt.subplots(nrows, 2, figsize=(17, 7 * nrows))
    axes_dunn = axes_dunn.flatten()
    
    alpha = 0.005 

    for idx, metric in enumerate(dunn_keys):
        dunn_df = dunn_results_dict[metric].astype(float)
        annot = dunn_df.applymap(lambda p: f"{p:.3f}*" if p < alpha else f"{p:.3f}")
        ax = axes_dunn[idx]
        sns.heatmap(dunn_df, annot=annot, fmt='', cmap="viridis", vmin=0, vmax=1, ax=ax)
        ax.set_title(f"Dunn's p-values for {display_names[metric]}")
        ax.set_xlabel("Spatial Distribution")
        ax.set_ylabel("Spatial Distribution")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    
    for j in range(idx + 1, len(axes_dunn)):
        fig_dunn.delaxes(axes_dunn[j])
    
    plt.tight_layout()
    plt.savefig("Figure_S_Dunn_Results.png", dpi=300, bbox_inches="tight")
    plt.show()
else:
    print("No Dunn's test results to plot.")

# -------------------- Figure 7B: Spatial Targeted Data Boxplots--------------------
df_box = pd.read_csv("spatial_targeted.csv", on_bad_lines='skip', index_col=False)
final_time = df_box["TimeStep"].max()
df_box = df_box[df_box["TimeStep"] == final_time].copy()

# Force "GridIndex" to be integers then convert to strings.
df_box["GridIndex"] = df_box["GridIndex"].astype(int).astype(str)

# Create a blue palette with 6 colors and assign string keys.
palette = sns.light_palette("#4C72B0", n_colors=6, reverse=False)
palette_dict = {str(i): to_hex(color) for i, color in enumerate(palette)}

# Distribution labels will be displayed on the x-axis.
distribution_labels = ["Uniform", "Random", "Clumped", "Radial", "Inverse Radial", "Gradient"]

# Define metrics and their display names.
metrics = ["TumorCellCount", "FractalDimension", "Lacunarity", "Eccentricity"]
display_names = {
    "TumorCellCount": "Tumor Cell Count",
    "FractalDimension": "Fractal Dimension",
    "Lacunarity": "Lacunarity",
    "Eccentricity": "Eccentricity"
}

# Create subplots 
fig, axes = plt.subplots(1, 4, figsize=(20, 6), sharey=False)

for ax, metric in zip(axes, metrics):
    sns.boxplot(
        x="GridIndex",
        y=metric,
        data=df_box,
        palette=palette_dict,          
        order=["0", "1", "2", "3", "4", "5"], 
        width=0.5,
        showfliers=False,
        ax=ax
    )
    sns.stripplot(
        x="GridIndex",
        y=metric,
        data=df_box,
        color="black",
        size=3,
        jitter=True,
        order=["0", "1", "2", "3", "4", "5"],
        ax=ax
    )
    ax.set_title(display_names[metric], fontsize=20)
    ax.set_ylabel(display_names[metric], fontsize=18)
    ax.set_xlabel("")  
    ax.set_xticks(np.arange(len(distribution_labels)))
    ax.set_xticklabels(distribution_labels, rotation=40, ha="right", fontsize=12)
    ax.tick_params(axis="both", labelsize=18)

plt.tight_layout()
plt.savefig("Figure_7A_Spatial_Targeted_Boxplots.png", dpi=600, bbox_inches="tight")
plt.show()

# -------------------- Figure 7C: Spatial PSA Boxplots--------------------
df_psa = pd.read_csv("spatial_psa.csv", on_bad_lines='skip')
df_psa = df_psa.reset_index().rename(columns={"index": "ParamSetIndex"})

# Define desired columns and assign new column names
desired_columns = [
    "ParamSetIndex", "Replication", "TimeStep", "GridIndex", 
    "effectAntiMet", "effectProMet", "conversionThreshold", 
    "S2", "S4", "effectPerTumorCell", "TumorCellCount", 
    "FractalDimension", "Lacunarity", "Eccentricity"
]
df_psa = df_psa.iloc[:, :len(desired_columns)]
df_psa.columns = desired_columns

# Rename sensitivity columns for clarity
df_psa = df_psa.rename(columns={"S2": "switchSensitivity", "S4": "divisionSensitivity"})

# Keep only the data from the final time step
final_time = df_psa["TimeStep"].max()
df_final = df_psa[df_psa["TimeStep"] == final_time].copy()

# Require a minimum TumorCellCount 
df_final = df_final[df_final["TumorCellCount"] >= 3000].copy()

# Define regimes via tertile split
df_ref = df_final[df_final["GridIndex"] == 1].copy()

# Group by ParamSetIndex and compute average TumorCellCount for the reference condition.
df_ranked = df_ref.groupby("ParamSetIndex")["TumorCellCount"].mean().reset_index()
df_ranked = df_ranked.sort_values("TumorCellCount").reset_index(drop=True)

# Calculate the number of samples and the size of one tertile.
n = len(df_ranked)
third = n // 3

# Label the bottom tertile as "Inhibitory", middle as "Neutral", and top as "Promoting"
df_ranked["Regime"] = (["Inhibitory"] * third +
                       ["Neutral"] * third +
                       ["Promoting"] * (n - 2 * third))

# Keep only "Inhibitory" and "Promoting" regimes
df_ranked = df_ranked[df_ranked["Regime"].isin(["Inhibitory", "Promoting"])]

# Merge the regime assignments back into the full PSA dataframe.
df_final = df_final.merge(df_ranked[["ParamSetIndex", "Regime"]], on="ParamSetIndex", how="inner")

# Create heatmap matrices for each metric
grid_indices = [0, 1, 2, 3, 4, 5]    
regimes = ["Inhibitory", "Promoting"]  

# Define the metrics and their display names
metrics = ["TumorCellCount", "FractalDimension", "Lacunarity", "Eccentricity"]
display_names = {
    "TumorCellCount": "Tumor Cell Count",
    "FractalDimension": "Fractal Dimension",
    "Lacunarity": "Lacunarity",
    "Eccentricity": "Eccentricity"
}

# Initialize a dictionary to store heatmap matrices
heatmap_matrices = {}
for metric in metrics:
    # Create an empty DataFrame with regimes as rows and grid_indices as columns
    heat_data = pd.DataFrame(index=regimes, columns=grid_indices)
    # For each grid index and for each regime, compute the mean of the given metric.
    for g in grid_indices:
        subset = df_final[df_final["GridIndex"] == g]
        for reg in regimes:
            mean_val = subset[subset["Regime"] == reg][metric].mean()
            heat_data.loc[reg, g] = mean_val
    heatmap_matrices[metric] = heat_data.astype(float)

# Define x-axis labels for the heatmaps.
distribution_labels = ["Uniform", "Random", "Clumped", "Radial", "Inverse Radial", "Gradient"]

# Create subplots for each metric 
fig, axes = plt.subplots(1, 4, figsize=(20, 4), sharey=True)
for i, metric in enumerate(metrics):
    ax = axes[i]
    heatmap_matrix = heatmap_matrices[metric].copy()
    heatmap_matrix.columns = distribution_labels
    sns.heatmap(
        heatmap_matrix,
        annot=False,   
        fmt=".2f",
        cmap="Blues",
        cbar=True,
        cbar_kws={'shrink': 0.9, 'aspect': 30, 'pad': 0.02},
        square=False,
        ax=ax
    )
    ax.set_title(display_names[metric], fontsize=20)
    ax.set_xticks(np.arange(len(distribution_labels)) + 0.5)
    ax.set_xticklabels(distribution_labels, rotation=45, ha="right", fontsize=18)
    ax.tick_params(axis="y", rotation=0, labelsize=18)

plt.tight_layout()
plt.savefig("Figure_7B_Spatial_PSA_Heatmaps.png", format="png", dpi=600, bbox_inches="tight")
plt.show()
