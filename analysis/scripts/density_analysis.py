import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats  
import statsmodels.formula.api as smf
import statsmodels.api as sm
import scikit_posthocs as sp 
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# ---------------- Global Style Settings ---------------
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "axes.titlesize": 20,
    "legend.fontsize": 18
})

# ---------------- Density Targeted Statistical Analysis ----------------
df = pd.read_csv("density_targeted.csv")

# Map GridIndex values to astrocyte density percentages
astrocyte_density_map = {0: 0, 1: 10, 2: 20, 3: 30, 4: 40, 5: 50}
df["AstrocyteDensity"] = df["GridIndex"].map(astrocyte_density_map)

# Keep only the final time step and drop rows missing key morphology metrics
last_time = df["TimeStep"].max()
df = df[df["TimeStep"] == last_time].copy().dropna(subset=["FractalDimension", "Eccentricity"])

metrics = ["TumorCellCount", "FractalDimension", "Lacunarity", "Eccentricity"]
results = []

# Loop through each metric and perform tests based on data normality and variance
for metric in metrics:
    normality_p_values = []
    # Check normality for each density level (if sample size is sufficient)
    for density in astrocyte_density_map.values():
        subset = df[df["AstrocyteDensity"] == density][metric]
        if len(subset) >= 3:
            _, p_val = stats.shapiro(subset)
            normality_p_values.append(p_val)
          
    # Perform Levene's test for homogeneity of variances
    levene_stat, levene_p = stats.levene(*[df[df["AstrocyteDensity"] == d][metric]
                                            for d in astrocyte_density_map.values()])
    
    # If data across density levels is normally distributed and variances are equal, use ANOVA
    if all(p > 0.05 for p in normality_p_values) and levene_p > 0.05:
        test_used = "ANOVA"
        model = smf.ols(f"{metric} ~ C(AstrocyteDensity)", data=df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        tukey_results = pairwise_tukeyhsd(df[metric], df["AstrocyteDensity"]).summary()
        results.append({
            "Metric": metric,
            "Test Used": test_used,
            "ANOVA p-value": anova_table["PR(>F)"][0],
            "Post-hoc Test": "Tukey HSD",
            "Post-hoc Results": tukey_results
        })
    else:
        test_used = "Kruskal-Wallis"
        kw_stat, kw_p = stats.kruskal(*[df[df["AstrocyteDensity"] == d][metric]
                                        for d in astrocyte_density_map.values()])
        dunn_results = sp.posthoc_dunn(df, val_col=metric, group_col="AstrocyteDensity", p_adjust="bonferroni")
        results.append({
            "Metric": metric,
            "Test Used": test_used,
            "Kruskal-Wallis H-statistic": kw_stat,
            "Kruskal-Wallis p-value": kw_p,
            "Post-hoc Test": "Dunn’s test",
            "Post-hoc Results": dunn_results
        })
      
# Display Statistical Analysis Results
results_df = pd.DataFrame(results)
print("=== Density Targeted Statistical Analysis ===\n")
for _, row in results_df.iterrows():
    print(f"Metric: {row['Metric']}")
    if row["Test Used"] == "ANOVA":
        print(f"  {row['Test Used']} p-value: {row['ANOVA p-value']:.4g}")
    else:
        print(f"  {row['Test Used']} H-statistic: {row['Kruskal-Wallis H-statistic']:.4g}, "
              f"p-value: {row['Kruskal-Wallis p-value']:.4g}")
    print(f"  Post-hoc Test: {row['Post-hoc Test']}")
    print(f"  Post-hoc Results:\n{row['Post-hoc Results']}\n")

# ---------------- Figure 6B: Density Targeted Box Plots ----------------
df_density = pd.read_csv("density_targeted.csv")
last_time_density = df_density["TimeStep"].max()
df_density = df_density[df_density["TimeStep"] == last_time_density].copy()
metrics_density = ["TumorCellCount", "FractalDimension", "Lacunarity", "Eccentricity"]
disp_names = {"TumorCellCount": "Tumor Cell Count",
              "FractalDimension": "Fractal Dimension",
              "Lacunarity": "Lacunarity",
              "Eccentricity": "Eccentricity"}

# Generate a palette for density levels
palette = sns.light_palette("#4C72B0", n_colors=6, reverse=False)
palette_dict = {i: color for i, color in enumerate(palette)}

# Create subplots for each metric
fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharex=True)
for ax, metric in zip(axes, metrics_density):
    # Drop missing values for the current metric
    df_met = df_density[["GridIndex", metric]].dropna()
    # Create boxplots with a strip overlay for detailed distribution
    sns.boxplot(x="GridIndex", y=metric, data=df_met, hue="GridIndex",
                palette=palette_dict, showfliers=False, width=0.5, legend=False, ax=ax)
    sns.stripplot(x="GridIndex", y=metric, data=df_met, color="black", size=3, jitter=True, ax=ax)
    ax.set_title(disp_names[metric], fontsize=16)
    ax.set_xticks(np.arange(6))
    ax.set_xticklabels(["0%", "10%", "20%", "30%", "40%", "50%"], fontsize=16)
    ax.set_xlabel("")
    ax.tick_params(axis='y', labelsize=18)
    ax.set_ylabel(disp_names[metric], fontsize=16)
plt.tight_layout()
plt.savefig("Density_Targeted.png", dpi=600, bbox_inches="tight")
plt.show()

# ----------- Figure 6C: Density PSA Analysis Box Plots -----------
df_psa = pd.read_csv("density_psa.csv", on_bad_lines='skip')
df_psa = df_psa.reset_index().rename(columns={"index": "ParamSetIndex"})

desired_columns = [
    "ParamSetIndex", "Replication", "TimeStep", "GridIndex", 
    "effectAntiMet", "effectProMet", "conversionThreshold", 
    "S2", "S4", "effectPerTumorCell", "TumorCellCount", 
    "FractalDimension", "Lacunarity", "Eccentricity"
]
df_psa = df_psa.iloc[:, :len(desired_columns)]
df_psa.columns = desired_columns
df_psa = df_psa.rename(columns={"S2": "switchSensitivity", "S4": "divisionSensitivity"})

# Keep only the final time step
final_time = df_psa["TimeStep"].max()
df_final = df_psa[df_psa["TimeStep"] == final_time].copy()

# Define regimes based on tumor cell count for a 10% density condition 
# For GridIndex==1 (10%), group by ParamSetIndex and compute the average TumorCellCount
df_10 = df_final[df_final["GridIndex"] == 1]
df_ranked = df_10.groupby("ParamSetIndex")["TumorCellCount"].mean().reset_index()
df_ranked = df_ranked.sort_values("TumorCellCount").reset_index(drop=True)
n = len(df_ranked)
third = n // 3

# Label the bottom tertile as "Inhibitory", middle as "Neutral", and top as "Promoting"
df_ranked["Regime"] = (["Inhibitory"] * third +
                       ["Neutral"] * third +
                       ["Promoting"] * (n - 2 * third))

# Keep only "Inhibitory" and "Promoting" regimes
selected_params = df_ranked[df_ranked["Regime"].isin(["Inhibitory", "Promoting"])]["ParamSetIndex"]

# Merge the regime labels back into the final dataframe and subset to relevant GridIndex values
df_final = df_final.merge(df_ranked[["ParamSetIndex", "Regime"]], on="ParamSetIndex", how="inner")

# Only consider GridIndex 1 (10%) and 5 (50%) for plotting
df_plot = df_final[df_final["GridIndex"].isin([1, 5])].copy()
df_plot = df_plot[df_plot["Regime"].isin(["Inhibitory", "Promoting"])].copy()

# Map GridIndex to display labels ("10%" or "50%")
df_plot["GridLabel"] = df_plot["GridIndex"].apply(lambda ix: "10%" if ix == 1 else "50%")

# Define metrics to plot and their display names
metrics = ["TumorCellCount", "FractalDimension", "Lacunarity", "Eccentricity"]
display_names = {
    "TumorCellCount": "Tumor Cell Count",
    "FractalDimension": "Fractal Dimension",
    "Lacunarity": "Lacunarity",
    "Eccentricity": "Eccentricity"
}

# Define custom colors for regimes
color_dict = {
    "Inhibitory": "#4C72B0",  # Blue
    "Promoting": "#8172B3"    # Purple
}

# Define x-axis positions for each combination of GridLabel and Regime
xpos = {
    ("10%", "Inhibitory"): 0 - 0.15,
    ("10%", "Promoting"):  0 + 0.15,
    ("50%", "Inhibitory"): 1 - 0.15,
    ("50%", "Promoting"):  1 + 0.15
}

fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharey=False)

for ax, metric in zip(axes, metrics):
    # Loop through each combination of GridLabel and Regime 
    for (gridlbl, regime), x_pos in xpos.items():
        data = df_plot[(df_plot["GridLabel"] == gridlbl) & (df_plot["Regime"] == regime)][metric].dropna()
        # Plot a boxplot for the current dataset
        bplot = ax.boxplot(
            data.values,
            positions=[x_pos],
            widths=0.25,
            patch_artist=True,
            showfliers=True
        )
        # Customize box appearance for each regime
        for box in bplot["boxes"]:
            box.set_facecolor(color_dict[regime])
            box.set_edgecolor("black")
        for whisker in bplot["whiskers"]:
            whisker.set_color("black")
        for cap in bplot["caps"]:
            cap.set_color("black")
        for median in bplot["medians"]:
            median.set_color("black")
        # Customize outlier markers
        for flier in bplot["fliers"]:
            flier.set(marker='o', alpha=0.4, markerfacecolor='black', markersize=4)
    
    # Set x-axis labels and title for each subplot
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["10%", "50%"], fontsize=16)
    ax.grid(False, axis='x')
    ax.set_title(display_names[metric], fontsize=16)
    ax.set_ylabel(display_names[metric], fontsize=16)
    ax.tick_params(axis='y', labelsize=16)

plt.tight_layout(rect=[0, 0, 0.95, 1])
plt.savefig("manual_boxplot.png", dpi=300, bbox_inches="tight")
plt.show()
