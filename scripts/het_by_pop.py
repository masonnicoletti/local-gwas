import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# PLINK's per-sample heterozygosity results
het = pd.read_csv("results/qc/chr22_varqc_samples.het", sep="\t")
het = het.rename(columns={"#IID": "IID"})

# Panel file with population information
panel = pd.read_csv(
    "data/raw/integrated_call_samples_v3.20130502.ALL.panel",
    sep="\t", usecols=[0, 1, 2, 3], header=0,
    names=["IID", "pop", "super_pop", "sex"],
)

# Merge on samples
df = het.merge(panel, on="IID", how="left")
print(df["pop"].isna().sum(), "samples without a panel match (expect 0)")

# Group by ancestry groups
print(df.groupby("super_pop")["F"].agg(["count", "mean", "std"]).round(3))

# Outliers within populations
df["z"] = df.groupby("pop")["F"].transform(lambda x: (x - x.mean()) / x.std())
print(df.loc[df["z"].abs() > 4, ["IID", "pop", "F", "z"]])

# Boxplot of F by population, grouped by super-population
pops = df.sort_values(["super_pop", "pop"])["pop"].unique()
fig, ax = plt.subplots(figsize=(12, 4))
ax.boxplot([df.loc[df["pop"] == p, "F"] for p in pops])
ax.set_xticks(range(1, len(pops) + 1))
ax.set_xticklabels(pops, rotation=90)
ax.axhline(0, color="grey", lw=0.5)
ax.set_ylabel("F (inbreeding coefficient)")
plt.tight_layout()
plt.savefig("results/qc/het_by_pop.png", dpi=150)

# Save outliers to CSV
df.loc[df["z"].abs() > 4, ["IID", "pop", "super_pop", "F", "z"]].to_csv("results/qc/het_outliers.csv", index=False)