import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Read eigenvectors
pcs = pd.read_csv("results/pca/chr22_pca.eigenvec", sep="\t")
pcs = pcs.rename(columns={"#IID": "IID"})

# Get population labels from panel file
panel = pd.read_csv(
    "data/raw/integrated_call_samples_v3.20130502.ALL.panel",
    sep="\t", usecols=[0, 1, 2], header=0,
    names=["IID", "pop", "super_pop"],
)
df = pcs.merge(panel, on="IID", how="left")

# Read outliers from the heterozygosity step
out = pd.read_csv("results/qc/het_outliers_pruned.csv")
df["het_outlier"] = df["IID"].isin(out["IID"])

colors = {"AFR": "tab:orange", "AMR": "tab:green", "EAS": "tab:red",
          "EUR": "tab:blue", "SAS": "tab:purple"}

# Plot PCA
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (a, b) in zip(axes, [("PC1", "PC2"), ("PC3", "PC4")]):
    for sp, c in colors.items():
        sub = df[df["super_pop"] == sp]
        ax.scatter(sub[a], sub[b], s=6, c=c, label=sp, alpha=0.7)
    flagged = df[df["het_outlier"]]
    ax.scatter(flagged[a], flagged[b], s=50, facecolors="none",
               edgecolors="black", label="het outlier")
    ax.set_xlabel(a)
    ax.set_ylabel(b)
axes[0].legend(markerscale=2)
plt.tight_layout()
plt.savefig("results/pca/pca_by_superpop.png", dpi=150)
print(df["super_pop"].value_counts())

# Plot PCA for AMR populations
amr = df[df["super_pop"] == "AMR"]
fig, ax = plt.subplots(figsize=(6, 5))
for p, sub in amr.groupby("pop"):
    ax.scatter(sub["PC3"], sub["PC4"], s=10, label=p)
ax.set_xlabel("PC3")
ax.set_ylabel("PC4")
ax.legend()
plt.tight_layout()
plt.savefig("results/pca/pca_amr_pop.png", dpi=150)