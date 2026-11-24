import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

THRESH = 5e-8
RUNS = {
    "no covariates": "results/assoc/chr22_nocov.trait.glm.linear",
    "with PC1-PC10": "results/assoc/chr22_pcs.trait.glm.linear",
}

# Read association results
def load(path):
    df = pd.read_csv(path, sep="\t").rename(columns={"#CHROM": "CHROM"})
    df = df[df["TEST"] == "ADD"].copy()
    df = df.dropna(subset=["P"])
    df["logp"] = -np.log10(df["P"].clip(lower=1e-300))
    return df

# Calculate genomic inflation factor (lambda)
def genomic_lambda(p):
    return np.median(chi2.isf(p, df=1)) / chi2.ppf(0.5, df=1)

results = {name: load(path) for name, path in RUNS.items()}

# Read causal variants
key = pd.read_csv("data/processed/causal_key.csv")
causal_pos = key["variant"].str.split(":").str[1].astype(int)

# Estimate association values
lams = {}
for name, df in results.items():
    lams[name] = genomic_lambda(df["P"].to_numpy())
    sig = df.loc[df["P"] < THRESH, "POS"].sort_values().to_numpy()
    loci = 0 if len(sig) == 0 else 1 + int((np.diff(sig) > 250_000).sum())
    print(f"\n== {name} ==")
    print(f"variants tested: {len(df)}")
    print(f"lambda: {lams[name]:.3f}")
    print(f"significant variants (P < {THRESH:g}): {len(sig)} in about {loci} locus(es); true causal variants: {len(key)}")

# Compare estimates with key
for name, df in results.items():
    m = key.merge(df[["ID", "A1", "BETA", "SE", "P"]],
                  left_on="variant", right_on="ID", how="left")
    flip = np.where(m["A1"] == m["counted_allele"], 1.0, -1.0)
    m["beta_hat"] = m["BETA"] * flip
    m["z_diff"] = (m["beta_hat"] - m["beta_per_allele"]) / m["SE"]
    m["detected"] = m["P"] < THRESH
    cols = ["variant", "counted_freq", "var_explained", "beta_per_allele",
            "beta_hat", "SE", "z_diff", "P", "detected"]
    print(f"\n-- causal variants, {name} --")
    print(m[cols].round(4).to_string(index=False))

# Produce Manhattan plot
fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True, sharey=True)
for ax, (name, df) in zip(axes, results.items()):
    for pos in causal_pos:
        ax.axvline(pos / 1e6, color="orange", lw=1, alpha=0.7, zorder=1)
    ax.scatter(df["POS"] / 1e6, df["logp"], s=3, c="tab:blue", alpha=0.6, zorder=2)
    ax.axhline(-np.log10(THRESH), color="red", ls="--", lw=1)
    ax.set_ylabel("-log10(P)")
    ax.set_title(f"{name}  (lambda = {lams[name]:.3f})")
axes[-1].set_xlabel("chr22 position (Mb, GRCh37)")
plt.tight_layout()
plt.savefig("results/assoc/manhattan.png", dpi=150)

# Produce QQ plot
fig, ax = plt.subplots(figsize=(5, 5))
for name, df in results.items():
    p = np.sort(df["P"].to_numpy())
    n = len(p)
    obs = -np.log10(np.clip(p, 1e-300, None))
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    ax.scatter(exp, obs, s=3, label=name, alpha=0.7)
lim = ax.get_xlim()[1]
ax.plot([0, lim], [0, lim], color="grey", lw=1)
ax.set_xlabel("expected -log10(P)")
ax.set_ylabel("observed -log10(P)")
ax.legend()
plt.tight_layout()
plt.savefig("results/assoc/qq.png", dpi=150)