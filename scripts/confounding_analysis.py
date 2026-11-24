import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

THRESH = 5e-8
RUNS = {
    "no covariates": "results/assoc/conf_nocov.trait_conf.glm.linear",
    "with PC1-PC10": "results/assoc/conf_pcs.trait_conf.glm.linear",
}

# Read association results from confounding analysis
def load(path):
    df = pd.read_csv(path, sep="\t").rename(columns={"#CHROM": "CHROM"})
    df = df[df["TEST"] == "ADD"].dropna(subset=["P"]).copy()
    df["logp"] = -np.log10(df["P"].clip(lower=1e-300))
    return df

# Calculate genomic inflation factor (lambda)
def lam(p):
    return np.median(chi2.isf(np.asarray(p), df=1)) / chi2.ppf(0.5, df=1)

results = {name: load(path) for name, path in RUNS.items()}
lams = {}

# Estimate association statistics
for name, df in results.items():
    lams[name] = lam(df["P"])
    sig = df[df["P"] < THRESH]
    print(f"\n== {name} ==")
    print(f"lambda: {lams[name]:.3f}")
    print(f"significant variants: {len(sig)} "
          f"(expected by chance: {len(df) * THRESH:.3f}); every one is a false positive")
    print(f"smallest P: {df['P'].min():.2e}")
    print(f"share of variants with P < 0.05: {(df['P'] < 0.05).mean():.3f} (chance: 0.050)")

# Manhattan plot
fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True, sharey=True)
for ax, (name, df) in zip(axes, results.items()):
    ax.scatter(df["POS"] / 1e6, df["logp"], s=3, c="tab:blue", alpha=0.6)
    ax.axhline(-np.log10(THRESH), color="red", ls="--", lw=1)
    ax.set_ylabel("-log10(P)")
    ax.set_title(f"{name}  (lambda = {lams[name]:.3f})")
axes[-1].set_xlabel("chr22 position (Mb, GRCh37)")
plt.tight_layout()
plt.savefig("results/assoc/conf_manhattan.png", dpi=150)

# QQ plot
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
plt.savefig("results/assoc/conf_qq.png", dpi=150)