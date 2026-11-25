#!/usr/bin/env python3
"""Summarize two GWAS runs (no covariates vs with PCs): lambda, significant
variants, Manhattan and QQ plots, and (if a key is given) a comparison with
the causal-variant answer key."""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

THRESH = 5e-8

ap = argparse.ArgumentParser()
ap.add_argument("--nocov", required=True)
ap.add_argument("--pcs", required=True)
ap.add_argument("--key", default=None)
ap.add_argument("--prefix", default="assoc")
args = ap.parse_args()

# Read pipeline results
def load(path):
    df = pd.read_csv(path, sep="\t").rename(columns={"#CHROM": "CHROM"})
    df = df[df["TEST"] == "ADD"].dropna(subset=["P"]).copy()
    df["logp"] = -np.log10(df["P"].clip(lower=1e-300))
    return df

# Calculate genomic inflation factor (lambda)
def lam(p):
    return float(np.median(chi2.isf(np.asarray(p), df=1)) / chi2.ppf(0.5, df=1))

# Count loci
def count_loci(positions, gap=250_000):
    pos = np.sort(positions)
    return 0 if len(pos) == 0 else 1 + int((np.diff(pos) > gap).sum())

results = {"no covariates": load(args.nocov), "with PCs": load(args.pcs)}

key = pd.read_csv(args.key) if args.key else None
if key is not None:
    causal_pos = key["variant"].str.split(":").str[1].astype(int)
else:
    causal_pos = []

# Summary statistics for runs
summary = []
lams = {}
for name, df in results.items():
    lams[name] = lam(df["P"])
    sig = df[df["P"] < THRESH]
    summary.append({
        "run": name,
        "variants_tested": len(df),
        "lambda": round(lams[name], 3),
        "significant": len(sig),
        "loci": count_loci(sig["POS"].to_numpy()),
    })
summary = pd.DataFrame(summary)
summary.to_csv(f"{args.prefix}_summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False))

# Compare to causal variants
if key is not None:
    rows = []
    for name, df in results.items():
        m = key.merge(df[["ID", "A1", "BETA", "SE", "P"]],
                      left_on="variant", right_on="ID", how="left")
        flip = np.where(m["A1"] == m["counted_allele"], 1.0, -1.0)
        m["beta_hat"] = m["BETA"] * flip
        m["run"] = name
        m["detected"] = m["P"] < THRESH
        rows.append(m[["run", "variant", "counted_freq", "var_explained",
                       "beta_per_allele", "beta_hat", "SE", "P", "detected"]])
    causal = pd.concat(rows)
    causal.to_csv(f"{args.prefix}_causal.tsv", sep="\t", index=False)
    print(causal.groupby("run")["detected"].sum().to_string())

# Manhattan plot
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
plt.savefig(f"{args.prefix}_manhattan.png", dpi=150)

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
plt.savefig(f"{args.prefix}_qq.png", dpi=150)