import numpy as np
import pandas as pd
from scipy.stats import chi2

THRESH = 5e-8
WINDOW = 500_000

RUNS = {
    "no covariates": "results/assoc/chr22_nocov.trait.glm.linear",
    "with PC1-PC10": "results/assoc/chr22_pcs.trait.glm.linear",
}

key = pd.read_csv("data/processed/causal_key.csv")
causal_pos = key["variant"].str.split(":").str[1].astype(int).to_numpy()

def lam(p):
    return np.median(chi2.isf(np.asarray(p), df=1)) / chi2.ppf(0.5, df=1)

for name, path in RUNS.items():
    df = pd.read_csv(path, sep="\t")
    df = df[df["TEST"] == "ADD"].dropna(subset=["P"]).copy()

    # Calculate distance to nearest causal variant
    df["dist"] = np.abs(df["POS"].to_numpy()[:, None] - causal_pos[None, :]).min(axis=1)
    near = df["dist"] <= WINDOW
    sig = df["P"] < THRESH

    print(f"\n== {name} ==")
    print(f"lambda, all variants:           {lam(df['P']):.3f}")
    print(f"lambda, outside causal windows: {lam(df.loc[~near, 'P']):.3f}  ({(~near).sum()} variants)")
    print(f"significant near a causal variant (<= {WINDOW // 1000} kb): {(sig & near).sum()}")
    print(f"significant far from every causal variant:      {(sig & ~near).sum()}")
    print(df.loc[sig & ~near, ["ID", "POS", "dist", "P"]].sort_values("P").head(10).to_string(index=False))

    # Report P-values for causal variants
    m = key[["variant"]].merge(df[["ID", "P"]], left_on="variant", right_on="ID")
    m["P"] = m["P"].map(lambda x: f"{x:.2e}")
    print(m[["variant", "P"]].to_string(index=False))