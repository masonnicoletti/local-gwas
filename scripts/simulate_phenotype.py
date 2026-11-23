import numpy as np
import pandas as pd

SEED = 45
H2 = 0.30

# Read genotypes
raw = pd.read_csv("data/processed/causal_genotypes.raw", sep=r"\s+")
raw = raw.rename(columns={"#IID": "IID"})
geno_cols = [c for c in raw.columns if ":" in c]
G = raw[geno_cols].to_numpy(dtype=float)
n, m = G.shape
assert not np.isnan(G).any(), "missing genotypes found"

# Standardize each variant
mean = G.mean(axis=0)
sd = G.std(axis=0)
Z = (G - mean) / sd

rng = np.random.default_rng(SEED)

# Random effect sizes
beta_std = rng.normal(size=m)
beta_std *= np.sqrt(H2 / (Z @ beta_std).var())
genetic = Z @ beta_std

# Add noise
noise = rng.normal(0, np.sqrt(1 - H2), size=n)
y = genetic + noise

# Create key for causal variants
key = pd.DataFrame({
    "variant": [c.rsplit("_", 1)[0] for c in geno_cols],
    "counted_allele": [c.rsplit("_", 1)[1] for c in geno_cols],
    "counted_freq": mean / 2,
    "beta_per_allele": beta_std / sd,
    "var_explained": beta_std ** 2,
})
key.to_csv("data/processed/causal_key.csv", index=False)

pd.DataFrame({"#IID": raw["IID"], "trait": y}).to_csv(
    "data/processed/pheno.txt", sep="\t", index=False)

print(n, "samples,", m, "causal variants")
print("var(genetic) / var(trait) =", round(genetic.var() / y.var(), 3))
print("sum of per-variant var_explained =", round(key["var_explained"].sum(), 3))
print(key.round(4).to_string(index=False))