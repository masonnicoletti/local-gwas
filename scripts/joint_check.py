import numpy as np
import pandas as pd

raw = pd.read_csv("data/processed/causal_genotypes.raw", sep=r"\s+").rename(columns={"#IID": "IID"})
pheno = pd.read_csv("data/processed/pheno.txt", sep="\t").rename(columns={"#IID": "IID"})
pcs = pd.read_csv("results/pca/chr22_pca.eigenvec", sep="\t").rename(columns={"#IID": "IID"})
df = raw.merge(pheno, on="IID").merge(pcs, on="IID")

gcols = [c for c in raw.columns if ":" in c]
pccols = [f"PC{i}" for i in range(1, 11)]
y = df["trait"].to_numpy()

def ols(X):
    X1 = np.column_stack([np.ones(len(X)), X])
    return np.linalg.lstsq(X1, y, rcond=None)[0]

key = pd.read_csv("data/processed/causal_key.csv")
key["joint"] = ols(df[gcols].to_numpy(float))[1:11]
key["joint_with_PCs"] = ols(df[gcols + pccols].to_numpy(float))[1:11]
print(key[["variant", "beta_per_allele", "joint", "joint_with_PCs"]].round(4).to_string(index=False))