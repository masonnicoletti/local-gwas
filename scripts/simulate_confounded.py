import numpy as np
import pandas as pd

SEED = 6
C2 = 0.20

# Sample QC and 1000 Genomes panel file
psam = pd.read_csv("data/processed/chr22_sampleqc.psam", sep="\t").rename(columns={"#IID": "IID"})
panel = pd.read_csv(
    "data/raw/integrated_call_samples_v3.20130502.ALL.panel",
    sep="\t", usecols=[0, 1, 2], header=0, names=["IID", "pop", "super_pop"],
)
df = psam[["IID"]].merge(panel, on="IID", how="left")
assert df["pop"].notna().all(), "samples missing from the panel file"

rng = np.random.default_rng(SEED)

# Random environmental shift per population
pops = sorted(df["pop"].unique())
shift = pd.Series(rng.normal(size=len(pops)), index=pops)

# Apply population shift to sampples
env = df["pop"].map(shift).to_numpy()
env = env - env.mean()
env = env * np.sqrt(C2 / env.var())

# Add noise
noise = rng.normal(0, np.sqrt(1 - C2), size=len(df))
y = env + noise

pd.DataFrame({"#IID": df["IID"], "trait_conf": y}).to_csv(
    "data/processed/pheno_confounded.txt", sep="\t", index=False)

# Output key
key = pd.DataFrame({"pop": pops, "shift": shift.to_numpy() * np.sqrt(C2 / df["pop"].map(shift).var())})
key = key.merge(panel[["pop", "super_pop"]].drop_duplicates(), on="pop")
key.sort_values(["super_pop", "pop"]).to_csv("data/processed/confounder_key.csv", index=False)

df["trait"] = y
print(df.groupby("super_pop")["trait"].agg(["count", "mean", "std"]).round(3))
print("overall variance of trait:", round(y.var(), 3))