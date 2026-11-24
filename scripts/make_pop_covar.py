import pandas as pd

# Sample QC and 1000 Genomcs panel
psam = pd.read_csv("data/processed/chr22_sampleqc.psam", sep="\t").rename(columns={"#IID": "IID"})
panel = pd.read_csv("data/raw/integrated_call_samples_v3.20130502.ALL.panel",
                    sep="\t", usecols=[0, 1], header=0, names=["IID", "pop"])
df = psam[["IID"]].merge(panel, on="IID", how="left")

# Encode column per population
dummies = pd.get_dummies(df["pop"], prefix="pop", drop_first=True).astype(int)
out = pd.concat([df[["IID"]].rename(columns={"IID": "#IID"}), dummies], axis=1)
out.to_csv("data/processed/covar_pop.txt", sep="\t", index=False)
print(out.shape)