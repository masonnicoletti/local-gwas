import numpy as np
import pandas as pd
from scipy.stats import chi2

THRESH = 5e-8
RUNS = {
    "none":               "results/assoc/conf_nocov.trait_conf.glm.linear",
    "PC1-10":             "results/assoc/conf_pcs.trait_conf.glm.linear",
    "PC1-20":             "results/assoc/conf_pcs20.trait_conf.glm.linear",
    "PC1-30":             "results/assoc/conf_pcs30.trait_conf.glm.linear",
    "population dummies": "results/assoc/conf_pop.trait_conf.glm.linear",
}

# Get association statistics for different PC covariate sets
rows = []
for name, path in RUNS.items():
    df = pd.read_csv(path, sep="\t")
    p = df.loc[df["TEST"] == "ADD", "P"].dropna().to_numpy()
    rows.append({
        "covariates": name,
        "lambda": round(float(np.median(chi2.isf(p, df=1)) / chi2.ppf(0.5, df=1)), 3),
        "n_sig": int((p < THRESH).sum()),
        "min_P": f"{p.min():.2e}",
        "share_P<0.05": round(float((p < 0.05).mean()), 3),
    })
print(pd.DataFrame(rows).to_string(index=False))