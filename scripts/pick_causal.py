import numpy as np
import pandas as pd

SEED = 45
N_CAUSAL = 10
MIN_MAF = 0.05

# Read frequency file
freq = pd.read_csv("results/pheno/chr22_pruned_freq.afreq", sep="\t")

# Determine minor allele frequency (MAF)
freq["MAF"] = np.minimum(freq["ALT_FREQS"], 1 - freq["ALT_FREQS"])

# Keep variants that meet the MAF threshold
candidates = freq[freq["MAF"] >= MIN_MAF]
print(f"{len(freq)} pruned variants; {len(candidates)} with MAF >= {MIN_MAF}")

# Keep 10 randomly selected causal variants
picked = candidates.sample(n=N_CAUSAL, random_state=SEED).copy()
picked["POS"] = picked["ID"].str.split(":").str[1].astype(int)
picked = picked.sort_values("POS")

picked["ID"].to_csv("data/processed/causal_ids.txt", index=False, header=False)
print(picked[["ID", "MAF"]].to_string(index=False))