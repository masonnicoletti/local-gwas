import numpy as np
import pandas as pd

new = pd.read_csv("results/pipeline/pca/chr22_pca.eigenvec", sep="\t").set_index("#IID")
old = pd.read_csv("results/pca/chr22_pca.eigenvec", sep="\t").set_index("#IID")

print("same sample IDs:", set(new.index) == set(old.index))
old = old.loc[new.index]

new_ev = np.loadtxt("results/pipeline/pca/chr22_pca.eigenval")
old_ev = np.loadtxt("results/pca/chr22_pca.eigenval")
print("max eigenvalue difference:", f"{np.abs(new_ev - old_ev).max():.2e}")

for pc in new.columns:
    a = new[pc].to_numpy()
    b = old[pc].to_numpy()
    s = np.sign(a @ b)
    label = "flipped" if s < 0 else "same sign"
    print(pc, label, "max difference after alignment:", f"{np.abs(a - s * b).max():.2e}")