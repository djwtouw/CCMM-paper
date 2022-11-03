"""
Code to visualize the differences between clusterspaths for different types of
weights. LaTeX is required.

"""

import numpy as np
import matplotlib.pyplot as plt
from src.Python.plot_clusterpath import plot_clusterpath

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Computer Modern Sans Serif"]})


# %%

res1 = np.genfromtxt("Methodology/Output/Dense_vs_Sparse/res_dense.csv")
res2 = np.genfromtxt("Methodology/Output/Dense_vs_Sparse/res_sparse.csv")
res3 = np.genfromtxt("Methodology/Output/Dense_vs_Sparse/res_unit.csv")
y = np.genfromtxt("Methodology/Data/Dense_vs_Sparse/2half_moons_y.csv")


plt.figure(figsize=[2, 1.46], dpi=300)
plot_clusterpath(res1, y, s=2.5, lw=0.4, fs_label=10, r_tol=3)
plt.savefig("Methodology/Figures/Dense_2Moons.pdf", bbox_inches="tight")

plt.figure(figsize=[2, 1.46], dpi=300)
plot_clusterpath(res2, y, s=2.5, lw=0.4, fs_label=10, r_tol=3)
plt.savefig("Methodology/Figures/Sparse_2Moons.pdf", bbox_inches="tight")

plt.figure(figsize=[2, 1.46], dpi=300)
plot_clusterpath(res3, y, s=2.5, lw=0.4, fs_label=10, r_tol=3, alpha=0.7)
plt.savefig("Methodology/Figures/Unit_2Moons.pdf", bbox_inches="tight")
