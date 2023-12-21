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

# Colorblind safe colors
colors = [(252, 141, 98),
          (141, 160, 203)]
for c_i in range(len(colors)):
    colors[c_i] = [c / 255 for c in colors[c_i]]

# Figure size
figsize = [2, 1.46]

plt.figure(figsize=figsize, dpi=300)
plot_clusterpath(res1, y, s=2.5, lw=0.4, fs_label=10, r_tol=3, colors=colors)
# Modify y limits a bit
ylim = list(plt.gca().get_ylim())
yrange = ylim[1] - ylim[0]
ylim[0] -= yrange * 0.08
ylim[1] += yrange * 0.08
plt.gca().set_ylim(ylim)
# Continue to saving the image
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Dense_vs_Sparse/Dense_2Moons.pdf", 
            bbox_inches="tight")

plt.figure(figsize=figsize, dpi=300)
plot_clusterpath(res2, y, s=2.5, lw=0.4, fs_label=10, r_tol=3, colors=colors)
# Modify y limits a bit
plt.gca().set_ylim(ylim)
# Continue to saving the image
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Dense_vs_Sparse/Sparse_2Moons.pdf", 
            bbox_inches="tight")

plt.figure(figsize=figsize, dpi=300)
plot_clusterpath(res3, y, s=2.5, lw=0.4, fs_label=10, r_tol=3, alpha=0.7, 
                 colors=colors)
# Modify y limits a bit
plt.gca().set_ylim(ylim)
# Continue to saving the image
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Dense_vs_Sparse/Unit_2Moons.pdf", 
            bbox_inches="tight")
