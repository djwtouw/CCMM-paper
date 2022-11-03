"""
Code to plot a clusterpath and dendrogram for a small data set. LaTeX is
required.

"""

import matplotlib.pyplot as plt
import numpy as np

from src.Python.plot_clusterpath import plot_clusterpath
from src.Python.plot_dendr import load_dendr, plot_dendr


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Computer Modern Sans Serif"]})


# %%

path = "Methodology/Output/Dendrogram/"
merge, height, order = load_dendr(path)

fs = 10

fig = plt.figure(figsize=[2.75, 2.2], dpi=500)
ax = fig.add_subplot(111)
plot_dendr(load_dendr(path),
           leaf_len=height[0], lc="black", lw=0.8, text_offset=-0.50, fs=fs)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticks((0, 1, 2, 3, 4, 5, 6))
ax.set_yticklabels(labels=["$0$", "$1$", "$2$", "$3$", "$4$", "$5$", "$6$"],
                   fontsize=fs)

plt.ylabel(r"$\lambda^{1/2}$", fontsize=fs, rotation=0, labelpad=9)
plt.ylim(0, int(height[-1]) + 0.6)

plt.savefig("Methodology/Figures/Dendrogram.pdf", bbox_inches="tight")


# %%

res = np.genfromtxt(path + "res.csv")
y = np.array([0, 0, 0, 0, 0, 0, 0])


plt.figure(figsize=[3, 2.2], dpi=300)
plot_clusterpath(res, y, s=2.5, lw=0.4, fs_label=10, r_tol=2,
                 labels=["$1$", "$2$", "$3$", "$4$", "$5$", "$6$", "$7$"],
                 lab_x_off=[+0.00, +0.00, +0.00, +0.00, +0.00, +0.02, +0.00],
                 lab_y_off=[+0.08, +0.08, +0.08, +0.08, +0.08, +0.08, +0.08])
plt.ylim([-1.4, 2.05])
plt.savefig("Methodology/Figures/Dendrogram_cp.pdf", bbox_inches="tight")
