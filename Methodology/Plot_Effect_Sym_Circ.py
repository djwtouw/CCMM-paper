"""
Code to visualize the effect of using a symmetric circulant for constructing a
convex clustering weight matrix. LaTeX is required.

"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src.Python.plot_clusterpath import plot_clusterpath

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Computer Modern Sans Serif"]})


def dense_W(dok_w, n):
    result = np.zeros([n, n])

    for i in range(dok_w.shape[0]):
        result[int(dok_w[i, 0] - 1), int(dok_w[i, 1] - 1)] = dok_w[i, 2]
        result[int(dok_w[i, 1] - 1), int(dok_w[i, 0] - 1)] = dok_w[i, 2]

    return result


# %%

res1 = np.genfromtxt("Methodology/Output/Sym_Circ/res_nSC.csv")
res2 = np.genfromtxt("Methodology/Output/Sym_Circ/res_SC.csv")
y = np.genfromtxt("Methodology/Data/Sym_Circ/Sym_Circ_y.csv")

plt.figure(figsize=[3, 2], dpi=300)
plot_clusterpath(res1, y, s=2.5, lw=0.4, fs_label=10, r_tol=3)
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Discon_Path.pdf", bbox_inches="tight")

plt.figure(figsize=[3, 2], dpi=300)
plot_clusterpath(res2, y, s=2.5, lw=0.4, fs_label=10, r_tol=3)
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Con_Path.pdf", bbox_inches="tight")


# %%

X = res1[:75, 2:]
Wc = np.genfromtxt("Methodology/Output/Sym_Circ/knn_weights_SC.csv")
Wd = np.genfromtxt("Methodology/Output/Sym_Circ/knn_weights_nSC.csv")

Wc = dense_W(Wc, X.shape[0])
Wd = dense_W(Wd, X.shape[0])

y_cols = []
for i in range(len(y)):
    if y[i] == 0:
        y_cols.append(sns.color_palette("muted")[0])
    if y[i] == 1:
        y_cols.append(sns.color_palette("muted")[1])
    if y[i] == 2:
        y_cols.append(sns.color_palette("muted")[2])

plt.figure(dpi=300, figsize=[3, 2])
plt.scatter(X[:, 0], X[:, 1], c=y_cols, zorder=1, s=3)
for i in range(Wc.shape[0]):
    for j in range(i):
        if Wd[i, j] > 0:
            plt.plot([X[i, 0], X[j, 0]], [X[i, 1], X[j, 1]], zorder=-1,
                     color="black", linewidth=0.4, linestyle="--",
                     alpha=0.8)
plt.tick_params(axis="x", bottom=False, labelbottom=False)
plt.tick_params(axis="y", left=False, labelleft=False)
plt.xlabel("$x_1$", fontsize=10)
plt.ylabel("$x_2$", fontsize=10, rotation=0, labelpad=9)
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Discon_Graph.pdf", bbox_inches="tight")


plt.figure(dpi=300, figsize=[3, 2])
plt.scatter(X[:, 0], X[:, 1], c=y_cols, zorder=1, s=3)
for i in range(Wc.shape[0]):
    for j in range(i):
        if Wc[i, j] > 0 and Wd[i, j] == 0:
            plt.plot([X[i, 0], X[j, 0]], [X[i, 1], X[j, 1]], zorder=-1,
                     color="darkgrey", linewidth=0.4, linestyle="--",
                     alpha=0.5)
for i in range(Wd.shape[0]):
    for j in range(i):
        if Wd[i, j] > 0:
            plt.plot([X[i, 0], X[j, 0]], [X[i, 1], X[j, 1]], zorder=-1,
                     color="black", linewidth=0.4, linestyle="--",
                     alpha=0.8)
plt.tick_params(axis="x", bottom=False, labelbottom=False)
plt.tick_params(axis="y", left=False, labelleft=False)
plt.xlabel("$x_1$", fontsize=10)
plt.ylabel("$x_2$", fontsize=10, rotation=0, labelpad=9)
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Con_Graph.pdf", bbox_inches="tight")
