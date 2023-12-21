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

res1 = np.genfromtxt("Methodology/Output/Connectedness/res_DC.csv")
res2 = np.genfromtxt("Methodology/Output/Connectedness/res_SC.csv")
res3 = np.genfromtxt("Methodology/Output/Connectedness/res_MST.csv")
y = np.genfromtxt("Methodology/Data/Connectedness/y.csv")

# Plot settings
figsize = [3, 2]
s = 3

# Colorblind safe colors
colors = [(102, 194, 165),
          (252, 141, 98),
          (141, 160, 203)]
for c_i in range(len(colors)):
    colors[c_i] = [c / 255 for c in colors[c_i]]

plt.figure(figsize=figsize, dpi=300)
plot_clusterpath(res1, y, s=s, lw=0.4, fs_label=10, r_tol=3, colors=colors)
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Connectedness/DC_Path.pdf",
            bbox_inches="tight")

plt.figure(figsize=figsize, dpi=300)
plot_clusterpath(res2, y, s=s, lw=0.4, fs_label=10, r_tol=3, colors=colors)
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Connectedness/SC_Path.pdf",
            bbox_inches="tight")

plt.figure(figsize=figsize, dpi=300)
plot_clusterpath(res3, y, s=s, lw=0.4, fs_label=10, r_tol=3, colors=colors)
plt.gca().set_aspect("equal")
plt.savefig("Methodology/Figures/Connectedness/MST_Path.pdf",
            bbox_inches="tight")


# %%

X = res1[:75, 2:]
Ws = np.genfromtxt("Methodology/Output/Connectedness/knn_weights_SC.csv")
Wd = np.genfromtxt("Methodology/Output/Connectedness/knn_weights_DC.csv")
Wm = np.genfromtxt("Methodology/Output/Connectedness/knn_weights_MST.csv")

Ws = dense_W(Ws, X.shape[0])
Wd = dense_W(Wd, X.shape[0])
Wm = dense_W(Wm, X.shape[0])

y_cols = []
for i in range(len(y)):
    if y[i] == 0:
        y_cols.append(colors[0])
    if y[i] == 1:
        y_cols.append(colors[1])
    if y[i] == 2:
        y_cols.append(colors[2])

# DISCONNECTED
plt.figure(figsize=figsize, dpi=300)
plt.scatter(X[:, 0], X[:, 1], c=y_cols, zorder=1, s=s)
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
plt.savefig("Methodology/Figures/Connectedness/DC_Graph.pdf",
            bbox_inches="tight")

# SYMMETRIC CIRCULANT
plt.figure(figsize=figsize, dpi=300)
plt.scatter(X[:, 0], X[:, 1], c=y_cols, zorder=1, s=s)
for i in range(Ws.shape[0]):
    for j in range(i):
        if Ws[i, j] > 0 and Wd[i, j] == 0:
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
plt.savefig("Methodology/Figures/Connectedness/SC_Graph.pdf",
            bbox_inches="tight")

# MINIMUM SPANNING TREE
plt.figure(figsize=figsize, dpi=300)
plt.scatter(X[:, 0], X[:, 1], c=y_cols, zorder=1, s=s)
for i in range(Ws.shape[0]):
    for j in range(i):
        if Wm[i, j] > 0 and Wd[i, j] == 0:
            plt.plot([X[i, 0], X[j, 0]], [X[i, 1], X[j, 1]], zorder=-1,
                     color="darkgrey", linewidth=0.4, linestyle="--",
                     alpha=0.5)
            plt.plot([X[i, 0], X[j, 0]], [X[i, 1], X[j, 1]], zorder=-2,
                     color="coral", linewidth=2.0, alpha=0.3)
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
plt.savefig("Methodology/Figures/Connectedness/MST_Graph.pdf",
            bbox_inches="tight")
