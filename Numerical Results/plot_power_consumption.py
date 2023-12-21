"""
Code to create images for the results of using convex clustering on the
individual household electric power consumption data set. LaTeX is required.

"""

import PyPDF2
from matplotlib.patches import Patch
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Computer Modern Sans Serif"]})


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


def pc_hist(data, clusters, cluster, smooth, smoothness, color):
    # Get the times of the day (in minutes) and make sure they start at 0
    times = data["TimeOfDay"][clusters == cluster].values
    times -= times.min()

    # Set the binsize and create a vector with the start and endpoints of each bin
    binsize = 15
    bins = np.arange(0, 1440, binsize)
    if bins[-1] != 1440:
        bins = np.concatenate([bins, np.array([1440])])

    # Initialize a vector with counts for each bin
    counts = np.zeros(bins.size - 1, dtype=int)

    # Count the occurence of each time of the day in the cluster
    for t in times:
        counts[t // binsize] += 1

    # Prepare the count vector for smoothing
    counts_smooth = np.empty(2 * counts.size + 1)
    bins_smooth = np.empty(2 * bins.size - 1)

    # Smoothed histogram
    for i in range(bins_smooth.size):
        if i % 2 == 0:
            bins_smooth[i] = bins[i >> 1]
        else:
            bins_smooth[i] = 0.5 * (bins[i >> 1] + bins[(i >> 1) + 1])

    for i in range(counts_smooth.size):
        if i == 0:
            counts_smooth[i] = counts[i]
        elif i == counts_smooth.size - 1:
            counts_smooth[i] = counts[-1]
        elif i % 2 == 1:
            counts_smooth[i] = counts[i >> 1]
        else:
            counts_smooth[i] = 0.5 * (counts[i >> 1] + counts[(i >> 1) - 1])

    for _ in range(smoothness):
        for i in range(1, counts_smooth.size - 1):
            counts_smooth[i] = 0.25 * counts_smooth[i - 1]\
                + 0.50 * counts_smooth[i]\
                + 0.25 * counts_smooth[i + 1]

    poly_x = np.concatenate((np.array([0]),
                             bins_smooth,
                             np.array([bins_smooth[-1]])))
    poly_y = np.concatenate((np.array([0]),
                             counts_smooth,
                             np.array([0])))

    # Exact histogram
    x = []
    y = []
    b_max = bins.size
    for i in range(0, b_max):
        if i == 0:
            x.append(0)
            x.append(bins[i])
            y.append(0)
            y.append(counts[i])
        elif i == b_max - 1:
            x.append(bins[i])
            x.append(bins[i])
            y.append(counts[i - 1])
            y.append(0)
        else:
            x.append(bins[i])
            x.append(bins[i])
            y.append(counts[i - 1])
            y.append(counts[i])

    if not smooth:
        plt.fill(x, y, color=color, linewidth=0)
    else:
        plt.fill(poly_x, poly_y, color=color, linewidth=0, alpha=0.25)
        plt.plot(bins_smooth, counts_smooth, color="black", linewidth=0.5)

    return None


def plot_dendr(dendr, leaf_len, lc, lw, text_offset, fs, draw_labels=True,
               leaf_points=None, linecols=None):
    merge, height, order = dendr

    leaf_points_map = {}
    height_map = {}

    for lp, o in zip(leaf_points, order):
        leaf_points_map[-(o + 1)] = lp

    for o in order:
        height_map[-(o + 1)] = 0

    for i in range(19):
        p0 = merge[i, 0]
        p1 = merge[i, 1]

        # Horizontal line
        if i not in [6, 7, 14, 15]:
            x = [leaf_points_map[p0], leaf_points_map[p1]]
            y = [height[i], height[i]]
            plt.plot(x, y, c=linecols[i], lw=lw)
        elif i == 7:
            x = [leaf_points_map[4], leaf_points_map[p1]]
            y = [height[i], height[i]]
            plt.plot(x, y, c=linecols[i], lw=lw)
        elif i == 15:
            x = [leaf_points_map[14], leaf_points_map[p1]]
            y = [height[i], height[i]]
            plt.plot(x, y, c=linecols[i], lw=lw)

        leaf_points_map[i + 1] = sum(x) / len(x)

        # Vertical line 1
        x = [leaf_points_map[p0]] * 2
        y = [height_map[p0], height[i]]
        if i != 7:
            plt.plot(x, y, c=linecols[i], zorder=-i, lw=lw)

        # Vertical line 2
        x = [leaf_points_map[p1]] * 2
        y = [height_map[p1], height[i]]
        plt.plot(x, y, c=linecols[i], zorder=-i, lw=lw)

        height_map[i + 1] = height[i]

    return None


# Colorblind safe colors
colors = [(102, 194, 165),
          (252, 141, 98),
          (141, 160, 203)]
for c_i in range(len(colors)):
    colors[c_i] = [c / 255 for c in colors[c_i]]


# %% Load data

data = pd.read_csv("Numerical Results/Data/UCI/power_consumption.csv")
clusters = np.genfromtxt("Numerical Results/Output/power_consumption_" +
                         "clusters_1M.csv", dtype=int, delimiter=",")
data = data.iloc[:clusters.shape[1], :]


# %% Histogram for clusters

plt.figure(figsize=[2.5, 2], dpi=500)
pc_hist(data, clusters[0, :], 1, True, 5, colors[1])
pc_hist(data, clusters[0, :], 0, True, 5, colors[2])

plt.xticks([0, 480, 960, 1440], ["00:00", "08:00", "16:00", "24:00"],
           fontsize=7)
plt.yticks([2500, 5000, 7500, 10000], ["2,500", "5,000", "7,500", "10,000"],
           fontsize=7)
plt.ylabel("Frequency", fontsize=10)
plt.xlabel("Time of day", fontsize=10)
plt.ylim([0, 12500])
plt.xlim([0, 1440])


legend_elements = [Patch(facecolor=(colors[2][0], colors[2][1], colors[2][2],
                                    0.25), edgecolor="black",
                         label="Cluster A", lw=0.5),
                   Patch(facecolor=(colors[1][0], colors[1][1], colors[1][2],
                                    0.25), edgecolor="black",
                         label="Cluster B", lw=0.5)]
plt.legend(handles=legend_elements, loc=2, fontsize=7, frameon=False, ncol=2)

plt.savefig("Numerical Results/Figures/power_consumption_12.pdf",
            bbox_inches="tight")

# %% Preparation for heatmap

clusters = clusters[::-1, :]
clusters += 1

U = np.zeros([clusters.shape[0] + 1, 20], dtype=int)
U[0, :] = np.arange(1, 21)
U[-1, :] = 1

k = 0
clusters20 = []
for i in range(clusters[k].max()):
    indices = np.argwhere(clusters[k, :] == i + 1)
    clusters20.append(set(indices.flatten()))

for k in range(U.shape[0] - 2):
    clusters19 = []
    for i in range(clusters[k + 1].max()):
        indices = np.argwhere(clusters[k + 1, :] == i + 1)
        clusters19.append(set(indices.flatten()))

    for clust19 in clusters19:
        indices = []

        for i, clust20 in enumerate(clusters20):
            if clust20.issubset(clust19):
                indices.append(i)

        if U[k + 1, indices[0]] == 0:
            for idx in indices:
                U[k + 1, idx] = U[k, indices[0]]

# Get the merge table and height vector for the dendrogram
merge = np.zeros([U.shape[1] - 1, 3], dtype=int)
merge[:, 2] = np.arange(1, U.shape[1])
height = np.zeros(U.shape[1] - 1)

cluster_id = -np.arange(1, U.shape[1] + 1, dtype=int)
idx = 0
for k in range(1, U.shape[0]):
    indices = set()
    for i in range(20):
        if U[k, i] != U[k - 1, i]:
            temp = (U[k, i], U[k - 1, i])
            indices = indices.union({temp})

    for ind in indices:
        merge[idx, 0] = cluster_id[ind[0] - 1]
        merge[idx, 1] = cluster_id[ind[1] - 1]

        height[idx] = k

        cluster_id[ind[0] - 1] = idx + 1
        cluster_id[ind[1] - 1] = idx + 1

        idx += 1

merge = merge[:, :2]

print(merge)
print(height)

# Get order for the dendrogram
d = dict()
for i in range(merge.shape[0] + 1):
    d[-(i + 1)] = [i]

for i, entry in enumerate(merge):
    d[i + 1] = d.pop(entry[0]) + d.pop(entry[1])

order = d[merge.shape[0]]
order = np.array(order)

print(order)

del(clust19, clust20, cluster_id, clusters19, clusters20, i, idx, ind,
    indices, k, temp, U)

# Change order for a better looking heatmap
order[:3] = order[:3][::-1]
order[14:16] = order[14:16][::-1]


# %% Plotting heatmap

# Size of colored slabs
counts = np.zeros(20)
for j in range(20):
    counts[j] = (clusters[0, :] == order[j] + 1).sum()
counts = np.log(counts + 5)
counts = counts / counts.sum() * 100


leaf_points = np.zeros(20)
for j in range(1, 20):
    leaf_points[j] = leaf_points[j - 1] + 0.5 * (counts[j - 1] + counts[j])

# Determining the line colors
linecols = [""] * 19


def setcol(index, color, merge, linecols):
    linecols[index] = color

    if merge[index, 0] > 0:
        setcol(merge[index, 0] - 1, color, merge, linecols)
    if merge[index, 1] > 0:
        setcol(merge[index, 1] - 1, color, merge, linecols)


setcol(18, "black", merge, linecols)
setcol(17, colors[1], merge, linecols)
setcol(11, colors[2], merge, linecols)

fs = 10
fig = plt.figure(figsize=[2.06, 2.95], dpi=500)
ax = fig.add_subplot(111)
plot_dendr((merge[:, :], height * 0.75, order), leaf_len=height[0], lc="black",
           lw=0.5, text_offset=-0.40, fs=fs, draw_labels=False,
           leaf_points=leaf_points - 0.5, linecols=linecols)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticks([])

color_grid = np.zeros([7, 20])
data_scaled = data.iloc[:, 2:9]
for j in range(7):
    data_scaled.iloc[:, j] -= data_scaled.iloc[:, j].mean()
    data_scaled.iloc[:, j] /= data_scaled.iloc[:, j].var()**0.5

for i, c in enumerate(order):
    for j in range(7):
        color_grid[j, i] = data_scaled.iloc[clusters[0, :] == c + 1, j].mean()

color_grid

for i in range(7):
    color_grid[i, :] -= color_grid[i, :].min()
    color_grid[i, :] /= color_grid[i, :].max()

cmap = sns.color_palette("Greens", 120 + 1)[10:-10]
width = 8

for i in range(7):
    for j in range(20):
        x1 = leaf_points[j] - counts[j] * 0.5
        x2 = leaf_points[j] + counts[j] * 0.5
        x = [x1, x2, x2, x1]
        y = [-i * width - 0.5,
             -i * width - 0.5,
             -(i + 1) * width - 0.5,
             -(i + 1) * width - 0.5]
        plt.fill(x, y, color=cmap[round(color_grid[i, j] * 100)])

colnames = ["GAP", "GRP", "V", "GI", r"ESM\textsubscript{1}",
            r"ESM\textsubscript{2}", r"ESM\textsubscript{3}"]
for i in range(7):
    plt.text(-0.80 * counts[0], -(i + 0.5) * width - 1, colnames[i],
             rotation=270, ha="center", va="center",  fontsize=6)
plt.subplots_adjust(left=0.15, right=1.0, top=1.0, bottom=0.05)

plt.savefig("Numerical Results/Figures/power_consumption_heatmap.pdf")

# Rotate image 90 degrees
pdf_in = open("Numerical Results/Figures/power_consumption_heatmap.pdf", 'rb')
pdf_reader = PyPDF2.PdfFileReader(pdf_in)
pdf_writer = PyPDF2.PdfFileWriter()

for pagenum in range(pdf_reader.numPages):
    page = pdf_reader.getPage(pagenum)
    page.rotateClockwise(270)
    pdf_writer.addPage(page)

pdf_out = open("Numerical Results/Figures/power_consumption_heatmap_r.pdf",
               'wb')
pdf_writer.write(pdf_out)
pdf_out.close()
pdf_in.close()


# %% Descriptive statistics of the data and clusters

df0 = data.iloc[:, 2:9]
df1 = data.iloc[clusters[-1, :] == 2, 2:9]
df2 = data.iloc[clusters[-1, :] == 1, 2:9]

m0 = df0.mean()
m1 = df1.mean()
m2 = df2.mean()

table = ""
for i, row in enumerate(m0.index):
    table += "{} \t && ".format(row)
    table += "{:.2f} \t & ".format(m0[i])
    table += "{:.2f} \t & ".format(m2[i])
    table += "{:.2f} \\\\ \n".format(m1[i])
table += "\midrule\n"
table += "Number of measurements && "
table += "{:,} \t & ".format(df0.shape[0])
table += "{:,} \t & ".format(df2.shape[0])
table += "{:,} \\\\ \n ".format(df1.shape[0])

print(table)
