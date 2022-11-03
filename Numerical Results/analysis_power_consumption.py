"""
Script for analyzing the individual household power consumption data. The data
set is too large for github, to download it, use the download_hpc() function in
src/R/download_hpc.R.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ccmmpy import SparseWeights, CCMM


def save_value(path, value):
    with open(path, "w") as f:
        f.write(f"{value}")

    return None


# %%

# Load data
df = pd.read_csv("Numerical Results/Data/UCI/power_consumption.csv")
df.TimeOfDay -= 1

# Select timeframe
year1 = np.where((df.Date == "16/12/2007") & (df.Time == "17:23:00"))[0][0]
year2 = np.where((df.Date == "16/12/2008") & (df.Time == "17:23:00"))[0][0]
df = df.iloc[:year2, :]
# df = df.iloc[:(1440 * 60), :]

# Select variables
X = df.iloc[:, 2:9].values

# Standardize the columns
m = X.mean(axis=0)
s = X.std(axis=0)
for j in range(X.shape[1]):
    X[:, j] = (X[:, j] - m[j]) / s[j]

# Settings
k = 15
phi = 0.5

# Compute weights
W = SparseWeights(X, k, phi)

# Look for a range of clusters to find the lambda for which 2 clusters are left
res = CCMM(X, W).convex_clustering(2, 20, verbose=1)

print(res.phase_1_instances)
print(res.phase_2_instances)

# Save time and number of instances solved
path = "Numerical Results/Output/power_consumption_time_1M.csv"
save_value(path, res.elapsed_time)
path = "Numerical Results/Output/power_consumption_ph1_iterations_1M.csv"
save_value(path, res.phase_1_instances)
path = "Numerical Results/Output/power_consumption_ph2_iterations_1M.csv"
save_value(path, res.phase_2_instances)

# Save clusterings
path = "Numerical Results/Output/power_consumption_clusters_1M.csv"

clusters = np.empty((res.info["clusters"].shape[0], X.shape[0]), dtype=int)
for i, c in enumerate(res.info["clusters"].values):
    clusters[i, :] = res.clusters(c)
clusters = clusters[::-1, :]

np.savetxt(path, clusters, delimiter=",", fmt="%i")

# Plot histograms
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
_ = ax1.hist(df.TimeOfDay.values[clusters[0, :] == 0],
             bins=np.arange(0, 1440 + 1e-6, 20), color="orangered")
_ = ax2.hist(df.TimeOfDay.values[clusters[0, :] == 1],
             bins=np.arange(0, 1440 + 1e-6, 20), color="cornflowerblue")
