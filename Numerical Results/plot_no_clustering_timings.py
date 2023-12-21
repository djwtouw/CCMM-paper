import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from itertools import permutations
from matplotlib import rc
from scipy import stats

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Computer Modern Sans Serif"]})

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


def SpearmanR(x, y, n=10000):
    # Save Spearman rank correlation coefficients
    spcor = np.empty(n)
    for i in range(n):
        xi = np.random.permutation(x)
        spcor[i] = stats.spearmanr(y, xi).correlation

    # Spearman rank correlation coefficient of the sample
    correlation = stats.spearmanr(y, x).correlation

    # Compute p-value
    p_value = sum(abs(np.array(spcor)) >= correlation) / n

    return correlation, p_value


# %%

timings_ssnal = pd.read_csv(
    "Numerical Results/Output/no_clustering_timings_ccmm_ssnal.csv"
)

# Check how many were completed and select
idx = timings_ssnal.iloc[:, 0] > 0
timings_ssnal = timings_ssnal[idx]

# The values for the x axis
N_ssnal = np.arange(500, sum(idx) * 500 + 1, 500)

# Ratio
ratio_ssnal = timings_ssnal["SSNAL"] / timings_ssnal["CCMM"]

# Color palette
colors = sns.color_palette("muted", 3)

plt.figure(dpi=400, figsize=[6, 2])
plt.subplot(1, 2, 1)
plt.scatter(N_ssnal, ratio_ssnal, color=colors[0], s=5)
plt.xlabel("Number of objects $n$", fontsize=10)
plt.ylabel("Speedup", fontsize=10)
plt.xticks([0, 5000, 10000, 15000, 20000],
           ["0", "5,000", "10,000", "15,000", "20,000"],
           fontsize=7)
plt.yticks([50, 60, 70, 80, 90, 100, 110],
           fontsize=7)
plt.ylim([45, 115])
plt.xlim([-800, 20800])

plt.savefig("Numerical Results/Figures/" +
            "timings_ratio_ssnal.pdf", bbox_inches="tight")

# Compute the Spearman rank correlation coefficient
test = SpearmanR(ratio_ssnal, N_ssnal)
print("Spearman's Rank Correlation")
print("correlation: {:.4f}".format(test[0]))
print("p-value:     {:.4f}".format(test[1]))

# Load the values of the loss function
losses_ssnal = pd.read_csv(
    "Numerical Results/Output/no_clustering_losses_ccmm_ssnal.csv"
)

# Find the maximum deviation from the loss attained by SSNAL in the wrong
# direction: where CCMM attains a higher loss
print(
    "Maximum relative deviation in the loss: {:.7f}"
    .format(max(losses_ssnal["CCMM"] / losses_ssnal["SSNAL"]) - 1)
)


# %%

timings_ama = pd.read_csv(
    "Numerical Results/Output/no_clustering_timings_ccmm_ama.csv"
)

# Check how many were completed and select
idx = timings_ama.iloc[:, 0] > 0
timings_ama = timings_ama[idx]

# The values for the x axis
N_ama = np.arange(500, sum(idx) * 500 + 1, 500)

# Ratio
ratio_ama = timings_ama["AMA"] / timings_ama["CCMM"]

plt.figure(dpi=400, figsize=[6, 2])
plt.subplot(1, 2, 1)
plt.scatter(N_ama, ratio_ama, color=colors[0], s=5)
plt.xlabel("Number of objects $n$", fontsize=10)
plt.ylabel("Speedup", fontsize=10)
plt.xticks([0, 2500, 5000, 7500, 10000],
           ["0", "2,500", "5,000", "7,500", "10,000"],
           fontsize=7)
plt.yticks([0, 100, 200, 300, 400, 500, 600],
           fontsize=7)
plt.ylim([-50, 650])
plt.xlim([-400, 10400])

plt.savefig("Numerical Results/Figures/" +
            "timings_ratio_ama.pdf", bbox_inches="tight")

# Compute the Spearman rank correlation coefficient
test = SpearmanR(ratio_ama, N_ama)
print("Spearman's Rank Correlation")
print("correlation: {:.4f}".format(test[0]))
print("p-value:     {:.4f}".format(test[1]))

# Load the values of the loss function
losses_ama = pd.read_csv(
    "Numerical Results/Output/no_clustering_losses_ccmm_ama.csv"
)

# Find the maximum deviation from the loss attained by AMA in the wrong
# direction: where CCMM attains a higher loss
print(
    "Maximum relative deviation in the loss: {:.7f}"
    .format(max(losses_ama["CCMM"] / losses_ama["AMA"]) - 1)
)
