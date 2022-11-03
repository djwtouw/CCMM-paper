"""
Code to visualize the results of several analyses performed using CCMM, AMA, 
and SSNAL to perform convex clustering. LaTeX is required.

"""

from matplotlib import rc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Computer Modern Sans Serif"]})

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


# %% Scaling in n for CCMM

timings = pd.read_csv("Numerical Results/Output/scaling_n_ccmm.csv")

labels = ["CCMM", "AMA", "SSNAL"]
colors = sns.color_palette("muted", 3)
lw = 0.75

plt.figure(dpi=400, figsize=[6, 2])
plt.subplot(1, 2, 1)

plt.plot(timings["n"], timings["time"] / timings["iterations"] / 5,
         label=labels[0], color=colors[0], lw=lw)
plt.xlabel("Number of objects $n$", fontsize=10)
plt.ylabel("Time (s)", fontsize=10)
plt.xticks([6000, 9500, 13000, 16500, 20000],
           ["6,000", "9,500", "13,000", "16,500", "20,000"],
           fontsize=7)
plt.yticks([0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040],
           fontsize=7)
plt.ylim([0.0007, 0.0043])

plt.savefig("Numerical Results/Figures/" +
            "scaling_in_n.pdf", bbox_inches="tight")


# %% CPU times for CCMM, AMA, and SSNAL for n=1000..5000 without SC

timings = pd.read_csv("Numerical Results/Output/n1000_to_5000_timings.csv")
timings = timings.drop("completed", axis=1) / 10
timings_log = timings.applymap(np.log)

labels = ["CCMM", "AMA", "SSNAL"]
colors = sns.color_palette("muted", 3)
lw = 0.75

plt.figure(dpi=400, figsize=[6, 2])
plt.subplot(1, 2, 1)

plt.plot(timings_log["CCMM"], label=labels[0], linewidth=lw,
         c=colors[0])
plt.plot(timings_log["AMA"], label=labels[1], linewidth=lw,
         c=colors[1])
plt.plot(timings_log["SSNAL"], label=labels[2], linewidth=lw,
         c=colors[2])

plt.xlabel("Number of objects $n$", fontsize=10)
plt.ylabel("Time (s)", fontsize=10)
plt.xticks([0, 1, 2, 3, 4], ["1,000", "2,000", "3,000", "4,000", "5,000"],
           fontsize=7)
plt.yticks([np.log(0.1), np.log(1), np.log(10), np.log(100), np.log(1000),
            np.log(10000), np.log(100000)],
           ["$10^{-1}$", "$10^{0}$", "$10^{1}$", "$10^{2}$", "$10^{3}$",
            "$10^{4}$", "$10^{5}$"],
           ha="left", fontsize=7)
plt.tick_params(axis='y', pad=17)
plt.ylim([np.log(0.1) - (np.log(100) - np.log(10)) * 0.6,
          np.log(1e5) + (np.log(100) - np.log(10)) * 0.6])
plt.legend(loc=2, fontsize=7, frameon=False, ncol=2)


plt.savefig("Numerical Results/Figures/n1000_to_5000_timings.pdf",
            bbox_inches="tight")


# %% Relative losses w.r.t. SSNAL for CCMM and AMA for n=1000

losses = pd.read_csv("Numerical Results/Output/n1000_losses.csv").drop(0)
ratios = losses.copy()
ratios["CCMM"] /= losses["SSNAL"]
ratios["AMA"] /= losses["SSNAL"]
ratios["SSNAL"] /= losses["SSNAL"]

lambdas = np.arange(0.2, 110.1, 0.2)
M = 120

plt.figure(dpi=400, figsize=[6, 2])
plt.subplot(1, 2, 1)

plt.plot(lambdas[lambdas <= M], ratios["CCMM"][lambdas <= M],
         label=labels[0], linewidth=lw, c=colors[0])
plt.plot(lambdas[lambdas <= M], ratios["AMA"][lambdas <= M],
         label=labels[1], linewidth=lw, c=colors[1])
plt.plot(lambdas[lambdas <= M], ratios["SSNAL"][lambdas <= M],
         label=labels[2], linewidth=lw, c=colors[2])


plt.xticks(fontsize=7)
plt.yticks([1.0000, 1.00005, 1.00010, 1.00015, 1.00020, 1.00025],
           ["1.00000", "1.00005", "1.00010", "1.00015", "1.00020", "1.00025"],
           fontsize=7)
plt.ylim([0.999985, 1.000275])

plt.legend(loc=2, fontsize=7, frameon=False, ncol=2)

plt.xlabel("Regularization parameter $\lambda$", fontsize=10)
plt.ylabel("Relative loss", fontsize=10)


plt.savefig("Numerical Results/Figures/n1000_losses.pdf", bbox_inches="tight")


# %% Relative losses w.r.t. SSNAL for CCMM and AMA for n=5000

losses = pd.read_csv("Numerical Results/Output/n5000_losses.csv").drop(0)
ratios = losses.copy()
ratios["CCMM"] /= losses["SSNAL"]
ratios["AMA"] /= losses["SSNAL"]
ratios["SSNAL"] /= losses["SSNAL"]

lambdas = np.arange(0.2, 110.1, 0.2)

plt.figure(dpi=400, figsize=[6, 2])
plt.subplot(1, 2, 1)

plt.plot(lambdas[lambdas <= M], ratios["CCMM"][lambdas <= M],
         label=labels[0], linewidth=lw, c=colors[0])
plt.plot(lambdas[lambdas <= M], ratios["AMA"][lambdas <= M],
         label=labels[1], linewidth=lw, c=colors[1])
plt.plot(lambdas[lambdas <= M], ratios["SSNAL"][lambdas <= M],
         label=labels[2], linewidth=lw, c=colors[2])

plt.xticks(fontsize=7)
plt.yticks([1.0000, 1.0008, 1.0016, 1.0024, 1.0032, 1.0040],
           ["1.00000", "1.00080", "1.00160", "1.00240", "1.00320", "1.00400"],
           fontsize=7)
plt.ylim([1 - 0.0044 * 0.054545, 1.0044])

plt.legend(loc=2, fontsize=7, frameon=False, ncol=2)

plt.xlabel("Regularization parameter $\lambda$", fontsize=10)
plt.ylabel("Relative loss", fontsize=10)


plt.savefig("Numerical Results/Figures/n5000_losses.pdf", bbox_inches="tight")
