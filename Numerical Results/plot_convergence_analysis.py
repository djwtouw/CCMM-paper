import struct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.lines import Line2D


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Computer Modern Sans Serif"]})


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


def load_vector(file_path):
    # Open the binary file for reading in binary mode
    with open(file_path, "rb") as file:
        # Read the binary data from the file
        binary_data = file.read()

    # Define the format of each double in the binary data
    # 'd' represents a double-precision floating-point number
    double_format = 'd'

    # Calculate the number of doubles in the binary data
    num_doubles = len(binary_data) // struct.calcsize(double_format)

    # Use struct.unpack to convert the binary data to a list of doubles
    doubles = struct.unpack(f'{num_doubles}{double_format}', binary_data)

    # Now, 'doubles' contains the vector of doubles from the binary file
    return np.array(doubles)


# %%

# Load orders
orders = load_vector("Numerical Results/Output/convergence_analysis.bin")

# Range for the x axis histogram
hist_min = 0.0
hist_max = 2.0

# Set the binsize and create a vector with the start and endpoints of each bin
binsize = 0.01
bins = np.arange(hist_min, hist_max + binsize / 2, binsize)

# Initialize a vector with counts for each bin
counts = np.zeros(bins.size - 1, dtype=int)

# Count the occurence of each time of the day in the cluster
for o in orders:
    if o >= hist_min and o <= hist_max:
        counts[int((o - hist_min) // binsize)] += 1

# Exact histogram
x = []
y = []
b_max = bins.size
for i in range(0, b_max):
    x.append(bins[i])
    x.append(bins[i])

    if i == 0:
        y.append(0)
        y.append(counts[i])
    elif i == b_max - 1:
        y.append(counts[i - 1])
        y.append(0)
    else:
        y.append(counts[i - 1])
        y.append(counts[i])

# Colorblind safe colors
colors = [(102, 194, 165),
          (252, 141, 98),
          (141, 160, 203)]
for c_i in range(len(colors)):
    colors[c_i] = [c / 255 for c in colors[c_i]]

# Plot histogram
plt.figure(figsize=[2.5, 2], dpi=500)
plt.fill(x, y, linewidth=0, color=colors[2])
plt.xticks(np.arange(hist_min, hist_max + 1e-10, (hist_max - hist_min) / 8),
           fontsize=7)
plt.yticks([5000, 10000, 15000, 20000, 25000],
           ["5,000", "10,000", "15,000", "20,000", "25,000"], fontsize=7)
plt.ylabel("Frequency", fontsize=10)
plt.xlabel("Estimated order of convergence $\\hat{q}$", fontsize=10)
plt.ylim([0, 30000])
plt.xlim([hist_min, hist_max])
plt.plot([orders.mean(), orders.mean()], plt.gca().get_ylim(), lw=0.75,
         color=colors[0], ls="--", label="Mean")
plt.plot([np.median(orders), np.median(orders)], plt.gca().get_ylim(), lw=0.75,
         color=colors[1], ls="--", label="Median")
plt.legend(loc="upper right", fontsize=7, frameon=False)

plt.savefig("Numerical Results/Figures/convergence_1.pdf",
            bbox_inches="tight")


# %%


# Load data
data = pd.read_csv("Numerical Results/Output/convergence_analysis_losses.csv",
                   header=None)
losses = data.values

# Colorblind safe colors
colors = [(102, 194, 165),
          (252, 141, 98),
          (141, 160, 203)]
for c_i in range(len(colors)):
    colors[c_i] = [c / 255 for c in colors[c_i]]

# Plot losses
plt.figure(figsize=[2.5, 2], dpi=500)
for i in range(losses.shape[1]):
    plt.plot(losses[:, i], alpha=0.10, color=colors[1], lw=0.25)
plt.plot(losses.mean(axis=1), lw=0.5, c="black")
plt.xticks([1, 26, 51, 76, 101], [0, 25, 50, 75, 100], fontsize=7)
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=7)
plt.ylabel("Scaled loss", fontsize=10)
plt.xlabel("Iteration $t$", fontsize=10)

# Get content of legend
handles, labels = plt.gca().get_legend_handles_labels()

# Create manual entries
line1 = Line2D([0], [0], label="Loss", color=colors[1], lw=0.5)
line2 = Line2D([0], [0], label="Mean loss", color="black", lw=0.5)

# Add manual symbols to auto legend
handles.extend([line1, line2])
plt.legend(handles=handles, loc="upper right", fontsize=7, frameon=False)

plt.savefig("Numerical Results/Figures/convergence_2.pdf",
            bbox_inches="tight")
