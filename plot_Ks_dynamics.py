"""Plots salinity, sodicity, and hydraulic conductivity for different soil types. Only means plotted. Output saved in soil_type_comparision.pdf."""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# soils to plot
soils = ["class_1_rev", "class_1_irr",
         "class_2_rev", "class_2_irr",
         "class_3_rev", "class_3_irr",
         "hydrus"]

# colors
colors = ['#fdbf6f', '#ff7f00',
          "xkcd:green", "xkcd:dark green",
          "xkcd:purple", "xkcd:royal purple",
          "xkcd:dark grey"]

colormaps = [plt.cm.Oranges, plt.cm.Oranges,
             plt.cm.Greens, plt.cm.Greens,
             plt.cm.Purples, plt.cm.Purples,
             plt.cm.Greys]

# labels
labels = ["Weight function A1",
          "Weight function A2",
          "Weight function B1",
          "Weight function B2",
          "Weight function C1",
          "Weight function C2",
          "McNeal function"]

alphas = [1.0, 1.0, 0.3, 0.3, 0.3, 0.3, 1]

# simulation years and rainy days for shading
years = 2
rainy_days = 130

# Begin plotting
plt.ion()
fig = plt.figure(1)
plt.clf()
fig, ax = plt.subplots(num = 1,
                       nrows=7,
                       ncols=1)

fig.subplots_adjust(left=0.11,
                    right=0.97,
                    top=0.98,
                    bottom=0.1,
                    hspace=0.15,
                    wspace=0.15)


# add shading regions and text labeling irrigation and rain
rain_peak = 365 - rainy_days/2
rain_start = rain_peak - rainy_days/2
rain_end = rain_peak + rainy_days/2
for i in range(years):
    rain_peak = 365*(i+1) - rainy_days/2
    rain_start = rain_peak - rainy_days/2
    rain_end = rain_peak + rainy_days/2
    for jj in range(7):
        ax[jj].axvspan(rain_start,
                         rain_end,
                         alpha=0.15,
                         color='xkcd:azure',
                         lw = 0)

# label rainfall irrigation periods
ax[0].annotate('irrigation',
               (25, 0.80),)
ax[0].annotate('rainfall',
               (255, 0.80),)
ax[0].annotate('irrigation',
               (380, 0.80),)
ax[0].annotate('rainfall',
               (620, 0.80),)

# panel a-c
for ii in range(len(soils)):

    if ii % 2 == 0:
        lower = 0.2
        upper = 0.5
    else:
        lower = 0.6
        upper = 0.8

    # stochastic
    input_data_Ks = pd.read_csv("./output_csvs/deg_"+soils[ii]+"_Ks.csv",
                                header = None)
    input_data_Ks.index = input_data_Ks.index/10
    vals = np.linspace(lower,upper,1000)
    np.random.shuffle(vals)
    cmap = plt.cm.colors.ListedColormap(colormaps[ii](vals))
    input_data_Ks.plot(kind='line',
                       ax=ax[ii],
                       colormap=cmap,
                       legend = False,
                       alpha = alphas[ii],
                       lw = 0.2)

    # ensemble mean
    data = pd.read_csv("./output_csvs/deg_"+soils[ii]+".csv")
    ax[ii].plot(data['time'],
                data['hyd_cond'],
                color=colors[ii],
                lw=1.3)

    # label panel
    ax[ii].text(0.01, 0.65,
                labels[ii],
                transform=ax[ii].transAxes,
                horizontalalignment='left',
                verticalalignment='top',
                fontweight='bold')


# tick ranges and axis labels
for ii in range(7):
    ax[ii].set(xlim = [0, 1.5*365],
                ylim = [0.75, 1.05],
                xticks = [],
                ylabel = "Rel. $K_{s}$")

ax[6].set(xlabel = "Time (years)",
            xticks = [0, 365, 730],
            xticklabels = ['0', '1', '2'],
            ylim = [0.0, 1.15])


# save figure
fig.savefig("./figure_Ks.pdf")
