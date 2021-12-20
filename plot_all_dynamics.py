"""Plots salinity, sodicity, and hydraulic conductivity for different soil types. Only means plotted. Output saved in soil_type_comparision.pdf."""

import matplotlib.pyplot as plt
import pandas as pd

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

# line styles
styles = ['-', '--',
          '-', '--',
          '-', '--',
          '-']

# simulation years and days needed for shading
years = 1.5
rainy_days = 130

# plot settings

# Begin plotting
plt.ion()
fig = plt.figure(1)
plt.clf()
fig, ax = plt.subplots(num = 1,
                       nrows=4,
                       ncols=1,
                       sharex=True)
fig.subplots_adjust(left=0.1,
                    right=0.98,
                    top=0.98,
                    bottom=0.18,
                    hspace=0.2,
                    wspace=0.15)

# add shading regions and text labeling irrigation and rain
rain_peak = 365 - rainy_days/2
rain_start = rain_peak - rainy_days/2
rain_end = rain_peak + rainy_days/2

for i in range(int(years)):
    rain_peak = 365*(i+1) - rainy_days/2
    rain_start = rain_peak - rainy_days/2
    rain_end = rain_peak + rainy_days/2
    ax[0].axvspan(rain_start, rain_end, alpha=0.15, color='xkcd:azure', lw = 0)
    ax[1].axvspan(rain_start, rain_end, alpha=0.15, color='xkcd:azure', lw = 0)
    ax[2].axvspan(rain_start, rain_end, alpha=0.15, color='xkcd:azure', lw = 0)
    ax[3].axvspan(rain_start, rain_end, alpha=0.15, color='xkcd:azure', lw = 0)

leg_lines = []
for ii in range(len(soils)):
    # import data
    data = pd.read_csv("./output_csvs/deg_"+soils[ii]+".csv",
                       index_col = 0)
    lin = ax[0].plot(data['time'],
                     data['water'],
                     color=colors[ii],
                     lw=1,
                     ls=styles[ii])
    ax[1].plot(data['time'],
               data['salinity'],
               color=colors[ii],
               lw=1,
               ls=styles[ii])
    ax[2].plot(data['time'],
               data['sodicity'],
               color=colors[ii],
               lw=1,
               ls=styles[ii])
    ax[3].plot(data['time'],
               data['hyd_cond'],
               color=colors[ii],
               lw=1,
               ls=styles[ii])

    leg_lines += [lin]

# tick ranges and labels
ax[0].set_xlim([0,years*365])
ax[1].set_xlim([0,years*365])
ax[2].set_xlim([0,years*365])
ax[3].set_xlim([0,years*365])
#ax[3].set_ylim([0.9,1.03])
ax[3].set_xticks([0, 365])
ax[3].set_xticklabels(['0', '1'])


# axis labels
fig.text(0.5, 0.1, "Time (years)", ha='center')
ax[0].set_ylabel("$s$")
ax[0].text(0.01, 0.22, "Water content", transform=ax[0].transAxes,
         horizontalalignment='left', verticalalignment='top')
ax[1].set_ylabel("$C_s$")
ax[1].text(0.01, 0.93, "Salinity", transform=ax[1].transAxes,
         horizontalalignment='left', verticalalignment='top')
ax[2].set_ylabel("$E_x$")
ax[2].text(0.01, 0.93, "Sodicity", transform=ax[2].transAxes,
         horizontalalignment='left', verticalalignment='top')
ax[3].set_ylabel("Rel. $K_{s}$")
ax[3].text(0.01, 0.6, "Relative $K_{s}$", transform=ax[3].transAxes,
         horizontalalignment='left', verticalalignment='top')

# legend
leg = fig.legend([leg_lines[0][0], leg_lines[1][0], leg_lines[2][0],
                 leg_lines[3][0], leg_lines[4][0], leg_lines[5][0],
                 leg_lines[6][0]],
                  ["A1 weights", "A2 weights", "B1 weights", "B2 weights",
                   "C1 weights", "C2 weights", "McNeal"],
                  loc='upper center', bbox_to_anchor=(0.5, 0.08), ncol=4,
                  labelspacing = 0)
for line in leg.get_lines():
    line.set_linewidth(2.0)

fig.savefig("./figure_all_dynamics.pdf")
