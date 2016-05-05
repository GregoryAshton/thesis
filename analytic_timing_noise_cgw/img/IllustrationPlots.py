#!/usr/bin/python

import matplotlib.pyplot as plt
import TNtools as TN
import numpy as np
import sys
from matplotlib.patches import Rectangle

TN.PlotDefaults()

# Set up the axes and data
ax1 = plt.subplot(111)
N = 9
time = np.arange(N)
dt = time[1] - time[0]
kappa = 0.0
#F1hat = np.random.normal(0, 0.5, N); print F1hat
F1hat = [-0.06673121, -0.16184903, 0.07196745, 0.29913356,
          0.43034329, -0.3241269,  -0.48430681, -0.1279246, 0.0132123]

# Add the 0 line and f1dot
ax1.axhline(0, color="k")
DeltaF1 = np.cumsum(F1hat)
p0, = ax1.step(time, DeltaF1, lw=2, color="k", where="post",
               )

# Add the reference times
#for i in range(len(time)):
#    ax1.plot([time[i]+kappa*dt, time[i]+kappa*dt,], [0, DeltaF1[i]],
#            ls="--", color="k", lw=0.8)
#p1, = ax1.plot(time[:-1]+kappa*dt, DeltaF1[:-1], "ro")
#
#

# Fill up integration
def custom_fill(ax, x, y, dx):
    ax.fill_between([x, x+dx],
                     y1=0, y2=y,
                     alpha=.5, color="g")

end = N-2
for i in range(end):
    custom_fill(ax1, time[i], DeltaF1[i], dt)

custom_fill(ax1, time[end], DeltaF1[end], kappa*dt)



# Add annotations
arr_off = 0.03
text_off = 0.05

sign = 1.0
plt.annotate('', xy=(1, sign*arr_off), xycoords = 'data',
             xytext = (2, sign*arr_off), textcoords = 'data',
             arrowprops = {'arrowstyle':'<->'})
plt.text(1.5, sign*text_off, "$\Delta T$", horizontalalignment='center')

#sign = 1.0
#text_off = 0.05
#plt.annotate('', xy=(end, sign*arr_off), xycoords = 'data',
             #xytext = (end+kappa*dt, sign*arr_off), textcoords = 'data',
             #arrowprops = dict(arrowstyle='<->'))
#plt.text(end+.5*kappa*dt, sign*text_off, "$\kappa \Delta T$", horizontalalignment='center')

# Frequency
#ax2 = ax1.twinx()
#DeltaF0 = dt * np.array([np.sum(DeltaF1[:i]) + .5 * DeltaF1[i]
#                                      for i in range(N)])
#p2, = ax2.step(time, DeltaF0, where="post",
#         color="g", lw=2)

# Set ticks and labels
#ax2_ymax = 1.1*max(abs(DeltaF0.min()), abs(DeltaF0.max()))
#ax2.set_yticks(np.linspace(-ax2_ymax, ax2_ymax, 3))
#ax2.set_yticklabels([])

ax1_ymax = 1.1*max(abs(DeltaF1.min()), abs(DeltaF1.max()))
ax1.set_yticks(np.linspace(-ax1_ymax, ax1_ymax, 3))

ax1.set_xlim(0, N-1.1)
ax1.set_xlabel("time")

n = len(time)-2
reftimes = [time[i]+kappa*dt for i in range(n)]
reftimes.append(end+kappa*dt)
xticklabels = [r"$t_\mathrm{{ref}}^{{i-{}}}$".format(n-i) for i in range(n)]
xticklabels.append("$t_\mathrm{ref}^i = t^i$")
ax1.tick_params(direction='in', pad=10)
ax1.set_xticks(reftimes)
ax1.set_xticklabels(xticklabels)
ax1.set_yticks([0])
ax1.set_yticklabels(["0"])

ax1.set_ylabel("spin-down rate", rotation="vertical",
              labelpad=18.0)
r = Rectangle((0, 0), 1, 1, color="g", alpha=.5) # creates rectangle patch for legend use.
plt.legend([p0, r],
           [r"$\Delta\dot{f}$",
            #r"$t_{\textrm{ref}}$",
            r"$\Delta f(t)= \int \Delta\dot{f}(t)dt$",
            #r"$\Delta f_{i}$"
            ],
            loc=2, frameon=False)
plt.tight_layout()
plt.savefig("Illustration_F1_int.pdf")
