import matplotlib.pyplot as plt
import numpy as np

plt.style.use('thesis')
fig, ax = plt.subplots()

F0 = 1
F1 = -0.04
time = np.linspace(-1, 1, 1000)
tglitch = 0
tau = 0.05
dF0 = 0.03
dF0R = 0.05

F = F0 + F1*time
d = dF0 + np.exp(-(time - tglitch)/tau)*dF0R


F[time > tglitch] += d[time > tglitch]


ax.plot(time, F, "-k")
ax.set_xticks([0])
ax.set_xticklabels(['$t_{g}$'])

xval = -0.05
xoff = -0.2
ax.annotate("", xy=(xval, F0), xytext=(xval, F0+dF0+0.005),
            arrowprops=dict(arrowstyle="<|-|>", fc='k' ))
ax.annotate(r"$\Delta\nu$", xy=(xval+xoff, F0+0.5*dF0),
            xytext=(xval+xoff, F0+0.5*dF0))
ax.annotate("", xy=(xval, F0+dF0), xytext=(xval, F0+dF0+dF0R),
            arrowprops=dict(arrowstyle="<|-|>", fc='k' ))
ax.annotate(r"$\Delta\nu_\mathrm{t}$", xy=(xval+xoff,
            F0+dF0+0.5*dF0R), xytext=(xval+xoff, F0+dF0+0.5*dF0R))
ax.set_yticklabels([])
ax.set_ylabel(r"Spin-frequency $\nu$")

ax.set_ylim(0.95, 1.1)
plt.savefig("img/glitch_sketch.png")
