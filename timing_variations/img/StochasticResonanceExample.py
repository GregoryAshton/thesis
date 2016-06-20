import numpy as np
import matplotlib.pyplot  as plt
import matplotlib.ticker as ticker

plt.style.use('thesis')


def eom(x, t, Lambda=1, epsilon=0.001, D=1.5, omega0=1e2):
    # white noise has a mean of zero and a std of 1
    #noise = D*np.random.normal(0, 1)
    noise = np.random.normal(0, D)
    return Lambda * x - pow(x, 3) + noise + epsilon * np.cos(omega0 * t)

dt = 0.1
time = np.arange(0, 500, dt)

fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(6, 2))

for ax, (j, D) in zip([ax1, ax2, ax3], enumerate([0.1, 1.5, 6.6])):
    x = [1.0]
    for i, t in enumerate(time[:-1]):
        dx = eom(x[i], t, D=D)
        x.append(x[i]+dx*dt)
    ax.plot(time, x, color='k')
    ax.set_ylim(-2, 2)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.set_xticklabels([])
    ax.set_title("$\sigma$={}".format(D))

ax1.set_xlabel("time")
ax2.set_xlabel("time")
ax3.set_xlabel("time")

ax1.set_ylabel('$x$-position of particle')

plt.tight_layout()
plt.savefig("Stochastic_resonance.png")
