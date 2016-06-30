import matplotlib.pyplot as plt
import numpy as np
import sys

plt.style.use('thesis')
plt.rcParams['ytick.labelsize'] = 6
plt.rcParams['xtick.labelsize'] = 6

R = 1e6
c = 3e10
N = 100
I0 = 1e45
Bmax = 16
Bmin = 10
omega_0 = 1e4
Bs_list = np.logspace(Bmin, Bmax, N)
pre = pow(R, 5) / (2 * I0 * c * c)
epsA_list = [pre*pow(Bs_list[i], 2) for i in range(N)]
pre = (2*R/(3*c))*omega_0
B23 = [pre * epsA_list[i] for i in range(N)]

def base_plot(ax):
    # plot a line of tau_p = tau_s
    ax.loglog(Bs_list, B23, label=r"$\tau_{\mathrm{P}}=\tau_{\mathrm{S}}$",
               dashes=(3, 1.5), color="k")

    # Add labels
    ax.set_xlabel(r"$B_{0}$ [Gauss]")
    ax.set_ylabel(r"$\epsilon_{\mathrm{I}}$", rotation="horizontal")

    zero = np.zeros(N) + 1e-23
    ax.fill_between(Bs_list, zero, B23, facecolor="cyan", alpha=0.2)
    #ax.fill_between(Bs_list, epsA_list, np.zeros(N)+1, facecolor="white",
    #                alpha=1)

    # Region A Neutron star
    Bs = 1.3416407865e+13
    epsI = 1e-9
    ax.scatter(Bs, epsI, marker="v", s=10, label=r"Neutron star A", color="r",
               edgecolor="k")
    # Region B Neutron star
    Bs = 1.3416407865e+13
    epsI = 4.0e-11
    ax.scatter(Bs, epsI, marker="s", s=10, label=r"Neutron star B", color="r",
               edgecolor="k", zorder=10)
    # Region C Neutron star
    Bs = 1.3416407865e+13
    epsI = 1e-15
    ax.scatter(Bs, epsI, marker="o", s=10, label=r"Neutron star C", color="r",
               edgecolor="k")

    ax.legend(loc=2, scatterpoints=1, fontsize=8, frameon=False)
    ax.set_yticks(ax.get_yticks()[2:-1:2])
    ax.set_xlim(1e10, 1e16)
    ax.set_ylim(1e-22, 1e-1)
    plt.tight_layout()

    return ax

fig, ax = plt.subplots(figsize=(4, 3))
# plot a line of eps_I = eps_A (given by epsA list!)
ax.loglog(Bs_list, epsA_list, label=r"$\tau_{\mathrm{P}}=\tau_{\mathrm{A}}$",
          ls="-", color="k")

# Add text to denote regions
ax.text(5e12, 1e-5,
        r"$\tau_{\mathrm{S}} > \tau_{\mathrm{A}} > \tau_{\mathrm{P}}$",
        size=8, rotation=0)
ax.text(5e12, 1e-19,
        r"$\tau_{\mathrm{P}} > \tau_{\mathrm{S}} > \tau_{\mathrm{A}}$",
        size=8, rotation=0)
ax.fill_between(Bs_list, B23, epsA_list, facecolor="0.6", alpha=0.7)
ax = base_plot(ax)
fig.savefig("img/phase_space.png", dpi=500)


# Phase space without anomalous torque
fig, ax = plt.subplots(figsize=(4, 3))
ax.text(5e12, 1e-5, r"$\tau_{\mathrm{S}} > \tau_{\mathrm{P}}$",
        size=8, rotation=0)
ax.text(5e12, 1e-19, r"$\tau_{\mathrm{P}} > \tau_{\mathrm{S}} $",
        size=8, rotation=0)
ax = base_plot(ax)
fig.savefig("img/phase_space_no_anom.png", dpi=500)
