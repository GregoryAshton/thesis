import matplotlib.pyplot as plt
import numpy as np

plt.style.use("thesis")


def Beta_Function(epsI, epsA, chi):
    a = epsA*epsA+epsI*epsI-2*epsA*epsI*np.cos(2*chi)
    return np.arctan((epsI-epsA*np.cos(2*chi)-np.sqrt(a)) /
                     (2*epsA*np.sin(chi)*np.cos(chi)))

# Set some values
N = 100
chi = 30.0*np.pi/180

# Plot chi etc.
fig = plt.figure(figsize=(4,  3))
ax = plt.subplot(111)
plt.axhline(chi*180/np.pi, ls="-", color="b")
plt.axhline(chi*180/np.pi-90, ls="--", color="r")
ax.text(1.5e-3,  chi*180/np.pi*1.15, "$\chi$", color="b")
ax.text(1.5e-3,  (chi*180/np.pi*1.15-90), "$\chi-90^{\circ}$", color="r")

# Plot the line $\tau_{A}=\tau_{P}$
plt.axvline(1.0, ls="--", color="k")
plt.text(1.1, 75, r"$\tau_{\textrm{A}}=\tau_{\textrm{P}}$",
         color="k", fontsize=8)

# Plot the dafault beta for prolate and oblate cases

j = 0
labels = [r"$\beta(\epsilon_\mathrm{I}<0)$",
          r"$\beta (\epsilon_\mathrm{I}>0)$"]
styles = ["--", "-"]

for epsI in [-1e-6, 1e-6]:
    epsA_list = np.logspace(-9, -3, N)
    ratio_list = []
    beta_list_G = []
    beta_list_J = []
    for i in range(N):
        ratio_list.append(abs(epsA_list[i]/epsI))
        beta = Beta_Function(epsI, epsA_list[i], chi)
        beta_list_G.append((beta)*180/np.pi)

    plt.semilogx(ratio_list,  beta_list_G,  color="k",  ls=styles[j],
                 label=labels[j])
    j += 1

plt.yticks(np.arange(-90, 105, 15))

plt.grid(True,  linewidth=0.1,  linestyle="-")

plt.xlabel(r"$|\epsilon_\mathrm{A}/\epsilon_\mathrm{I}|$")
plt.ylabel(r"$\beta = z \; \angle \; z' \; $ [degs]")

leg = ax.legend(loc=1, numpoints=1)
leg.draggable()
leg.draw_frame(False)
plt.xlim(pow(10, -3), pow(10, 3))

plt.tight_layout()
plt.savefig("img/beta.png", dpi=500)

