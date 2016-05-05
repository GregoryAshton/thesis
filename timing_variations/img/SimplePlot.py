import matplotlib.pylab as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.integrate import cumtrapz
#from nsmod import Plot
#import TNtools as TN
plt.style.use("thesis")

nu_dot_1 = -365.68e-15
nu_dot_2 = 0.71 * nu_dot_1

R = 2.0
D = 1e-10
D = 0.3
T = 150 * 24 * 3600

t1 = 100 * 24 * 3600
t2 = R * t1

N = 100

t1_list = np.random.normal(t1, D*t1, size=N)
t2_list = np.random.normal(t2, D*t2, size=N)

flip_markers = [sum(np.dstack((t1_list, t2_list)).flat[0:i])
                for i in range(2*N)]

t_list = np.linspace(0, 15 * t1, 4000)

nu_dot_list = []
current_nu_dot, other_nu_dot = nu_dot_1, nu_dot_2
j = 0
for t in t_list:
    if t < flip_markers[j]:
        nu_dot_list.append(current_nu_dot)
    if t >= flip_markers[j]:
        current_nu_dot, other_nu_dot = other_nu_dot, current_nu_dot
        j += 1
        nu_dot_list.append(current_nu_dot)

# Plotting
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True, figsize=(4.5, 3.0))
ax1.plot(t_list, nu_dot_list, "-k")
ax1.set_ylabel(r"$\dot{\nu}$")
ax1.set_xticklabels([])


nu = 2 + cumtrapz(y=nu_dot_list, x=t_list, initial=0)
phase = 2*np.pi*cumtrapz(y=nu, x=t_list, initial=0)

dn = np.argmin(np.abs(t_list - T))
ilist = np.arange(0, len(phase)-dn, 1)
nudot_ave = [np.polyfit(t_list[i:i+dn], phase[i:i+dn], 2)[0]/(np.pi) for i in ilist]
t_mids = [np.mean(t_list[i:i+dn]) for i in ilist]
ax2.plot(t_mids, nudot_ave, "-k")
ax2.set_ylabel(r"$\langle\dot{\nu}\rangle$")

for ax in [ax1, ax2]:
    ax.set_yticks([nu_dot_1, nu_dot_2])
    ax.set_yticklabels([r"$\dot{\nu}_{1}$", r"$\dot{\nu}_{2}$"])


# Fit polynomial to phase of order order
coefs = np.polyfit(t_list, phase, 3)
# poly1d returns the polyn we then evaluate this at time giving the fitted phi
phase_fit = np.poly1d(coefs)(t_list)

# Subtract the two to get a residual
phase_res = (phase - phase_fit) / nu / (2*np.pi)

ax3.plot(t_list, phase_res, "-k")
ax3.set_ylabel("Phase residual \n [cycles]")
ax3.set_xlabel(r"time")
ax3.set_yticks(ax3.get_yticks()[0:-1])
ax3.yaxis.set_major_locator(MaxNLocator(5))

ax1.set_xlim(0, t_list[-100])
ax3.set_xlim(0, t_list[-100])

#ax3.set_xticks([t1, t1+t2, 2*t1+t2])
#ax3.set_xticklabels(["$t_A$", "$t_A + t_B$", "$2t_A + t_B$"], rotation=30)
T = t1 + t2
ax3.set_xticks([0, t1] + [i*T for i in range(1, 5)])
ax3.set_xticklabels([0, "$t_\mathrm{A}$", "$t_\mathrm{A} + t_\mathrm{B}$"]
                    + ["${}(t_\mathrm{{A}} + t_\mathrm{{B}})$".format(i) for i in range(2, 6)])

if D < 1e-9:
    D = 0
ax1.set_title("$R={}$,    $D={}$".format(R, D))

plt.tight_layout()
plt.savefig("R_{}_D_{}.pdf".format(R, D))
