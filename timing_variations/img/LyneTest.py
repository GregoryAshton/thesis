import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz

plt.style.use("thesis")

def Perera(time, nudot1, nudot2, nuddot, T, tA, tB, tC, phi0):
    Spindown = np.zeros(len(time)) + nudot1

    ti_mod = np.mod(time + phi0*T, T)
    Spindown[(tA < ti_mod) & (ti_mod < tA+tB)] = nudot2
    Spindown[tA + tB + tC < ti_mod] = nudot2
    Spindown = Spindown + nuddot * time

    return Spindown


def PhaseFromSpindown(time, F1):

    F0 = cumtrapz(y=F1, x=time, initial=0)

    P0 = 2 * np.pi * cumtrapz(y=F0, x=time, initial=0)

    return P0


def SpindownFromPhase(sampling_times, time, Phase, DeltaT=100, frac=0.25):

    dt = time[1] - time[0]
    DeltaT_idx = int(DeltaT / dt)
    DeltaT_idx += DeltaT_idx % 2  # Make it even

    tref = time[0]
    timeprime = time - tref

    idxs = closest_idxs(time, sampling_times)

    vert_idx_list = idxs
    hori_idx_list = np.arange(-DeltaT_idx/2, DeltaT_idx/2)
    A, B = np.meshgrid(hori_idx_list, vert_idx_list)
    idx_array = A + B

    timeprime_array = timeprime[idx_array]

    Phase_array = Phase[idx_array]

    F1 = np.polynomial.polynomial.polyfit(timeprime_array[0],
                                          Phase_array.T, 2)[2, :]/np.pi

    return sampling_times, F1


def closest_idxs(array, points):
    """ Closests indexes of points in array

    Note: If the array is too large, this can quickly crash the system
    """
    points = points.reshape(len(points), 1)
    deltas = np.abs(points - array)
    return np.argmin(deltas, axis=1)


def SignalModel(time, nudot1, nudot2, nuddot, T, tA, tB, tC, phi0, DeltaT):
    sampling_times = time

    Nfine = 1000
    frac = 0.25  # Fraction to shift by

    time_fine = np.linspace(time[0]-DeltaT/2., time[-1]+DeltaT/2., Nfine)
    spindown = Perera(time_fine, nudot1, nudot2, nuddot, T, tA, tB, tC,
                      phi0)
    phase = PhaseFromSpindown(time_fine, spindown)
    time_ave, spindown_ave = SpindownFromPhase(sampling_times,
                                               time_fine,
                                               phase,
                                               DeltaT, frac)

    return spindown_ave

def MeasureHeightsAtShorterSwitch(DeltaT):
    UnderLyingSpindown = Perera(time, nudot1, nudot2, nuddot, T, tA, tB, tC, phi0)
    ObservedSpindown = SignalModel(time, nudot1, nudot2, nuddot, T, tA, tB, tC, phi0, DeltaT)
    return UnderLyingSpindown[idx_to_measure], ObservedSpindown[idx_to_measure]

nudot1 = -1e-13
nudot2 = -3e-13
T = 350*86400.
tA = 100*86400.
tB = 100*86400.
tC = 50*86400.
tD = T - (tA + tB + tC)
DeltaT = 100*86400
nuddot = 0
phi0 = 0

time = np.linspace(0, T, 50000)
t_to_measure = tA + tB + tC/2.
idx_to_measure = np.argmin(np.abs(time - t_to_measure))
UnderLyingSpindown = Perera(time, nudot1, nudot2, nuddot, T, tA, tB, tC, phi0)
ObservedSpindown = SignalModel(time, nudot1, nudot2, nuddot, T, tA, tB, tC, phi0, DeltaT)

obsH, actH = MeasureHeightsAtShorterSwitch(DeltaT)

fig, ax = plt.subplots()
ax.plot(time/86400, UnderLyingSpindown, color="k",
        label=r" $\dot{\nu}$")
ax.plot(time/86400, ObservedSpindown, ":", color="k",
        label=r"$\langle\dot{\nu}\rangle$")

ax.plot(time[idx_to_measure]/86400, obsH, "or")
ax.plot(time[idx_to_measure]/86400, actH, "ob")


Td = T / 86400
tAd = tA / 86400
tBd = tB / 86400
tCd = tC / 86400
tDd = tD / 86400
#ax.annotate("$T$", xy=(T*(0.5 - 0.02), 0.45*nudot1))
#ax.annotate("", xy=(0, 0.5*nudot1), xytext=(T, 0.5*nudot1),
#            arrowprops=dict(arrowstyle="<|-|>", facecolor='k', linewidth=1.0))
deltax = -12
dy = 0.9
dy2 = 0.8
ax.annotate("$t_{\mathrm{A}}$", xy=(.5*tAd + deltax, dy2*nudot1) )
ax.annotate("", xy=(0, dy*nudot1), xytext=(tAd, dy*nudot1),
            arrowprops=dict(arrowstyle="<|-|>, head_width=0.1, head_length=0.1"))
ax.annotate("$t_{\mathrm{B}}$", xy=(tAd + .5*tBd + deltax, dy2*nudot1) )
ax.annotate("", xy=(tAd, dy*nudot1), xytext=(tAd+tBd, dy*nudot1),
            arrowprops=dict(arrowstyle="<|-|>, head_width=0.1, head_length=0.1"))
ax.annotate("$t_{\mathrm{C}}$", xy=(tAd + tBd + .5*tCd + deltax, dy2*nudot1) )
ax.annotate("", xy=(tAd+tBd, dy*nudot1), xytext=(tAd+tBd+tCd, dy*nudot1),
            arrowprops=dict(arrowstyle="<|-|>, head_width=0.1, head_length=0.1"))
ax.annotate("$t_{\mathrm{D}}$", xy=(tAd + tBd + tCd + 0.5*tDd + deltax, dy2*nudot1) )
ax.annotate("", xy=(tAd+tBd+tCd, dy*nudot1), xytext=(Td, dy*nudot1),
            arrowprops=dict(arrowstyle="<|-|>, head_width=0.1, head_length=0.1"))


#ax.set_ylabel("Spindown rate [$\mathrm{s}^{-2}$]")
ax.set_xlabel("time [days]")
ax.set_yticks([nudot2, nudot1])
ax.set_yticklabels([r"$\dot{\nu}_{1}$", r"$\dot{\nu}_{2}$"])
ax.set_ylim(1.1*nudot2, 0.5*nudot1)
ax.legend(loc=3, frameon=False)
plt.tight_layout()
fig.savefig("TestLyneUnderlying.pdf")

DeltaT_list = np.linspace(10*86400, 150*86500, 100)
ratio = []
for DeltaT in DeltaT_list:
        obsH, actH = MeasureHeightsAtShorterSwitch(DeltaT)
        ratio.append(obsH/actH)

ax = plt.subplot(111)
ax.plot(DeltaT_list/86400, ratio, "-k")

ax.set_ylabel("$\mathcal{R}$", rotation="horizontal")
ax.set_xlabel("$T$ the time averaging baseline [days]")
plt.tight_layout()
plt.savefig("TestLyne.pdf")

