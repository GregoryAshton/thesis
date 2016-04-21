"""
Python script taken from 
/home/greg/Neutron_star_modelling/TwoStateSwitching/PrecessionInducedBySwitching

"""

import matplotlib.pyplot as plt
import numpy as np
from nsmod import Plot
import sys
from nsmod.manual_switching_torque_with_Euler import main
from nsmod import Plot
from nsmod.Physics_Functions import Beta_Function
from nsmod import File_Functions
plt.style.use("thesis")

R = 1e6
c = 3e10

def BetaPrime(epsI, epsA, Sa, chi):
    chi = np.radians(chi)
    epsA = (1-Sa) * epsA
    top = epsI - epsA * np.cos(2*chi) - np.sqrt(epsA**2 + epsI**2 - 2 * epsA * 
                                                epsI * np.cos(2 * chi))
    bottom = 2 * epsA * np.sin(chi) * np.cos(chi)

    return np.arctan2(top, bottom)


ax = plt.subplot(111)
chi = 65
epsI = 1e-3
epsA = np.logspace(-5, -1, 100)

for Sa, Sa_label in zip([0.2, 0.5, 0.8], 
                        [r"$\frac{1}{5}$", r"$\frac{1}{2}$", r"$\frac{4}{5}$"]):
    DeltaBeta = np.abs(BetaPrime(epsI, epsA, 0, chi) - 
                       BetaPrime(epsI, epsA, Sa, chi))

    ax.semilogx(epsA/epsI, np.degrees(DeltaBeta), 
                label=r"$S_{\mathrm{A}}=$" + Sa_label, rasterized=True)
    ax.set_xlabel(r"$\frac{\epsilon_{A}}{\epsilon_{I}}$")
    ax.set_ylabel(r"$\Delta \beta = |\beta' - \beta|$ [degrees]", rotation="vertical")


ax.legend(loc=1, frameon=False, fontsize=6)
plt.tight_layout()
plt.savefig("img/DeltaBetaPlot.pdf", dpi=270)

print "B1828-11"
nu = 2.47
nudot = -3.65e-13
c = 3e8
R = 10*1e3
epsA_B1828 = np.abs(c/(8*np.pi*R*nu**3)*nudot)
epsI_B1828 = np.logspace(-13, -9, 1000)
print "epsA = ", epsA_B1828
DeltaBeta = np.degrees(np.abs(BetaPrime(epsI_B1828, epsA_B1828, 0, chi) - 
                              BetaPrime(epsI_B1828, epsA_B1828, 0.0071, chi)))
print "Max Delta Beta = ", np.max(DeltaBeta)
print "Max epsI = ", epsI_B1828[np.argmax(DeltaBeta)]
fig, ax = plt.subplots()
ax.semilogx(np.abs(epsA_B1828/epsI_B1828), DeltaBeta)
plt.savefig("img/B1828.png")
