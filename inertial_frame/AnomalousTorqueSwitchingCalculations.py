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

R = 1e6
c = 3e10

def BetaPrime(epsI, epsA, upsilon, chi):
    chi = np.radians(chi)
    epsA = upsilon * epsA
    top = epsI - epsA * np.cos(2*chi) - np.sqrt(epsA**2 + epsI**2 - 2 * epsA * 
                                                       epsI * np.cos(2 * chi))
    bottom = 2 * epsA * np.sin(chi) * np.cos(chi)

    return np.arctan2(top, bottom)

if "BetaPlot" in sys.argv:

    ax = plt.subplot(111)
    chi = 65
    epsI = 1e-3
    epsA = np.logspace(-5, -1, 100)

    for upsilon, upsilon_label in zip([0.2, 0.5, 0.8], 
                                  [r"$\frac{1}{5}$", r"$\frac{1}{2}$", r"$\frac{4}{5}$"]):
        DeltaBeta = np.abs(BetaPrime(epsI, epsA, 1, chi) - 
                           BetaPrime(epsI, epsA, upsilon, chi))

        ax.semilogx(epsA/epsI, np.degrees(DeltaBeta), 
                    label="$\upsilon=$" + upsilon_label, rasterized=True)
        ax.set_xlabel(r"$\frac{\epsilon_{A}}{\epsilon_{I}}$")
        ax.set_ylabel(r"$\Delta \beta = |\beta' - \beta|$ [degrees]", rotation="vertical")

    ax.legend(loc=1, frameon=False, fontsize=20)
    plt.tight_layout()
    plt.savefig("img/DeltaBetaPlot.pdf", dpi=270)
    plt.show()
