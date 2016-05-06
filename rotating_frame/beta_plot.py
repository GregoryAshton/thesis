import matplotlib.pyplot as plt
import numpy as np
from math import pi
import sys

plt.style.use("thesis")

 # Set the defaults
plt.rcParams['axes.grid'] = False

## Beta function
def Beta_Function(epsI, epsA, chi):
    a=epsA*epsA+epsI*epsI-2*epsA*epsI*np.cos(2*chi)
    return np.arctan((epsI-epsA*np.cos(2*chi)-np.sqrt(a))/(2*epsA*np.sin(chi)*np.cos(chi)))

# Set some values
N = 100
chi = 30.0*pi/180

# Plot chi etc.
fig = plt.figure(figsize=(5, 3))
ax = plt.subplot(111)
plt.axhline(chi*180/pi,ls="-",color="b",lw=1.2)
plt.axhline(chi*180/pi-90,ls="--",color="r",lw=1.2)
ax.text(1.5e-3, chi*180/pi*1.15,"$\chi$",color="b")
ax.text(1.5e-3, (chi*180/pi*1.15-90),"$\chi-90^{\circ}$",color="r")

# Plot the line $\tau_{A}=\tau_{P}$
plt.axvline(1.0,ls="--",color="k",lw=1.2)
plt.text(1.1, 75, r"$\tau_{\textrm{A}}=\tau_{\textrm{P}}$", color="k", fontsize=8)

# Plot the dafault beta for prolate and oblate cases

j=0 ; labels=[r"$\beta(\epsilon_\mathrm{I}<0)$",r"$\beta (\epsilon_\mathrm{I}>0)$"] ;styles=["--","-"]
for epsI in [-1e-6 , 1e-6]:

    epsA_list=np.logspace(-9,-3,N)

    ratio_list=[] ; beta_list_G=[] ; beta_list_J=[]
    for i in range(N):
            ratio_list.append(abs(epsA_list[i]/epsI))
            beta = Beta_Function(epsI,epsA_list[i],chi)
            beta_list_G.append((beta)*180/pi)

    plt.semilogx(ratio_list, beta_list_G, color="k", ls=styles[j],
                label=labels[j],lw=2)
    j += 1


# Plot the Jones beta function if required
if "Jones_beta" in sys.argv:
    for epsI in [-1e-6 , 1e-6]:
        epsA_list=np.logspace(-9,-3,N)

        ratio_list=[] ; beta_list_G=[] ; beta_list_J=[]
        for i in range(N):
            ratio_list.append(abs(epsA_list[i]/epsI))
            beta_list_J.append(Beta_func_Jones(epsI,epsA_list[i],chi))

        plt.semilogx(ratio_list,beta_list_J,color="b",ls="--",label=r"$\beta_{2}$",lw=2.5)


plt.yticks(np.arange(-90,105,15))

# Horrible code to produce lines easily
lw_val=0.2
plt.axhline(0,ls="--",color="k",lw=lw_val)
plt.axhline(30,ls="--",color="k",lw=lw_val)
plt.axhline(60,ls="--",color="k",lw=lw_val)
plt.axhline(90,ls="--",color="k",lw=lw_val)
plt.axhline(-30,ls="--",color="k",lw=lw_val)
plt.axhline(-60,ls="--",color="k",lw=lw_val)
plt.axhline(-90,ls="--",color="k",lw=lw_val)

plt.axvline(pow(10,0),ls="--",color="k",lw=lw_val)
plt.axvline(pow(10,1),ls="--",color="k",lw=lw_val)
plt.axvline(pow(10,2),ls="--",color="k",lw=lw_val)
plt.axvline(pow(10,3),ls="--",color="k",lw=lw_val)
plt.axvline(pow(10,-1),ls="--",color="k",lw=lw_val)
plt.axvline(pow(10,-2),ls="--",color="k",lw=lw_val)
plt.axvline(pow(10,-3),ls="--",color="k",lw=lw_val)



if "Jones_beta" in sys.argv: plt.legend(loc=1)
plt.xlabel(r"$|\epsilon_\mathrm{A}/\epsilon_\mathrm{I}|$")
plt.ylabel(r"$\beta = z \; \angle \; z' \; $ [degs]")

#py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12,hspace=0.0)
leg = ax.legend(loc=1,numpoints=1)
leg.draggable()
leg.draw_frame(False)
plt.xlim(pow(10,-3),pow(10,3))

plt.tight_layout()
plt.savefig("img/beta.png")

