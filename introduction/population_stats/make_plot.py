import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.stats import norm
import BDATools as BDA
import matplotlib.ticker as ticker
import scipy.stats as ss
from get_data import df
plt.style.use("thesis")


def AddNormal(ax, data):
    mu = np.mean(data)
    std = np.std(data)
    dmin, dmax = data.min(), data.max()
    drange = dmax - dmin
    quants = np.linspace(dmin-0.2*drange, dmax+0.2*drange, 200)

    #kde = ss.gaussian_kde(data)
    #ax.plot(quants, kde.pdf(quants), "-b", label="KDE")

    s = "Norm. fit:\n $\mu$={}\n $\sigma$={}".format(
        BDA.Texify_Float(mu, 2), BDA.Texify_Float(std, 2))
    s = s.replace("\\text", "")
    ax.plot(quants, norm.pdf(quants, loc=mu, scale=std), "-", lw=1, color="r",
            label=s, dashes=(2, 1))

    return ax


def setup_hist(vals, ax, nbins=30):
    #ax.hist(vals, nbins, normed=True, histtype="step")
    ax.hist(vals, nbins, normed=True, color="grey", histtype="stepfilled",
            alpha=0.5, label="histogram")
    return ax

# timing properties
fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(6, 2.5))

logF0 = np.log10(df.F0[np.isfinite(df.F0.values)].values)
print("len(F0): {}".format(len(logF0)))
ax1 = setup_hist(logF0, ax1, 50)
ax1 = AddNormal(ax1, logF0)
ax1.set_xlabel(r"$\log_{10}(f/\textrm{s}^{-1})$")
ax1.set_ylabel("Normalised count")

F1 = np.abs(df.F1[np.isfinite(df.F1.values)]).values
print("len(F1): {}".format(len(F1)))
logF1 = np.log10(F1[F1 > 0])
ax2 = setup_hist(logF1, ax2, 50)
ax2 = AddNormal(ax2, logF1)
ax2.set_xlabel(r"$\log_{10}(\dot{f}/\textrm{s}^{-2})$")

logF2 = np.log10(np.abs(df.F2[np.isfinite(df.F2.values)]).values)
print("len(F2): {}".format(len(logF2)))
ax3 = setup_hist(logF2, ax3, 20)
ax3 = AddNormal(ax3, logF2)
ax3.set_xlabel(r"$\log_{10}(\ddot{f}/\textrm{s}^{-3})$")

for ax in (ax1, ax2, ax3):
    ax.set_yticks([])
    ax.legend(frameon=False, fontsize=5, loc=1)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))

fig.tight_layout()
plt.savefig("timing_distribution.pdf")


fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(6, 2.5))

tauAge = df.tauAge.values.astype(np.float64)
print("len(tauAge): {}".format(len(tauAge)))
tauAge = tauAge[~np.isnan(tauAge)] # Yrs
log10tauAge = np.log10(tauAge)
ax1 = setup_hist(log10tauAge, ax1, 30)
ax1 = AddNormal(ax1, log10tauAge)
ax1.set_xlabel(r"$\log_{10}(\tau_{\textrm{c}}/\textrm{yrs})$")
ax1.set_ylabel("Normalised count")

log10W10 = np.log10(df.W10[np.isfinite(df.W10)]) # s
print("len(W10): {}".format(len(log10W10)))
ax2 = setup_hist(log10W10, ax2, 50)
ax2 = AddNormal(ax2, log10W10)
ax2.set_xlabel(r"$\log_{10}(W_{10}/\textrm{s})$")

df['duty'] = df.W10 * df.F0
log10duty = np.log10(df.duty.values)
log10duty = log10duty[~np.isnan(log10duty)]
print("len(duty): {}".format(len(log10duty)))
ax3 = setup_hist(log10duty, ax3, 50)
ax3 = AddNormal(ax3, log10duty)
ax3.set_xlabel(r"$\log_{10}(W_{10} f)$")

for ax in (ax1, ax2, ax3):
    ax.set_yticks([])
    ax.legend(frameon=False, fontsize=5, loc=1)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))

fig.tight_layout()
plt.savefig("W10_and_age_distribution.pdf")

