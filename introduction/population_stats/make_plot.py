import matplotlib.pyplot as plt
import numpy as np
import sys, os
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

    statistic, pvalue = ss.mstats.normaltest(data)

    fit_text = "mean$=${}\n std. dev.$=${}".format(
        BDA.Texify_Float(mu, 2), BDA.Texify_Float(std, 2))
    fit_text = fit_text.replace("\\text", "")
    ax.annotate(fit_text, xy=(0.02, 0.89), xycoords="axes fraction", size=5)
    if pvalue < 0.01:
        pvaluepow = np.round(np.log10(pvalue))
        pvalue = r"\approx 10^{{ {:1.0f} }}".format(pvaluepow)
        pvalue = pvalue.replace("-", r"\textrm{-}")
    else:
        pvalue = "=" + BDA.Texify_Float(pvalue, 2)
    #s = "Norm. fit\n $p$-val.=\n${}$".format(BDA.Texify_Float(pvalue, 1))
    s = "Norm. fit\n $p$-val.${}$".format(pvalue)
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
ax1.set_xlabel(r"$\log_{10}(\nu/\textrm{s}^{-1})$")
ax1 = AddNormal(ax1, logF0)
ax1.set_ylabel("Normalised count")

F1 = np.abs(df.F1[np.isfinite(df.F1.values)]).values
print("len(F1): {}".format(len(F1)))
logF1 = np.log10(F1[F1 > 0])
ax2 = setup_hist(logF1, ax2, 50)
ax2.set_xlabel(r"$\log_{10}(\dot{\nu}/\textrm{s}^{-2})$")
ax2 = AddNormal(ax2, logF1)

logF2 = np.log10(np.abs(df.F2[np.isfinite(df.F2.values)]).values)
print("len(F2): {}".format(len(logF2)))
ax3 = setup_hist(logF2, ax3, 20)
ax3.set_xlabel(r"$\log_{10}(\ddot{\nu}/\textrm{s}^{-3})$")
ax3 = AddNormal(ax3, logF2)

for ax in (ax1, ax2, ax3):
    ax.set_yticks([])
    ax.set_ylim(0, 1.1*ax.get_ylim()[1])
    ax.legend(frameon=False, fontsize=5, loc=1)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))

fig.tight_layout()
fig.subplots_adjust(wspace=0.05)
plt.savefig("timing_distribution.pdf")


fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(6, 2.5))

tauAge = df.tauAge.values.astype(np.float64)
print("len(tauAge): {}".format(len(tauAge)))
tauAge = tauAge[~np.isnan(tauAge)] # Yrs
log10tauAge = np.log10(tauAge)
ax1 = setup_hist(log10tauAge, ax1, 30)
ax1.set_xlabel(r"$\log_{10}(\tau_{\mathrm{age}}/\textrm{yrs})$")
ax1 = AddNormal(ax1, log10tauAge)
ax1.set_ylabel("Normalised count")

log10W10 = np.log10(df.W10[np.isfinite(df.W10)]) # s
print("len(W10): {}".format(len(log10W10)))
ax2 = setup_hist(log10W10, ax2, 50)
ax2.set_xlabel(r"$\log_{10}(W_{10}/\textrm{s})$")
ax2 = AddNormal(ax2, log10W10)

df['duty'] = df.W10 * df.F0
log10duty = np.log10(df.duty.values)
log10duty = log10duty[~np.isnan(log10duty)]
print("len(duty): {}".format(len(log10duty)))
ax3 = setup_hist(log10duty, ax3, 50)
ax3.set_xlabel(r"$\log_{10}(W_{10} \nu)$")
ax3 = AddNormal(ax3, log10duty)

for ax in (ax1, ax2, ax3):
    ax.set_yticks([])
    ax.set_ylim(0, 1.1*ax.get_ylim()[1])
    ax.legend(frameon=False, fontsize=5, loc=1)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))

fig.tight_layout()
plt.savefig("W10_and_age_distribution.pdf")

