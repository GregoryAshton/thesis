import numpy as np
import BDATools as BDA


def m(Tobs, F0, F1):
    nu = F0 / 2.
    nudot = F1 / 2.
    P0 = 1 / nu
    P1 = -nudot / nu**2
    A = -1.4*np.log10(P0) + 0.8*np.log10(P1/1e-15) - 3.31
    SFNCRAB = 1.25e-22
    SFN = 10**(2*A) * SFNCRAB
    return SFN, np.pi**2/630 * SFN * Tobs**3

searches = {'S5 E@H all-sky': {'F0': 1190,
                               'F1': -2e-9,
                               'Tcoh': 25,
                               'Tobs': 694},
            'S5 E@H galactic center': {'F0': 496,
                                       'F1': -71e-9,
                                       'Tcoh': 11.5,
                                       'Tobs': 302},
            'S5 all-sky': {'F1': -0.89e-9,
                           'F0': 1000,
                           'Tcoh': 0.5,
                           'Tobs': 365},

            'VSR low-frequency all-sky': {'F1': -10e-9,
                                          'F0': 128,
                                          'Tcoh': 2.3,
                                          'Tobs': 185},
            'S5 supernova remnant (Cas A)': {'F1': -60.5e-9,
                                             'F0': 573,
                                              'Tcoh': np.nan,
                                              'Tobs': 8.4},
            }

table = r"& $\Tobs$ [days] & $S_{\textrm{FN}}$ & $E[\tilde{\mu}]$ \\ \hline" + "\n"

keys = ['S5 E@H all-sky', 'S5 E@H galactic center',
        'S5 all-sky', 'VSR low-frequency all-sky', 'S5 supernova remnant (Cas A)',
        ]
row = "{} & ${}$ & ${}$ & ${}$ \\\\\n"
for key in keys:
    d = searches[key]
    Tobs = d['Tobs'] * 86400
    Tcoh = d['Tcoh'] * 60 * 60
    F1 = d['F1']
    F0 = d['F0']

    SFN, mval = m(Tobs, F0, F1)
    table += row.format(key, Tobs/86400, BDA.Texify_Float(SFN), BDA.Texify_Float(mval))

table = table.rstrip("\\\\\n")
with open("past_searches_worst_cases.tex", "w+") as f:
    f.write(table)
