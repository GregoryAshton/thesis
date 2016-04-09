import numpy as np
import pandas as pd
from matplotlib import rc_file
rc_file("/home/greg/Neutron_star_modelling/matplotlibrc")

DATA_FILE = "ATNF_data_file.txt"

data = np.genfromtxt(DATA_FILE, skip_header=4, skip_footer=1, dtype=None)
name = data[:, 0]
F0 = np.genfromtxt(data[:, 1])
F1 = np.genfromtxt(data[:, 2])
F2 = np.genfromtxt(data[:, 3])
Binary = data[:, 4]  # If not "*" then it is a binary
Type = data[:, 5]
W10_ms = data[:, 6]
tauAge = data[:, 7]
tauAge[tauAge == "*"] = np.nan
tauAge = np.array(tauAge)

df = pd.DataFrame({'name': name,
                   'F0': F0,
                   'F1': F1,
                   'F2': F2,
                   'Binary': Binary,
                   'Type': Type,
                   'W10_ms': W10_ms,
                   'tauAge': tauAge
                   })

# Clean data
df = df[df.Binary == "*"]  # remve binarys
df = df[df.Type == "*"]  # only look at ordinary pulsars
df = df[df.F0 < 5e1]  # lazy cut of MSPs TBI

df = df.replace("*", np.nan)
df['W10_ms'] = df['W10_ms'].astype(float)
df['W10'] = df['W10_ms'] * 1e-3

