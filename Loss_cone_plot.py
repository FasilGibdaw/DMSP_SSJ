# simple loss cone versus L vlaue plot in a dipole field
import numpy as np
import matplotlib.pyplot as plt
L = np.arange(1.1, 12, 0.1)
# found at https://www.ucl.ac.uk/~ucapnac/AstrophysicalDiscsHomeworks.pdf
dd = np.sqrt(L**-3*1/(np.sqrt(3*(1-(1/L))+1)))
plt.plot(L, np.rad2deg(np.arcsin(dd)))
plt.xlabel('L')
plt.ylabel('Loss cone angle (deg.)')
