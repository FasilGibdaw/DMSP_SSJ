# use DMSP reader to plot polar map from many days
import dmsp_reader as dmsp
import numpy as np
import glob
files = glob.glob('*dmsp*.cdf')
files.sort()
#fname = files[-20]
DD = np.zeros([1, 4])[0]
for fname in files[0:30]:
    D = dmsp.dmsp_reader(fname, channel=5)
    DD = np.vstack([D, DD])
x, y, egrided_data, igrided_data = dmsp.dmsp_grid(DD)
dmsp.dmsp_polar_plot(x, y, egrided_data, igrided_data)
