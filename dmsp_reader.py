import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from datetime import datetime
from matplotlib.pyplot import text
import pandas as pd
import cdflib
import numpy as np
import glob
import requests
# files=glob.glob('*.cdf')
# files.sort()
# fname=files[100]


def main():
    file = 'https://cdaweb.gsfc.nasa.gov/pub/data/dmsp/dmspf18/ssj/precipitating-electrons-ions/2012/dmsp-f18_ssj_precipitating-electrons-ions_20120101_v1.1.1.cdf'
    r = requests.get(file, allow_redirects=True)
    fname = 'dmsp-f18_ssj_precipitating-electrons-ions_20120101_v1.1.1.cdf'
    open(fname, 'wb').write(r.content)
    D = dmsp_reader(fname, channel=9)
    x, y, grided_data = dmsp_grid(D)
    dmsp_polar_plot(x, y, grided_data)


def dmsp_reader(fname, channel=1):
    '''reads DMSP cdf file from https://cdaweb.gsfc.nasa.gov/pub/data/dmsp/'''
    cdf_file = cdflib.CDF(fname)
    flux = cdf_file.varget("ELE_DIFF_ENERGY_FLUX")
    ch_energy = cdf_file.varget("CHANNEL_ENERGIES")
    mlt = cdf_file.varget("SC_AACGM_LTIME")
    mlat = cdf_file.varget("SC_AACGM_LAT")
    mlon = cdf_file.varget("SC_AACGM_LON")
    epoch = cdf_file.varget("Epoch")
    ind = np.where(mlat > 40)
    # flux from specified channel (30keV, 20keV --- 30eV) order
    Data = flux[ind, channel]
    # MLT=mlt[ind]
    # MLAT=mlat[ind]
    MLON = mlon[ind]
    EPOCH = epoch[ind]
    return np.vstack([mlt[ind], mlat[ind], Data]).T


index_labels = ['MLT', 'MLAT', 'Prep']


def dmsp_grid(D):
    df = pd.DataFrame(D, columns=index_labels)
    df = df.round()
    # data=df.sort_values(by=['MLT','MLAT'])
    # df['MLT']=df['MLT'].round()
    # df['MLAT']=df['MLAT'].round()
    f = df.groupby(['MLT', 'MLAT'])['Prep'].apply(
        lambda g: g.mean(skipna=True))
    f = f.to_frame().reset_index()

    lati = np.arange(f['MLAT'].min(), f['MLAT'].max()+1)
    tt = np.arange(f['MLT'].min(), f['MLT'].max()+1)
    x, y = np.meshgrid(tt, lati)

    mlt = f['MLT']
    latitude = f['MLAT']
    Prep = f['Prep']
    grided_data = np.nan*(np.ones(x.shape))
    for i in range(len(x[0, :])):
        ii = np.where(np.in1d(np.array(mlt), x[0, i]))
        # np.array(latitude[ii[0]])
        mm = np.where(np.in1d(y[:, 0], np.array(latitude[ii[0]])))
        grided_data[mm, i] = Prep[ii[0]]
    return x, y, grided_data


def dmsp_polar_plot(x, y, grided_data):
    fig = plt.figure(figsize=[10, 5])
    #ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.SouthPolarStereo())
    ax2 = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
    fig.subplots_adjust(bottom=0.05, top=0.95,
                        left=0.04, right=0.95, wspace=0.02)
    # Limit the map to -60 degrees latitude and below.
    ax2.set_extent([-180, 180, 90, 40], ccrs.PlateCarree())
    # ax1.add_feature(cfeature.LAND)
    # ax1.add_feature(cfeature.OCEAN)
    # ax1.gridlines()
    ax2.gridlines()
    # ax2.add_feature(cfeature.LAND)
    # ax2.add_feature(cfeature.OCEAN)
    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    # ax2.pcolor(mlt['MLT'],la['La'],mat['avPrec'])
    c = ax2.pcolor(x*15, y, np.log10(grided_data),
                   transform=ccrs.PlateCarree(), cmap='jet', vmin=4, vmax=8)
    #c=ax2.pcolor(np.transpose(xv),np.transpose(yv),np.log10(pp),transform=ccrs.PlateCarree(),cmap='coolwarm',vmin=4, vmax=8)
    fig.colorbar(c)
    # ti = datetime.strptime(fname[-19:-11], '%Y%m%d').date()
    # sat = fname[0:8]
    # text(0.1, 1, sat, ha='center', va='center', transform=ax2.transAxes)
    # ax2.set_title(ti)
    plt.show()


if __name__ == '__main__':
    main()
