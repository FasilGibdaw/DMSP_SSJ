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
import warnings
# files=glob.glob('*.cdf')
# files.sort()
# fname=files[100]
warnings.filterwarnings("ignore")


def main():
    fname = 'dmsp-f18_ssj_precipitating-electrons-ions_20120101_v1.1.1.cdf'
    try:
        D = dmsp_reader(fname, channel=9)
        x, y, egrided_data, igrided_data = dmsp_grid(D)
        dmsp_polar_plot(x, y, egrided_data, igrided_data)
    except:
        file = 'https://cdaweb.gsfc.nasa.gov/pub/data/dmsp/dmspf18/ssj/precipitating-electrons-ions/2012/dmsp-f18_ssj_precipitating-electrons-ions_20120101_v1.1.1.cdf'
        r = requests.get(file, allow_redirects=True)
        open(fname, 'wb').write(r.content)
        D = dmsp_reader(fname, channel=9)
        x, y, egrided_data, igrided_data = dmsp_grid(D)
        dmsp_polar_plot(x, y, egrided_data, igrided_data)


def dmsp_reader(fname, channel=1):
    '''reads DMSP cdf file from https://cdaweb.gsfc.nasa.gov/pub/data/dmsp/'''
    cdf_file = cdflib.CDF(fname)
    eflux = cdf_file.varget("ELE_DIFF_ENERGY_FLUX")
    iflux = cdf_file.varget('ION_DIFF_ENERGY_FLUX')
    ch_energy = cdf_file.varget("CHANNEL_ENERGIES")
    mlt = cdf_file.varget("SC_AACGM_LTIME")
    mlat = cdf_file.varget("SC_AACGM_LAT")
    mlon = cdf_file.varget("SC_AACGM_LON")
    epoch = cdf_file.varget("Epoch")
    ind = np.where(mlat > 40)
    # flux from specified channel (30keV, 20keV --- 30eV) order
    eData = eflux[ind, channel]
    iData = iflux[ind, channel]

    # MLT=mlt[ind]
    # MLAT=mlat[ind]
    MLON = mlon[ind]
    EPOCH = epoch[ind]
    return np.vstack([mlt[ind], mlat[ind], eData, iData]).T


index_labels = ['MLT', 'MLAT', 'ePrep', 'iPrep']


def dmsp_grid(D):
    df = pd.DataFrame(D, columns=index_labels)
    df = df.round()
    # data=df.sort_values(by=['MLT','MLAT'])
    # df['MLT']=df['MLT'].round()
    # df['MLAT']=df['MLAT'].round()
    f = df.groupby(['MLT', 'MLAT'])['ePrep', 'iPrep'].apply(
        lambda g: g.mean(skipna=True))
    f = f.reset_index()

    lati = np.arange(f['MLAT'].min(), f['MLAT'].max()+1)
    tt = np.arange(f['MLT'].min(), f['MLT'].max()+1)
    x, y = np.meshgrid(tt, lati)

    mlt = f['MLT']
    latitude = f['MLAT']
    ePrep = f['ePrep']
    iPrep = f['iPrep']
    egrided_data = np.nan*(np.ones(x.shape))
    igrided_data = np.nan*(np.ones(x.shape))
    for i in range(len(x[0, :])):
        ii = np.where(np.in1d(np.array(mlt), x[0, i]))
        # np.array(latitude[ii[0]])
        mm = np.where(np.in1d(y[:, 0], np.array(latitude[ii[0]])))
        egrided_data[mm, i] = ePrep[ii[0]]
        igrided_data[mm, i] = iPrep[ii[0]]
    return x, y, egrided_data, igrided_data


def dmsp_polar_plot(x, y, egrided_data, igrided_data):
    fig = plt.figure(figsize=[10, 5])
    #ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.SouthPolarStereo())
    ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.NorthPolarStereo())
    fig.subplots_adjust(bottom=0.05, top=0.95,
                        left=0.04, right=0.95, wspace=0.02)
    # Limit the map to -60 degrees latitude and below.
    ax1.set_extent([-180, 180, 90, 40], ccrs.PlateCarree())
    ax1.gridlines()
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    c = ax1.pcolor(x*15, y, np.log10(egrided_data),
                   transform=ccrs.PlateCarree(), cmap='jet', vmin=4, vmax=8)
    fig.colorbar(c)

    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.NorthPolarStereo())
    fig.subplots_adjust(bottom=0.05, top=0.95,
                        left=0.04, right=0.95, wspace=0.02)
    # Limit the map to -60 degrees latitude and below.
    ax2.set_extent([-180, 180, 90, 40], ccrs.PlateCarree())
    ax2.gridlines()
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    c = ax2.pcolor(x*15, y, np.log10(igrided_data),
                   transform=ccrs.PlateCarree(), cmap='jet', vmin=4, vmax=8)
    fig.colorbar(c)
    # ax1.add_feature(cfeature.LAND)
    # ax1.add_feature(cfeature.OCEAN)
    # ax1.gridlines()

    # ax2.add_feature(cfeature.LAND)
    # ax2.add_feature(cfeature.OCEAN)
    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.

    # ax2.pcolor(mlt['MLT'],la['La'],mat['avPrec'])

    #c=ax2.pcolor(np.transpose(xv),np.transpose(yv),np.log10(pp),transform=ccrs.PlateCarree(),cmap='coolwarm',vmin=4, vmax=8)

    # ti = datetime.strptime(fname[-19:-11], '%Y%m%d').date()
    # sat = fname[0:8]
    # text(0.1, 1, sat, ha='center', va='center', transform=ax2.transAxes)
    # ax2.set_title(ti)
    plt.show()


if __name__ == '__main__':
    main()
