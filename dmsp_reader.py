import matplotlib.pyplot as plt
import matplotlib.path as mpath
# import cartopy.feature as cfeature
import cartopy.crs as ccrs
from datetime import datetime
from matplotlib.pyplot import text
import pandas as pd
import cdflib
import numpy as np
# import glob
import requests
import warnings
import matplotlib.ticker as mticker

warnings.filterwarnings("ignore")


ch_energies = [30000., 20400., 13900.,  9450.,  6460.,  4400.,  3000.,  2040.,
               1392.,   949.,   646.,   440.,   300.,   204.,   139.,    95.,
               65.,    44.,    30.]


def main():
    fname = './data/dmsp-f18_ssj_precipitating-electrons-ions_20120101_v1.1.1.cdf'
    ch = 1
    try:
        D = dmsp_reader(fname, channel=ch)
        x, y, egrided_data, igrided_data = dmsp_grid(D, 2)
        dmsp_polar_plot(x, y, egrided_data, igrided_data, ch=ch, savefig=True)
    except:
        file = 'https://cdaweb.gsfc.nasa.gov/pub/data/dmsp/dmspf18/ssj/precipitating-electrons-ions/2012/dmsp-f18_ssj_precipitating-electrons-ions_20120101_v1.1.1.cdf'
        r = requests.get(file, allow_redirects=True)
        open(fname, 'wb').write(r.content)
        D = dmsp_reader(fname, channel=1)
        x, y, egrided_data, igrided_data = dmsp_grid(D, 2)
        dmsp_polar_plot(x, y, egrided_data, igrided_data, ch=ch, savefig=True)


def dmsp_reader(fname, channel=1, hemisphere='North'):
    """_summary_

    Args:
        fname (str)): (location + name) of the DMSP data
        channel (int, optional): energy channel (0 to 18). Defaults to 1.

    Returns:
        array of N(observations) by 4 (magnetic local time, magnetic latitude, electron energy at channel, ion energy at channel)
    """
    cdf_file = cdflib.CDF(fname)
    eflux, iflux, ch_energy, mlt, mlat, mlon, epoch = [cdf_file[var][()] for var in
                                                       ["ELE_DIFF_ENERGY_FLUX",
                                                        "ION_DIFF_ENERGY_FLUX",
                                                        "CHANNEL_ENERGIES",
                                                        "SC_AACGM_LTIME",
                                                        "SC_AACGM_LAT",
                                                        "SC_AACGM_LON",
                                                        "Epoch"]]
    if hemisphere == 'South':
        ind = mlat < -40
        eData = eflux[ind, channel]
        iData = iflux[ind, channel]
        MLON = mlon[ind]
        EPOCH = epoch[ind]
    else:
        ind = mlat > 40
        eData = eflux[ind, channel]
        iData = iflux[ind, channel]
        MLON = mlon[ind]
        EPOCH = epoch[ind]

    return np.vstack([mlt[ind], mlat[ind], eData, iData]).T


def round_off(number, reciprocal=2):
    """Round a number to the closest half integer (reciprocal=2) or quarter (reciprocal=4)... etc
    """
    return round(number * reciprocal) / reciprocal


index_labels = ['MLT', 'MLAT', 'ePrep', 'iPrep']


def dmsp_grid(D, reciprocal):
    """_summary_

    Args:
        D (array)): array of N by M with 
        reciprocal(int): this can take only 1, 2, and 4... corresponding to (1 by 1), (0.5 by 0.5) and [0.25 by 0.25] (mlt, mlat) grid

    Returns:
        x, y grids and the grided data (electron and ion energies)
    """
    df = pd.DataFrame(D, columns=index_labels)
    df = df.apply(round_off, args=[reciprocal])
    f = df.groupby(['MLT', 'MLAT'])['ePrep', 'iPrep'].apply(
        lambda g: g.mean(skipna=True))
    f = f.reset_index()

    lati = np.arange(f['MLAT'].min(), f['MLAT'].max() +
                     (1/reciprocal), (1/reciprocal))
    tt = np.arange(f['MLT'].min(), f['MLT'].max() +
                   (1/reciprocal), (1/reciprocal))
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


def dmsp_polar_plot(x, y, egrided_data, igrided_data, hemisphere='North', ch=1, savefig=False):
    """This function plots ion and electron energy on a specific channel (ch) on a given grid.

    Args:
        x (float): gridded mlt
        y (_type_): gridded magnetic latitude
        egrided_data (_type_): electrons' energy on the grid (from channel ch)
        igrided_data (_type_): _description_
        ch (int, optional): energy channels of DMSP (0 to 18) Defaults to 1.
        savefig (bool, optional): save figure to the current directory. Defaults to False.
    """
    fig = plt.figure(figsize=[13, 5])

    if hemisphere == 'South':
        projection = ccrs.SouthPolarStereo()
    else:
        projection = ccrs.NorthPolarStereo()

    ax1 = fig.add_subplot(1, 2, 1, projection=projection)
    ax2 = fig.add_subplot(1, 2, 2, projection=projection)

    fig.subplots_adjust(bottom=0.05, top=0.95, left=0.04,
                        right=0.95, wspace=0.02)

    center, radius = [0.5, 0.5], 0.5
    theta = np.linspace(0, 2*np.pi, 100)
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    loc_x_mlt = [0.485, 0.86, 1.01, 0.86, 0.485, 0.1, -0.05, 0.1]
    loc_y_mlt = [-0.04, 0.11, 0.485, 0.86, 1.02, 0.86, 0.485, 0.1]
    loc_x_lat = [0.5] * 6
    loc_y_lat = [0.47, 0.4, 0.3, 0.2, 0.1, 0.]

    mlt_label = [str(elem) for elem in np.arange(0, 24, 3)]

    for ax in [ax1, ax2]:
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.set_extent([-180, 180, 40, 90] if hemisphere ==
                      'North' else [-180, 180, -40, -90], ccrs.PlateCarree())

        c = ax.pcolor(x * 15, y, np.log10(egrided_data if ax == ax1 else igrided_data),
                      transform=ccrs.PlateCarree(), cmap='jet', vmin=4, vmax=8)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                          linewidth=1, color='black', alpha=0.3, linestyle='--')

        xx = np.arange(-180, 180, 45)
        gl.xlocator = mticker.FixedLocator(xx)

        for xmlt, ymlt, label_mlt in zip(loc_x_mlt, loc_y_mlt, mlt_label):
            ax.text(xmlt, ymlt, label_mlt, transform=ax.transAxes)

        lat_label = [str(elem) for elem in np.arange(-90, -30, 10)]
        if hemisphere == 'South':
            yticks = list(np.arange(-40, -90, -15))
        else:
            yticks = list(np.arange(40, 90, 15))

        for x_lat, ylat, label_lat in zip(loc_x_lat, loc_y_lat, lat_label):
            ax.text(x_lat, ylat, label_lat, transform=ax.transAxes)

        ax.text(0.85, 0.95, str(
            ch_energies[ch] / 1000) + 'KeV', transform=ax.transAxes)
        # ax.set_title('Log10 (Electron diff. energy flux)' if ax ==
        #              ax1 else 'Log10 (Proton diff. energy flux)')

        fig.colorbar(c, label='Log10 (Proton diff. energy flux)' if ax ==
                     ax2 else 'Log10 (Electron diff. energy flux)')

    if savefig:
        plt.savefig(
            f'DMSP_at_{ch_energies[ch] / 1000}KeV_{hemisphere}.png', dpi=800)

    plt.show()


if __name__ == '__main__':
    main()
