from urllib.parse import urljoin
from bs4 import BeautifulSoup
import requests
from lompe.data_tools.dataloader import cross_track_los
from scipy import interpolate
import datetime as dt
import pandas as pd
import numpy as np
import xarray as xr
import os
import warnings
from tables import NaturalNameWarning

warnings.filterwarnings("ignore")


def read_ssj(event, sat='all', basepath='./', tempfile_path='./', **madrigal_kwargs):
    """Download and read DMSP SSJ data from Madrigal database. 
    Example usage:
        event = '2014-12-17'
        username = 'First Last'
        email = 'name@host.com'
        aff = 'University'
        basepath = '/Users/fasilkebede/Documents/LOMPE/Data/DMSP'
        tempfile_path ='/Users/fasilkebede/Documents/LOMPE/Data/DMSP'
        sat = 17 # DMSP satellite number
        madrigal_kwargs = {'user_fullname' : username, 'user_email' : email, 'user_affiliation' : aff}

    Args:
        event (str): 
            format 'YYYY-MM-DD'
        sat (int): 
            Satellite ID for DMSP Block 5D satellite (16-19), if not, all available data from satellites are downloaded. Defaults to all.
        basepath (str, optional): 
            path to raw files. currently is only for temporary storage of files from madrigal. Defaults to './'.
        tempfile_path (str, optional): 
            Path to dir where processed hdf files are placed. Defaults to './'.
        **madrigalkwargs (dict, optional):
            needed to download data from MadrigalWeb (see example usage above)

    Raises:
        ModuleNotFoundError: 
            if madrigalWeb module is not found ( if so install it using 'pip install madrigalWeb')
        RuntimeError: 
            if madrigal site is not working (sometimes is the case, try again later or mannual download)

    Returns:
        str: _description_
    Note:
        Based on dataloader in lompe package: Karl Magnus Laundal
    """
    if not tempfile_path.endswith('/'):
        tempfile_path += '/'
    if not basepath.endswith('/'):
        basepath += '/'
    if sat != 'all':
        savefile = basepath + \
            event.replace('-', '') + '_ssj_f' + str(sat) + '.nc'
        if os.path.isfile(savefile):
            return savefile
    warnings.filterwarnings('ignore', category=NaturalNameWarning)

    try:
        import madrigalWeb.madrigalWeb
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'read_ssj: Could not import MadrigalWeb module. Will not be able to download DMSP SSJ files.')

    madrigalUrl = 'http://cedar.openmadrigal.org'
    testData = madrigalWeb.madrigalWeb.MadrigalData(madrigalUrl)

    sTime = dt.datetime(int(event[0:4]), int(
        event[5:7]), int(event[8:10]), 0, 0)
    eTime = sTime + dt.timedelta(days=0)
    expList = testData.getExperiments(8100, sTime.year, sTime.month, sTime.day, sTime.hour, sTime.minute, 0,
                                      eTime.year, eTime.month, eTime.day, eTime.hour, eTime.minute, 0)

    savefiles = []
    for i in expList:
        fileList = testData.getExperimentFiles(i.id)
        filenames = [fname.name for fname in fileList]
        date_str = event.replace('-', '')

        if sat != 'all':
            ssj = str(
                [s for s in filenames if f'_{date_str}_' + str(sat) + 'e.' in s])
            savefile = tempfile_path + \
                event.replace('-', '') + '_ssj_f' + str(sat) + '.nc'
            if not os.path.isfile(savefile) and len(ssj) > 0:
                result = testData.downloadFile(
                    ssj[0], savefile, **madrigal_kwargs, format="netCDF4")
            savefiles.append(savefile)
        else:
            ssj = [s for s in filenames if 'e.' in s and f'_{date_str}_' in s]
            sats = [ssj[i].split('_')[-1][0:2] for i in range(len(ssj))]
            for idx, sat in enumerate(sats):
                savefile = tempfile_path + \
                    event.replace('-', '') + '_ssj_f' + str(sat) + '.nc'
                if not os.path.isfile(savefile) and len(ssj) > 0:
                    result = testData.downloadFile(
                        ssj[idx], savefile, **madrigal_kwargs, format="netCDF4")
                savefiles.append(savefile)
    return savefiles


def download_dmsp_ssj(event, sat, tempfile_path='./', **madrigal_kwargs):
    """ Download DMSP SSJ flux/energy data for a full day from Madrigal database.
        using FTP-like acess at http://cedar.openmadrigal.org/ftp/

    Args:
        event (str): 
            string on format 'yyyy-mm-dd' to specify time of event for model.
        sat (int): 
            Satellite ID for DMSP Block 5D satellite (16-19)
        tempfile_path (str, optional): 
            Path to dir where processed hdf files are placed. Default: './'
        **madrigalkwargs (dict): 
            needed to download data from MadrigalWeb (Madrigal needs user specifications through **madrigal_kwargs.)
    Example usage:
        event = '2014-07-13'
        sat = 17
        tempfile_path = './'
        madrigal_kwargs = {'user_fullname': 'First','user_email': 'name@host.com', 'user_affiliation': 'University'}

    Returns:
        savefile(str):
            Path to netCDF4-file containing SSJ data.
    """

    savefile = tempfile_path + 'dms_' + \
        event.replace('-', '') + '_' + str(sat) + 'e.001.nc'
    if os.path.exists(savefile):
        print(f'File: {savefile} already exists')
        return

    date_str = event.replace('-', '')
    year = event[:4]
    url_base = "http://cedar.openmadrigal.org"
    url = url_base + \
        f"/ftp/fullname/{madrigal_kwargs['user_fullname']}/email/{madrigal_kwargs['user_email']}/affiliation/{madrigal_kwargs['user_affiliation']}/kinst/8100/year/{year}/"

    response2 = requests.get(url)
    soup2 = BeautifulSoup(response2.content, 'html.parser')
    urls = [link.get('href') for link in soup2.find_all(
        'a', href=True) if '/format/netCDF4/' in link.get('href')]

    # url_ion_drift = url_base + [a['href'] for a in soup2.find_all(
    #     'a', href=True) if 'ion drift' in a.text and f'F{sat}' in a.text][0]
    # url_plasma_temp = url_base + [a['href'] for a in soup2.find_all(
    #     'a', href=True) if 'plasma temp' in a.text and f'F{sat}' in a.text][0]
    url_flux_energy = url_base + [a['href'] for a in soup2.find_all(
        'a', href=True) if 'flux/energy' in a.text and f'F{sat}' in a.text][0]

    # downloading the ion drift file
    flux_hdf5_file_url = url_flux_energy + '/format/netCDF4/'
    response2 = requests.get(flux_hdf5_file_url)
    soup2 = BeautifulSoup(response2.content, 'html.parser')
    url_target_file = url_base + [link.get('href') for link in soup2.find_all(
        'a', href=True) if '/format/netCDF4/fullFilename/' in link.get('href') and date_str in link.get('href')][0]

    response = requests.get(url_target_file, stream=True)
    # tempfile_path = '/Users/fasilkebede/Documents/'
    filename = (tempfile_path +
                url_target_file[-27:-1]).replace('hdf5', 'nc').replace('F', '')
    with open(filename, 'wb') as file:
        for chunk in response.iter_content(chunk_size=65536):
            if chunk:  # Filter out keep-alive new chunks
                file.write(chunk)
    return filename
