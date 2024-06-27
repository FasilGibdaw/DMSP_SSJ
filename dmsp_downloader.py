import datetime as dt
import pandas as pd
import numpy as np
import os


def read_ssj(event, sat, basepath='./', tempfile_path='./', forcenew=False, **madrigal_kwargs):
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
            Satellite ID for DMSP Block 5D satellite (16-19)
        basepath (str, optional): 
            path to raw files. currently is only for temporary storage of files from madrigal. Defaults to './'.
        tempfile_path (str, optional): 
            Path to dir where processed hdf files are placed. Defaults to './'.
        forcenew (bool, optional): 
            Force the function to download the data even if file exists. Defaults to False.
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
        
    savefile = tempfile_path + event.replace('-', '') + '_ssj_f' + str(sat) + '.h5'

    # do not download if file already exists
    if os.path.isfile(savefile) and not forcenew:
        return savefile

    # imports
    # silence NaturalNameWarnings
    import warnings
    from tables import NaturalNameWarning
    warnings.filterwarnings('ignore', category=NaturalNameWarning)

    # If file does not exist already, we need API to download
    try:
        import madrigalWeb.madrigalWeb  # API, http://cedar.openmadrigal.org/docs/name/rt_contents.html
    except ModuleNotFoundError:
        raise ModuleNotFoundError('read_ssies: Could not import MadrigalWeb module. Will not be able to download DMSP SSIES files.')

    madrigalUrl = 'http://cedar.openmadrigal.org'
    # madrigalUrl = 'http://madrigal.haystack.mit.edu/madrigal'
    try:
        testData = madrigalWeb.madrigalWeb.MadrigalData(madrigalUrl)
    except:
        raise RuntimeError('Madrigal site is not working. Try a manual download.')
        
    # specify one day download
    import calendar
    sTime = dt.datetime(int(event[0:4]), int(event[5:7]), int(event[8:10]), 0, 0)
    eTime = sTime + dt.timedelta(days = 1)
    usTime = calendar.timegm(sTime.utctimetuple())
    ueTime = calendar.timegm(eTime.utctimetuple())

    expList = testData.getExperiments(8100, sTime.year, sTime.month, sTime.day, sTime.hour, sTime.minute, 0, 
                                        eTime.year, eTime.month, eTime.day, eTime.hour, eTime.minute, 0)

    dmsp = pd.DataFrame()
    dmsp2 = pd.DataFrame()   # for the density fraction
    dmsp3 = pd.DataFrame()   # for the electron density

    for i in expList:        # get neccessary files
        fileList = testData.getExperimentFiles(i.id)
        filenames = []

        for fname in fileList:
            filenames.append(fname.name)

        ssj = str([s for s in filenames if '_' + str(sat) + 'e' in s][0])
        ssies = str([s for s in filenames if '_' + str(sat) + 's1' in s][0])
        temp_dens = str([s for s in filenames if '_' + str(sat) + 's4' in s][0])
        
        datafile = basepath + 'ssj_temp_' + event + '.hdf5'
        result = testData.downloadFile(ssj, datafile, **madrigal_kwargs, format = "hdf5")
        f = pd.read_hdf(datafile, mode = 'r', key = 'Data/Table Layout')
        
        tempdensfile = basepath + 'ssj_tempdens_data_' + event + '.hdf5'
        result = testData.downloadFile(temp_dens, tempdensfile, **madrigal_kwargs, format = "hdf5")
        f2 = pd.read_hdf(tempdensfile, mode = 'r', key = 'Data/Table Layout')
        
        tempssiesfile = basepath + 'ssies_temp_data_' + event + '.hdf5'
        result = testData.downloadFile(ssies, tempssiesfile, **madrigal_kwargs, format = "hdf5")
        f3 = pd.read_hdf(tempssiesfile, mode = 'r', key = 'Data/Table Layout')
        
        use = (f.ut1_unix >= usTime) & (f.ut1_unix < ueTime)
        use2 = (f2.ut1_unix >= usTime) & (f2.ut1_unix < ueTime)
        use3 = (f3.ut1_unix >= usTime) & (f3.ut1_unix < ueTime)
        temp = f[use]       # pd.DataFrame()
        temp2 = f2[use2]    # pd.DataFrame()
        temp3 = f3[use3]    # pd.DataFrame()

        dmsp = pd.concat([dmsp, temp])
        dmsp2 = pd.concat([dmsp2, temp2])
        dmsp3 = pd.concat([dmsp3, temp3])
        
    dmsp.index = np.arange(len(dmsp))
    dmsp2.index = np.arange(len(dmsp2))
    dmsp3.index = np.arange(len(dmsp3))

    # set datetime as index
    times = []
    for i in range(len(dmsp)):
        times.append(dt.datetime.utcfromtimestamp(int(dmsp['ut1_unix'][i])))
    dmsp.index = times

    #reindexing due to lower cadence measurements 
    times2 = []
    for i in range(len(dmsp2)):
        times2.append(dt.datetime.utcfromtimestamp(int(dmsp2['ut1_unix'][i])))
    dmsp2.index = times2
    dmsp2 = dmsp2.reindex(index = dmsp.index, method = 'nearest', tolerance = '2sec')

    times3 = []
    for i in range(len(dmsp3)):
        times3.append(dt.datetime.utcfromtimestamp(int(dmsp3['ut1_unix'][i])))
    dmsp3.index = times3
    dmsp3 = dmsp3.reindex(index = dmsp.index, method = 'nearest', tolerance = '2sec')

    dmsp.loc[:,'po+'] = dmsp2['po+']
    dmsp.loc[:,'te'] = dmsp2['te']
    dmsp.loc[:,'ti'] = dmsp2['ti']
    dmsp.loc[:,'ne'] = dmsp3['ne']


    # quality flag from Zhu et al 2020 page 8 https://doi.org/10.1029/2019JA027270
    # flag1 is good, flag2 is also usually acceptable
    flag1 = (dmsp['ne'] > 1e9) & (dmsp['po+'] > 0.85)
    flag2 = ((dmsp['po+'] > 0.85) & (dmsp['ne'] > 1e8) & (dmsp['ne'] < 1e9)) | ((dmsp['po+'] > 0.75) 
                                                                                & (dmsp['po+'] < 0.85) & (dmsp['ne'] > 1e8))
    flag3 = (dmsp['ne'] < 1e8) | (dmsp['po+'] < 0.75) 
    flag4 = dmsp['po+'].isna()

    dmsp.loc[:,'quality'] = np.zeros(dmsp.shape[0]) + np.nan
    dmsp.loc[flag1,'quality'] = 1
    dmsp.loc[flag2,'quality'] = 2
    dmsp.loc[flag3,'quality'] = 3
    dmsp.loc[flag4,'quality'] = 4

    dmsp.to_hdf(savefile, key = 'df', mode = 'w')
    print('DMSP SSJ file saved: ' + savefile)

    # remove temporary files
    os.remove(datafile)
    os.remove(tempdensfile)
    os.remove(tempssiesfile)
    return savefile