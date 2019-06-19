#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 23 18:03:34 2018

@author: Jacob Fausett
"""


import numpy as np
import pandas as pd
import robust as rb
import os, pickle, pdb, socket, inspect, sys
import matplotlib.pyplot as plt
from pyraf import iraf
import glob
from astropy.io import fits
from tqdm import tqdm
from jdcal import gcal2jd, jd2gcal
import datetime

"""
    
To Do:
        A lot: 
            
         
        Build logging routine similar to Thacher obs_summary with relevent info
        (check phot_pipe.py in Thacher/photometry or utilities)
        
        Slim everything down.  A lot of sloppy code that can be made more "Pythonic"

"""


########################################################################################################################
# Get appropriate paths (From Thachers Tabby_reduction)
def get_paths():
    """
    Description
    -----------
    Take account of the user and computer that is executing this script
    then return the appropriate data and outpath.

    Inputs
    ------
    None

    Outputs
    -------
    dpath = (string) path to raw data
    opath = (string) output path
    calpath = (string) path to archive
    execpath = (string) path to directory of this file
    backpath = (string) path to backup directory

    Example
    -------
    dpath,opath,calpath,execpath = get_paths()

    """

    # Environment variables
    user = os.environ['USER']
    home = os.environ['HOME']

    # Host not present on all operating systems
    try:
        host = os.environ['HOST']
    except:
        host = socket.gethostname()

    if host == 'Js-MacBook-Air' or host == 'Js-MacBook-Air.nv.charter.com' and user == 'jakrin':
        dpath = home + '/Tabby_Data/KIC8462852/'
        opath = home + '/Pipeline/pipe_out/'
        calpath = home + '/data/cal_files/'
        npath = home

    elif host == 'Astro-iMac' and user == 'jfausett':
        dpath = home + '/Dropbox/Tabby_Data/KIC8462852/'
        opath = home + '/PycharmProjects/gbo_tabby_pipeline/pipe_out/'
        calpath = home + '/Dropbox/Tabby_Calibration_Files/'
        npath = home + '/Dropbox/JFausett/KIC8462852/'

    else:
        print "Host and/or User don't match expected\nCheck code and add your info to create correct paths"

    # Add your local path here if desired...
    # if host == 'yourhost' and user == 'you':

    execpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/'

    return dpath, opath, calpath, execpath, npath, home


########################################################################################################################
#  Backup data and pipout to backup path using os.system and rsync
def backup(servers=['Dropbox', 'Box\ Sync', 'Google\ Drive']):
    '''

    Description
    -----------
    Uses Mac os built in "rsync" to backup latest data, pipeline, and pipe_out
    to the backup path specified in get_paths()

    Inputs
    ------
    None

    Outputs
    -------
    None

    Dependencies
    ------------
    ability to run rsync bash commands using os.system

    Example
    -------
    backup()

    '''

    dpath, opath, calpath, execpath, npath, home = get_paths()

    for server in servers:

        if server == 'Dropbox':
            print '\nSaving current pipeline to {}\n'.format(server)
            os.system('rsync -arv ' + execpath + ' ' + home + '/' + server + '/Synced_Pipeline')
            continue
        else:
            print '\nBacking up data directory from Dropbox to {}\n' .format(server)
            os.system('rsync -arv ' + dpath + ' ' + home + '/' + server + '/Tabby_Data/KIC8462852/')

            print '\nSaving current pipeline to {}\n' .format(server)
            os.system('rsync -arv ' + execpath + ' ' + home + '/' + server + '/Synced_Pipeline')


########################################################################################################################
#  Find new data and parse into data path
def move_new_files():

    dpath, opath, calpath, execpath, npath, home = get_paths()

    newdates = []

    print "Looking for new Data"

    for root, dirs, files in os.walk(npath):
        for file in files:
            if file.startswith('KIC') and file.endswith(".fts"):
                date = file[11:19]
                if not os.path.exists(dpath + date):
                    os.mkdir(dpath + date + '/')
                filename = npath + file
                os.rename(filename, filename.replace(npath, dpath + date + '/'))
                newdates.append(date)
    newdates = np.array(newdates)
    newdates.sort()
    newdates = np.unique(newdates)
    
    if len(newdates) >0:
        print '\nFound and moved data for these dates: \n{}'.format(newdates)
    else:
        print '\nNo new data found\n'


########################################################################################################################
# Return donedates obsdates and caldates
def get_dates(band, prefix=['IRAF_night', 'KIC']):

    dpath, opath, calpath, execpath, npath, home = get_paths()

    donedates = []
    for root, dirs, files in os.walk(opath):
        for file in files:
            if file.startswith(prefix[0]) and file.endswith(".pck"):
                donedates.append(root.split('/')[-1])
    donedates = np.array(donedates)
    donedates.sort()
    donedates = np.unique(donedates)

    obsdates = []
    baddates = []
    for root, dirs, files in os.walk(dpath.replace('Tabby_Data/KIC8462852', 'rejected_Tabby_data')):
        for file in files:
            if file.startswith(prefix[0]) and file.endswith(".pck"):
                baddates.append(root.split('/')[-1])
    for root, dirs, files in os.walk(dpath):
        for file in files:
            if file.startswith(prefix[1]) and file.endswith(".fts"):
                obsdates.append(root.split('/')[-1])
    baddates = np.array(baddates)
    baddates.sort()
    baddates = np.unique(baddates)
    obsdates = np.array(obsdates)
    obsdates = [x for x in obsdates if x not in baddates]
    obsdates.sort()
    obsdates = np.unique(obsdates)
    

    biasdates = []
    for root, dirs, files in os.walk(calpath + 'bias/'):
        for file in files:
            if file.endswith(".fts"):
                biasdates.append(root.split('/')[-1])
    biasdates = np.array(biasdates)
    biasdates.sort()
    biasdates = np.unique(biasdates)

    darkdates = []
    for root, dirs, files in os.walk(calpath + 'dark/'):
        for file in files:
            if file.endswith(".fts"):
                darkdates.append(root.split('/')[-1])
    darkdates = np.array(darkdates)
    darkdates.sort()
    darkdates = np.unique(darkdates)

    flatdates = []
    for root, dirs, files in os.walk(calpath + 'flat/' + band + '/'):
        for file in files:
            if file.endswith(".fts"):
                flatdates.append(root.split('/')[-1])
    flatdates = np.array(flatdates)
    flatdates.sort()
    flatdates = np.unique(flatdates)

    return donedates, baddates, obsdates, biasdates, darkdates, flatdates


########################################################################################################################
# Determine closest calibrationc file dates based on time from observation
def get_cal_dates(date, band):

    # Bug: Need to adjust dates for different runs
    # Not finding any flat for g data taken on 20180602 because there isn't a flat during that run.
    # Need to determin best way to parse dates for dark and flats

    donedates, baddates, obsdates, biasdates, darkdates, flatdates = get_dates(band=band)

    year, month, day = date[0:4], date[4:6], date[6:]

    JD = gcal2jd(year, month, day)

    biasdates = [gcal2jd(x[0:4], x[4:6], x[6:])[1] for x in biasdates]
    darkdates = [gcal2jd(x[0:4], x[4:6], x[6:])[1] for x in darkdates]
    flatdates = [gcal2jd(x[0:4], x[4:6], x[6:])[1] for x in flatdates]

    if JD[1] <= 58222: # Date = 20180414
        biasdates = [x for x in biasdates if x <= 58222]
        darkdates = [x for x in darkdates if x <= 58222]
        flatdates = [x for x in flatdates if x <= 58222]
    elif JD[1] <= 58314: # Date = 20180715
        biasdates = [x for x in biasdates if x <= 58314 and x >= 58223]
        darkdates = [x for x in darkdates if x <= 58314 and x >= 58223]
        flatdates = [x for x in flatdates if x <= 58314 and x >= 58223]
    elif JD[1] <= 58439: # Date = 20181117
        biasdates = [x for x in biasdates if x <= 58439 and x >= 58310]
        darkdates = [x for x in darkdates if x <= 58439 and x >= 58310]
        flatdates = [x for x in flatdates if x <= 58439 and x >= 58310]
    elif JD[1] >= 58440: # Date = 20181118
        biasdates = [x for x in biasdates if x >= 58440]
        darkdates = [x for x in darkdates if x >= 58440]
        flatdates = [x for x in flatdates if x >= 58440]

    biasdate = min(biasdates, key=lambda x: abs(int(x) - JD[1]))
    darkdate = min(darkdates, key=lambda x: abs(int(x) - JD[1]))
    flatdate = min(flatdates, key=lambda x: abs(int(x) - JD[1]))

    biasdate = [2400000.5, biasdate]
    biasdate = jd2gcal(biasdate[0],biasdate[1])
    year, month, day = biasdate[0], biasdate[1], biasdate[2]
    if month < 10:
        month = '0' + str(month)
    if day < 10:
        day = '0' + str(day)
    biasdate = str(year) + str(month) + str(day)

    darkdate = [2400000.5, darkdate]
    darkdate = jd2gcal(darkdate[0], darkdate[1])
    year, month, day = darkdate[0], darkdate[1], darkdate[2]
    if month < 10:
        month = '0' + str(month)
    if day < 10:
        day = '0' + str(day)
    darkdate = str(year) + str(month) + str(day)

    flatdate = [2400000.5, flatdate]
    flatdate = jd2gcal(flatdate[0], flatdate[1])
    year, month, day = flatdate[0], flatdate[1], flatdate[2]
    if month < 10:
        month = '0' + str(month)
    if day < 10:
        day = '0' + str(day)
    flatdate = str(year) + str(month) + str(day)

    return biasdate, darkdate, flatdate


########################################################################################################################
# Calibration frames (From Thachers Tabby_reduction)
def get_cal_frames(biasdate, darkdate, flatdate, setting=None, readnoise=False, band='r', flatten=False):
    """
    Description
    -----------
    Get the calibration frames for date supplied. Priority is for calibration frames
    from the given night of observing to be used. Else, master cal frames are used
    from the data calpath

    Inputs
    ------
    date = (string) date of observation
    setting = (int) camera setting
    readnoise = (boolean) compute readnoise from bias frames?
    band = (string) photometric band
    flatten = (boolean) use flattened flat fields (use with caution)

    Outputs
    -------
    bias = (float array) bias frame (counts)
    dark = (float array) dark frame (counts/sec)
    flat = (float array) flat field (relative sensitivity)

    Example
    -------
    bias,dark,flat = get_cal_frames('20171108',setting=1,readnoise=True,band='V')

    """

    # Get paths
    dpath, opath, calpath, execpath, npath, home = get_paths()

    # Make master bias from nightly calibrations, else use master in calpath
    biasfiles, bct = get_files(d=calpath + '/bias/' + biasdate + '/', prefix='Bias', tag='1X1', suffix='fts')
    if bct > 0:
        bias = master_bias(biasfiles, readnoise=readnoise, tag='_' + biasdate, outdir=calpath + '/master_cal/')
    if bct == 0 and setting == 1:
        try:
            bias, bh = fits.getdata(calpath + 'master_cal/master_bias_' + biasdate + '.fits', 0, header=True)
            print 'Using master biases'
        except:
            print 'No bias frame!'
            bias = None

    # Make master dark from nightly calibrations, else use master in calpath
    darkfiles, dct = get_files(d=calpath + '/dark/' + darkdate + '/', prefix='Dark', tag='1X1', suffix='fts')
    #    print dct
    if dct > 0:
        if bias is None:
            print ''
            print 'DATE: ' + darkdate
            print 'WARNING: creating dark frame with no bias!!!'
            pdb.set_trace()
        dark = master_dark(darkfiles, bias=bias, tag='_' + darkdate, outdir=calpath + '/master_cal/')
    if dct == 0 and setting == 1:
        try:
            dark, dh = fits.getdata(calpath + 'master_cal/master_dark_' + darkdate + '.fits', 0, header=True)
            print 'Using master dark'
        except:
            print 'No dark frame!'
            dark = None

    #   brought in from tabby_master
    flatfiles, fct = get_files(d=calpath + '/flat/' + band + '/' + flatdate + '/', prefix='Flat', suffix='fts')
    if fct > 0:
        flat = master_flat(flatfiles, bias=bias, dark=dark, tag='_' + flatdate, outdir=calpath + 'master_cal/',
                              band=band)
    if fct == 0:
        try:
            # flat = master_flat(None, bias=bias, dark=dark, tag=date, outdir=opath + 'master_flats/')
            flat, fh = fits.getdata(calpath + 'master_cal/master_flat_' + flatdate + '_' + band + '.fits', 0, header=True)
        except:
            if flatten:
                try:
                    print 'Using master_flat'
                    flat, fh = fits.getdata(calpath + "master_flat_" + band + "_flattened.fits", 0, header=True)
                except:
                    flat = None
            else:
                try:
                    flat, fh = fits.getdata(calpath + "master_flat_" + band + ".fits", 0, header=True)
                except:
                    print "No Flat Frame"
                    flat = None
    return bias, dark, flat


########################################################################################################################
# Get cal_frames and process reduction on all images for a night
def do_phot(date, filters=['g', 'r', 'i', 'z', 'V'], setting=None, source='iraf'):
    '''
        Chooses closest calibration dates to use for reduction
        
        Creates a calibrated image for iraf.phot to use then
        deletes image after done
    
    '''

    dpath, opath, calpath, execpath, npath, home = get_paths()

    if not os.path.exists(opath + date):
        os.makedirs(opath + date)

    for bands in tqdm(filters):

        biasdate, darkdate, flatdate = get_cal_dates(date, band=bands)
        print 'Using Bias, Dark, and Flat from: {} {} {}' .format(biasdate, darkdate, flatdate)
        print biasdate, darkdate, flatdate

        bias, dark, flat = get_cal_frames(biasdate, darkdate, flatdate, band=bands, setting=setting)
        if bands == 'V':
            photfiles = glob.glob(dpath + date + "/KIC*" + bands + ".fts")
        else:
            photfiles = glob.glob(dpath + date + "/KIC*" + bands + "'.fts")
        magfiles = glob.glob(opath + date + '/KIC*mag.1')
        if len(magfiles) > 0:
            print 'Photometry already done for: ' + date
            break
        for image in tqdm(photfiles):
            platesolve = None
            hdul = fits.open(image)
            try:
                platesolve = hdul[0].header['PLTSOLVD']
                hdul.close()
            except:
                print "There is no WCS information for this image"
                #print 'Moving to "rejected_data"'
                hdul.close()
            if not platesolve:
                #if not os.path.exists(dpath + 'rejected_data/no_plate_' + date):
                    #os.mkdir(dpath + 'rejected_data/no_plate_' + date)
                #os.rename(image, image.replace(date, 'rejected_data/no_plate_' + date, 1))
                continue
            outname = image.replace(str(dpath + date), str(opath + date))
            data, header = fits.getdata(image, header=True)
            exp = int(header['EXPTIME'])
            dataout = (data - bias - (dark * exp) / flat)
            calname = image.replace('.fts', '_cal.fts')
            fits.writeto(calname, dataout, header)
            if source == 'iraf':
                iraf.phot(calname, output=outname + ".mag.1")
            os.remove(calname)


########################################################################################################################
#  Read iraf .mag files and extract photometry data
def night_phot(date, write=True, source='iraf'):
    '''
        Read the IRAF phot mag files to get photometry stats
        Creates dataframe for every image, then one for every band
        and combined into one master dictionary of dictionarys for night
    '''

    # Get paths
    dpath, opath, calpath, execpath, npath, home = get_paths()

    # Define empty dictionary and dataframe for data
    nightphot = {}
    imagephot = pd.DataFrame()
    if source == 'iraf':
        # loop through all mag files produced by do_phot assigns values to keys
        photfiles = glob.glob(opath + date + "/KIC*.mag*")
        if len(photfiles) == 0:
            print 'No photometry files for ' + date
            print 'Try running do_phot(' + date + ')'
            return
        n = 0
        for files in photfiles:
            colnames = range(0, 8)
            image = pd.read_table(files, names=colnames, header=None, skiprows=77, delim_whitespace=True)
            name, JD, exposure = str(files).replace(dpath + date + '/', ''), float(image.at[1, 3]), float(
                image.at[1, 0])
            airmass, band, tabbyflux = float(image.at[1, 1]), image.at[1, 2], float(image.at[2, 3])
            ref1, ref2, ref3 = float(image.at[7, 3]), float(image.at[12, 3]), float(image.at[17, 3])
            tabbymag = float(image.at[2, 4])
            if tabbymag > 14.5:
                if not os.path.exists(dpath.replace('Tabby_Data/KIC8462852/', 'rejected_data/' + date)):
                    os.mkdir(dpath.replace('Tabby_Data/KIC8462852/', 'rejected_data/' + date))
                os.rename(files, files.replace(opath, dpath.replace('Tabby_Data/KIC8462852/', 'rejected_data/' + date + '/')))
                continue
            imagephot.at[n, 'name'], imagephot.at[n, 'JD'], imagephot.at[n, 'exposure'] = name.replace(
                opath + date + '/', ''), JD, exposure
            imagephot.at[n, 'airmass'], imagephot.at[n, 'filter'], imagephot.at[
                n, 'tabby'] = airmass, band, tabbyflux / exposure
            imagephot.at[n, 'ref1'], imagephot.at[n, 'ref2'], imagephot.at[
                n, 'ref3'] = ref1 / exposure, ref2 / exposure, ref3 / exposure
            imagephot.at[n, 'tabbymag'], imagephot.at[n, 'tabbyerr'] = tabbymag, float(image.at[2, 5])
            imagephot.at[n, 'tabbybkg'], imagephot.at[n, 'tabbyrms'] = float(image.at[0, 0]), float(image.at[0, 2])
            imagephot.at[n, 'tabbyflag'], imagephot.at[n, 'ref1mag'] = str(image.at[2, 7]), float(image.at[7, 4])
            imagephot.at[n, 'ref1magerr'], imagephot.at[n, 'ref1bkg'] = float(image.at[7, 5]), float(image.at[5, 0])
            imagephot.at[n, 'ref1rms'], imagephot.at[n, 'ref1flag'] = float(image.at[5, 2]), str(image.at[7, 7])
            imagephot.at[n, 'ref2magerr'], imagephot.at[n, 'ref2bkg'] = float(image.at[12, 5]), float(image.at[10, 0])
            imagephot.at[n, 'ref2rms'], imagephot.at[n, 'ref2flag'] = float(image.at[10, 2]), str(image.at[12, 7])
            imagephot.at[n, 'ref3magerr'], imagephot.at[n, 'ref3bkg'] = float(image.at[17, 5]), float(image.at[15, 0])
            imagephot.at[n, 'ref3rms'], imagephot.at[n, 'ref3flag'] = float(image.at[15, 2]), str(image.at[17, 7])
            imagephot.at[n, 'ref2mag'], imagephot.at[n, 'ref3mag'] = float(image.at[12, 4]), float(image.at[17, 4])
            n += 1
        # Loop through all bands and create individual dictionary for each
        filters = imagephot['filter'].unique()
        filters = [x.replace("'", "") for x in filters]
        for bands in filters:
            if bands == 'V':
                band = pd.DataFrame(imagephot.where(imagephot['filter'] == bands))
            else:
                band = pd.DataFrame(imagephot.where(imagephot['filter'] == bands + "'"))
            band.dropna(how='all', inplace=True)
            band = band.reset_index(drop=True)
            band.loc[:, 'relflux'] = band.loc[:, 'tabby'] / ((band.loc[:, 'ref1'] + \
                                                              band.loc[:, 'ref2'] + \
                                                              band.loc[:, 'ref3']) / 3)
            band.loc[:, 'ref1rel'] = band.loc[:, 'ref1'] / ((band.loc[:, 'ref2'] + band.loc[:, 'ref3']) / 2)
            band.loc[:, 'ref2rel'] = band.loc[:, 'ref2'] / ((band.loc[:, 'ref1'] + band.loc[:, 'ref3']) / 2)
            band.loc[:, 'ref3rel'] = band.loc[:, 'ref3'] / ((band.loc[:, 'ref2'] + band.loc[:, 'ref1']) / 2)
            band.at[0, 'median_JD'] = band['JD'].median(axis=0)
            band.at[0, 'mean'] = band['relflux'].mean(axis=0)
            band.at[0, 'error'] = band['relflux'].std(axis=0)
            sd = band.at[0, 'error']
            flux = [x for x in band.loc[:, 'relflux']]
            sigflux = [x for x in flux if x < (band.at[0, 'mean'] + 1.5 * sd)]
            sigflux = [x for x in sigflux if x > (band.at[0, 'mean'] - 1.5 * sd)]
            sigflux = np.array(sigflux)
            band.at[0, 'sigmean'] = sigflux.mean()
            band.at[0, 'sigerr'] = sigflux.std()
            band.at[0, 'tabbymagmean'] = band['tabbymag'].mean(axis=0)
            band.at[0, 'tabbymagerr'] = band['tabbymag'].std(axis=0)
            band.at[0, 'tabbyrmsmean'] = band['tabbyrms'].mean(axis=0)
            band.at[0, 'ref1mean'] = band['ref1rel'].mean(axis=0)
            band.at[0, 'ref2mean'] = band['ref2rel'].mean(axis=0)
            band.at[0, 'ref3mean'] = band['ref3rel'].mean(axis=0)
            band.at[0, 'ref1err'] = band['ref1rel'].std(axis=0)
            band.at[0, 'ref2err'] = band['ref2rel'].std(axis=0)
            band.at[0, 'ref3err'] = band['ref3rel'].std(axis=0)
            band.at[0, 'ref1magmean'] = band['ref1mag'].mean(axis=0)
            band.at[0, 'ref1magmeanerr'] = band['ref1mag'].std(axis=0)
            band.at[0, 'ref2magmean'] = band['ref2mag'].mean(axis=0)
            band.at[0, 'ref2magmeanerr'] = band['ref2mag'].std(axis=0)
            band.at[0, 'ref3magmean'] = band['ref3mag'].mean(axis=0)
            band.at[0, 'ref3magmeanerr'] = band['ref3mag'].std(axis=0)
            band.at[0, 'ref1rmsmean'] = band['ref1rms'].mean(axis=0)
            band.at[0, 'ref2rmsmean'] = band['ref2rms'].mean(axis=0)
            band.at[0, 'ref3rmsmean'] = band['ref3rms'].mean(axis=0)
            print band.at[0, 'mean']
            print band.at[0, 'sigmean']
            print band.at[0, 'error']
            print band.at[0, 'sigerr']
            nightphot[bands] = band.to_dict()

        # Write band photometry dictionaries to date file in pipe_out directory
        if write:
            print 'Writing out photometry files into: '
            print '  ' + opath + date + '/'
            fname = opath + date + '/IRAF_night_phot_' + date + '.pck'
            if not os.path.exists(opath + date):
                os.makedirs(opath + date)
            fout = open(fname, "w")
            pickle.dump(nightphot, fout)
            fout.close()


#    return  nightphot

########################################################################################################################
# Add dates to obs_summary and IRAF_all_photometry
def add_new_dates():
    dpath, opath, calpath, execpath, npath, home = get_paths()

    # Outputs an array of dates with images    
    donedates, baddates, obsdates, biasdates, darkdates, flatdates = get_dates(band='r')

    allphot = {}
    dates = []
    # For every date in obs, creates dictionary entry named date, into new allphot dictionary
    i = 0
    for i in tqdm(range(len(donedates))):
        date = donedates.item(i)
        date = str(date)
        dates.append(date)
        fname = open(opath + date + '/IRAF_night_phot_' + date + '.pck', 'r')
        allphot[date] = pickle.load(fname)
        fname.close()
    # Outputs the master dictionary for all photometry
    fname = open(opath + 'IRAF_all_phot.pck', 'w')
    pickle.dump(allphot, fname)
    fname.close()
    print 'There are ' + str(len(donedates)) + ' days of data'


########################################################################################################################
def make_summary():
    dpath, opath, calpath, execpath, npath, home = get_paths()
    donedates, baddates, obsdates, biasdates, darkdates, flatdates = get_dates(band='r')
    filters = ['g', 'r', 'i', 'z', 'V']

    allphot = pd.read_pickle(opath + 'IRAF_all_phot.pck')
    allphot = pd.DataFrame.from_dict(allphot)
    for bands in tqdm(filters):
        band = pd.DataFrame()
        for date in donedates:
            try:
                phot = band.from_dict(allphot[date][bands])
                night = phot.filter(['median_JD', 'sigmean', 'sigerr',
                                     'mean', 'error', 'tabbymagmean', 'tabbymagerr', 'tabbyrmsmean',
                                     'ref1mean', 'ref1err', 'ref1magmean', 'ref1magmeanerr', 'ref1rmsmean',
                                     'ref2mean', 'ref2err', 'ref2magmean', 'ref2magmeanerr', 'ref2rmsmean',
                                     'ref3mean', 'ref3err', 'ref3magmean', 'ref3magmeanerr', 'ref3rmsmean',
                                     ], axis=1)
                night.dropna(how='all', inplace=True)
                band = band.append(night, ignore_index=True)
            except:
                continue
        band.to_csv(opath + bands + '_summary.csv', index=False)

#######################################################################################################################


# Make plots
def make_plot(date='all', filters=['g', 'r', 'i', 'V'], sigma=False,
              write=False, sub=False, night=False, zoom=False,
              reference=False):
    dpath, opath, calpath, execpath, npath, home = get_paths()

    simJD = 2457800

    if date != 'all':
        plt.figure(figsize=(8, 4))
        plt.xlabel('JD 2457800')
        plt.ylabel('Relative Flux')
        plt.title(date)
        plt.minorticks_on()
        plt.tick_params(which='both', direction='in', right=True)

        if night:
            #filters.append('z')
            for bands in filters:
                if bands == 'g':
                    norm = 0.86353
                    color = 'blue'
                elif bands == 'r':
                    norm = 0.547102
                    color = 'green'
                elif bands == 'z':
                    norm = 0.3688151
                    color = 'red'
                band = pd.read_pickle(opath + date + '/IRAF_night_phot_' + date + '.pck')
                if date == '20180326':
                    plt.title('Evangeline Egress: Day 1')
                elif date == '20180327':
                    plt.title('Evangeline Egress: Day 2')
                elif date == '20180328':
                    plt.title('Evangeline Egress: Day 3')
                elif date == '20180330':
                    plt.title('Evangeline Egress: Day 5')
                try:
                    plt.ylim(.89, 1.01)
                    plt.axhline(y=1, linestyle='--', linewidth=1.5, color='black')
                    band = pd.DataFrame.from_dict(band[bands])
                    flux = band.loc[:, 'relflux']
                    y = flux/norm
                    plt.scatter(band.loc[:, 'JD'] - simJD, y, color = color, label = bands+"'")
                except:
                    continue


        else:
            plt.ylim(.94, 1.04)
            plt.axhline(y=1, linestyle='--')
            for bands in filters:
                if bands == 'g':
                    band = pd.read_pickle(opath + date + '/IRAF_night_phot_' + date + '.pck')
                    band = band['g']
                    plt.errorbar(band['median_JD'][0] - simJD,
                                 band['mean'][0] / (0.2893),
                                 yerr=band['error'][0] / (0.2893),
                                 fmt='o', label="g'", color='blue')
                elif bands == 'r':
                    band = pd.read_pickle(opath + date + '/IRAF_night_phot_' + date + '.pck')
                    band = band['r']
                    plt.errorbar(band['median_JD'][0] - simJD,
                                 band['mean'][0] / (0.1832),
                                 yerr=band['error'][0] / (0.1832),
                                 fmt='o', label="r'", color='green')
                elif bands == 'i':
                    band = pd.read_pickle(opath + date + '/IRAF_night_phot_' + date + '.pck')
                    band = band['i']
                    plt.errorbar(band['median_JD'][0] - simJD,
                                 band['mean'][0] / (0.1462),
                                 yerr=band['error'][0] / (0.1462),
                                 fmt='o', label="i'", color='orange')
                elif bands == 'z':
                    band = pd.read_pickle(opath + date + '/IRAF_night_phot_' + date + '.pck')
                    band = band['z']
                    plt.errorbar(band['median_JD'][0] - simJD,
                                 band['mean'][0] / (0.1232),
                                 yerr=band['error'][0] / (0.1232),
                                 fmt='o', label="z'", color='red')
                elif bands == 'V':
                    band = pd.read_pickle(opath + date + '/IRAF_night_phot_' + date + '.pck')
                    band = band['V']
                    plt.errorbar(band['median_JD'][0] - simJD,
                                 band['mean'][0] / (0.2309),
                                 yerr=band['error'][0] / (0.2309),
                                 fmt='o', label="V", color='limegreen')
        plt.legend(title='Filters', loc=4)
        if write:
            plt.savefig(opath + '/plots/tabbyflux_' + date + '.png', dpi=300)

    if date == 'all':

        if zoom != 'evangeline':
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig, ax = plt.subplots()
        ax.set_xlabel('Time (JD-2457800)', fontsize=14)
        ax.set_ylabel('Relative Flux', fontsize=14)
        ax.set_title('GBO Observations of KIC 8462852 \n 2018 March - current', fontsize=16)
        #        ax.set_xlim(378,670)
        if 'g' in filters:
            ax.set_ylim(.92, 1.04)
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))
        else:
            ax.set_ylim(.94, 1.04)
            ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.5)
        ax.axhline(y=1, linestyle='--', linewidth=1.5, color='black')
        ax.minorticks_on()
        ax.tick_params(which='major', labelsize=12, length=6, width=1.5, direction='in', right=True, top=True)
        ax.tick_params(which='minor', length=3, width=1.5, direction='in', right=True, top=True)
        if zoom != 'evangeline':
            ax.annotate('Evangeline', xy=(404, 1.002), xytext=(404, 1.025),
                        arrowprops=dict(facecolor='black', headwidth=6, headlength=5, width=1))
        fig.tight_layout(pad=2)

        if zoom == '2018':
            ax.set_title('GBO Observations of KIC 8462852 \n 2018 March - December')
            ax.set_xlim(390, 680)
            ax.set_ylim(.92,1.04)
        elif zoom == 'evangeline':
            ax.set_title('Evangeline')
            ax.set_xlim(400, 415)
            ax.set_ylim(0.91,1.01)
        elif zoom == False:
            ax.set_title('GBO Observations of KIC 8462852 \n 2018 March - current')
            ax.set_xlim(370, 840)

        nevadablue = '#003366'
        sfred = '#AA0000'
        gcolor = '#00c8ff'
        rcolor = '#ff7500'
        icolor = '#860000'
        Vcolor = '#99ff00'
        zcolor = '#610000'

        gcolor = 'blue'
        #rcolor = 'orange'
        #rcolor = 'green'
        # icolor = 'purple'
        Vcolor = 'green'

        for bands in filters:
            band = pd.read_csv(opath + bands + '_summary.csv')
            if bands == 'g':
                norm = band.loc[35:39, 'mean']
                norm = np.mean(norm)
                print 'g norm is {}'.format(norm)
                if sigma:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'sigmean'] / norm,
                                yerr=band.loc[:, 'sigerr'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=gcolor, label="g'")
                else:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'mean'] / norm,
                                yerr=band.loc[:, 'error'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=gcolor, label="g'")

            elif bands == 'r':
                norm = band.loc[35:39, 'mean']
                norm = np.mean(norm)
                print 'r norm is {}'.format(norm)
                if sigma:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'sigmean'] / norm,
                                yerr=band.loc[:, 'sigerr'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=rcolor, label="r'")
                else:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'mean'] / norm,
                                yerr=band.loc[:, 'error'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=rcolor, label="r'")

            elif bands == 'i':
                norm = band.loc[11:15, 'mean']
                norm = np.mean(norm)
                print 'i norm is {}'.format(norm)
                if sigma:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'sigmean'] / norm,
                                yerr=band.loc[:, 'sigerr'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=icolor, label="i'")
                else:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'mean'] / norm,
                                yerr=band.loc[:, 'error'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=icolor, label="i'")
            elif bands == 'z':
                norm = band.loc[34:38, 'mean']
                norm = np.mean(norm)
                print 'z norm is {}'.format(norm)
                if sigma:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'sigmean'] / norm,
                                yerr=band.loc[:, 'sigerr'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=zcolor, label="z'")
                else:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'mean'] / norm,
                                yerr=band.loc[:, 'error'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=zcolor, label="z'")
            elif bands == 'V':
                norm = band.loc[35:39, 'mean']
                norm = np.mean(norm)
                print 'V norm is {}'.format(norm)
                if sigma:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'sigmean'] / norm,
                                yerr=band.loc[:, 'sigerr'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=Vcolor, label="V")
                else:
                    ax.errorbar(band.loc[:, 'median_JD'] - simJD, band.loc[:, 'mean'] / norm,
                                yerr=band.loc[:, 'error'] / norm,
                                markersize=4, linewidth=1, fmt='o', color=Vcolor, label="V")
        ax.legend(title='Filters', loc=4)

        today = datetime.datetime.now().date()

        if write:
            if len(filters) > 1:
                if zoom == '2018' and filters.__contains__('z'):
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_all_2018_z.png', dpi=300)
                elif zoom == '2018':
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_all_2018.png', dpi=300)
                elif zoom == 'evangeline' and filters.__contains__('z'):
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_all_evangeline_z.png', dpi=300)
                elif zoom == 'evangeline':
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_all_evangeline.png', dpi=300)
                elif filters.__contains__('z'):
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_all_z.png', dpi=300)
                elif zoom == False:
                    plt.savefig(opath + 'plots/GBO_{}_Update.png'.format(today), dpi=300)
                if sigma:
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_all_sigma.png', dpi=300)

            else:
                if filters == 'g' and sigma:
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_g_sigma.png', dpi=300)
                if filters == 'g':
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_g.png', dpi=300)
                if filters == 'r' and sigma:
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_r_sigma.png', dpi=300)
                if filters == 'r':
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_r.png', dpi=300)
                if filters == 'i' and sigma:
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_i_sigma.png', dpi=300)
                if filters == 'i':
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_i.png', dpi=300)
                if filters == 'z' and sigma:
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_z_sigma.png', dpi=300)
                if filters == 'z':
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_z.png', dpi=300)
                if filters == 'V' and sigma:
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_V_sigma.png', dpi=300)
                if filters == 'V':
                    plt.savefig(opath + 'plots/tabbyflux_IRAF_V.png', dpi=300)

########################################################################################################################
def get_files(d="/",prefix='',tag='',suffix='.fits',clean=True):

    """
    Overview:
    ---------
    Returns list of files with a user defined prefix and suffix withing a
    specified directory


    Calling sequence:
    -----------------
    files = get_files('HATp33b',d='/home/users/bob/stuff/')

    """

    files = glob.glob(d+prefix+"*"+tag+"*"+suffix)

    fct = len(files)

    # Work around for difficult single quote and inconsistent file naming convention
    # due to filter names
    if clean:
        for file in files:
            inname  = file.replace("'","\\'")
            outname =  file.replace("'","")
            if inname != outname:
                mvcmd = "mv "+inname+" "+outname
                os.system(mvcmd)

        files = [file.replace("'","") for file in files]

        for file in files:
            inname  = file
            outname =  file.replace("p.fts",".fts")
            if inname != outname:
                mvcmd = "mv "+inname+" "+outname
                os.system(mvcmd)

        files = [file.replace("p.fts",".fts") for file in files]

    return files,fct



########################################################################################################################
def get_summary():
    """
    Description
    -----------
    Return the contents of the obs_summary file that lives in the tabby repo

    Inputs
    ------
    None

    Output
    ------
    Pandas dataframe with obs_summary information

    Example
    -------
    obs = get_summary()

    """

    # Path to this file
    dpath, opath, calpath, execpath, npath, home = get_paths()

    # Read obs_summary file in this directory
    obs = pd.read_csv(execpath + 'IRAF_obs_summary.csv')

    return obs


# ----------------------------------------------------------------------#
# master_bias:
# ----------------------------------------------------------------------#

def master_bias(files, write=True, outdir='/', readnoise=False, clobber=False, verbose=True,
                float32=True, tag='', median=True):
    """
   Overview:
    ---------
    Create master bias frame from series of biases (median filter).
    Returns a master_bias frame and writes FITS file to disk in specified
    directory.

    Optionally, the read noise is calculated from the variance of each
    pixel in the bias stack. This is *very* slow. So only use this option
    if you really need to. The readnoise image is also written to disk.

    Inputs:
    -------
    files       : List of flat field files from which a master bias will be created.
                  Must be provided, no default.

    Keyword arguments:
    ------------------
    write       : Toggle to write files to disk (default True)
    outdir      : Directory to which output files are written (default pwd)
    clobber     : Toggle to overwrite files if they already exist in outdir
                  (default False)
    readnoise   : Do readnoise calculation (very slow! default False)
    verbose     : Print out progress (default True)

    Calling sequence:
    -----------------
    master_bias = master_bias(biasfiles,write=True,readnoise=False,
                              outdir='/home/users/bob/stuff/')


    """

    # Don't redo master_bias unless clobber keyword set
    name = outdir + 'master_bias' + tag + '.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master bias already exists!")
        master_bias = fits.getdata(name, 0, header=False)
        return master_bias

    # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    ysz, xsz = image.shape
    stack = np.zeros((fct, ysz, xsz))
    temps = []

    # Load stack array and get CCD temperatures
    for i in np.arange(fct):
        output = 'Reading {}: frame {} of {} \r'.format(files[i].split('/')[-1], \
                                                        str(i + 1), str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        temps.append(header["CCD-TEMP"])
        stack[i, :, :] = image

    # Calculate read noise directly from bias frames if prompted
    if readnoise:
        rn = np.zeros((ysz, xsz))
        print("Starting readnoise calculation")
        pbar = tqdm(desc='Calculating readnoise', total=ysz, unit='rows')
        for i in np.arange(ysz):
            for j in np.arange(xsz):
                rn[i, j] = rb.std(stack[:, i, j])
            pbar.update(1)

        # Make a nice plot (after all that hard work)
        aspect = np.float(xsz) / np.float(ysz)
        plt.figure(39, figsize=(5 * aspect * 1.2, 5))
        plt.clf()
        sig = rb.std(rn)
        med = np.median(rn)
        mean = np.mean(rn)
        vmin = med - 2 * sig
        vmax = med + 2 * sig
        plt.imshow(rn, vmin=vmin, vmax=vmax, cmap='gist_heat', interpolation='nearest', origin='lower')
        plt.colorbar()
        plt.annotate(r'$\bar{\sigma}$ = %.2f cts' % mean, [0.95, 0.87], horizontalalignment='right',
                     xycoords='axes fraction', fontsize='large')
        #                    path_effects=[PathEffects.SimpleLineShadow(linewidth=3,foreground="w")])
        plt.annotate(r'med($\sigma$) = %.2f cts' % med, [0.95, 0.8], horizontalalignment='right',
                     xycoords='axes fraction', fontsize='large')
        #                    path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
        plt.annotate(r'$\sigma_\sigma$ = %.2f cts' % sig,
                     [0.95, 0.73], horizontalalignment='right',
                     xycoords='axes fraction', fontsize='large')
        #                    path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
        plt.title("Read Noise")
        plt.xlabel("pixel number")
        plt.ylabel("pixel number")

        if write:
            plt.savefig(outdir + 'readnoise' + tag + '.png', dpi=300)

    # Calculate master bias frame by median filter
    print('Calculating median of stacked frames...')
    if median:
        master_bias = np.median(stack, axis=0)
    else:
        master_bias = np.mean(stack, axis=0)

    # Make a plot
    aspect = np.float(xsz) / np.float(ysz)
    plt.figure(38, figsize=(5 * aspect * 1.2, 5))
    plt.clf()
    sig = rb.std(master_bias)
    med = np.median(master_bias)
    vmin = med - 2 * sig
    vmax = med + 2 * sig
    plt.imshow(master_bias, vmin=vmin, vmax=vmax, cmap='gist_heat', interpolation='nearest', origin='lower')
    plt.colorbar()
    plt.annotate('Bias Level = %.2f cts' % med, [0.95, 0.87], horizontalalignment='right',
                 xycoords='axes fraction', fontsize='large', color='k')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\sigma$ = %.2f cts' % sig, [0.95, 0.8], horizontalalignment='right',
                 xycoords='axes fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\langle T_{\rm CCD} \rangle$ = %.2f C' % np.median(temps),
                 [0.95, 0.73], horizontalalignment='right',
                 xycoords='axes fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.title("Master Bias")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

    # Write out bias, readnoise and plot
    if write:
        name = outdir + 'master_bias' + tag
        plt.savefig(name + '.png', dpi=300)

        hout = fits.Header()
        hout['CCDTEMP'] = (np.median(temps), "Median CCD temperature")
        hout["TEMPSIG"] = (np.std(temps), "CCD temperature RMS")
        hout["BIAS"] = (med, "Median bias level (cts)")
        hout["BIASSIG"] = (sig, "Bias RMS (cts)")
        if len(glob.glob(name + '.fits')) == 1:
            os.system('rm ' + name + '.fits')
        if float32:
            fits.writeto(name + '.fits', np.float32(master_bias), hout)
        else:
            fits.writeto(name + '.fits', master_bias, hout)

        if readnoise:
            name = outdir + 'readnoise' + tag
            if len(glob.glob(name + '.fits')) == 1:
                os.system('rm ' + name + '.fits')
            if float32:
                fits.writeto(name + '.fits', np.float32(rn), hout)
            else:
                fits.writeto(name + '.fits', rn, hout)

    return master_bias


# ----------------------------------------------------------------------#
# master_dark:
# ----------------------------------------------------------------------#

def master_dark(files, bias=None, write=True, outdir='/', clobber=False, float32=True, tag='',
                median=True):
    """
    Overview:
    ---------
    Create master dark frame from series of darks (median filter).
    Returns a master dark frame. If write is specified, a FITS file
    will be written to "outdir" (default is pwd).

    Inputs:
    -------
    files       : List of flat field files from which a master dark will be created.
                  Must be provided, no default.

    Keyword arguments:
    ------------------
    bias        : Master bias frame (default None)
    write       : Toggle to write files to disk (default True)
    outdir      : Directory to which output files are written (default pwd)
    clobber     : Toggle to overwrite files if they already exist in outdir
                  (default False)

    Calling sequence:
    -----------------
    master_dark = master_dark(darkfiles,bias=master_bias,write=True,
                              outdir='/home/users/bob/stuff/')

    """

    # Don't redo master_dark unless clobber keyword set
    name = outdir + 'master_dark' + tag + '.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master dark already exists!")
        master_dark = fits.getdata(name, 0, header=False)
        return master_dark

    # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    ysz, xsz = image.shape
    stack = np.zeros((fct, ysz, xsz))
    temps = []
    exps = []

    # Load stack array and get CCD temperatures
    for i in np.arange(fct):
        output = 'Reading {}: frame {} of {} \r'.format(files[i].split('/')[-1], \
                                                        str(i + 1), str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        exp = header["EXPOSURE"]
        exps.append(exp)
        temps.append(header["CCD-TEMP"])
        if length(bias) == 1:
            image = np.float(image) / exp
        else:
            image = (image - bias) / exp
        stack[i, :, :] = image

    # Obtain statistics for the master dark image header
    # Temperature
    tmax = np.max(temps)
    tmin = np.min(temps)
    tmean = np.mean(temps)
    tmed = np.median(temps)
    tsig = np.std(temps)
    # Exposure times
    expmax = np.max(exps)
    expmin = np.min(exps)
    print('')
    print("Minimum CCD Temp. %.2f C" % tmin)
    print("Maximum CCD Temp. %.2f C" % tmax)
    print("CCD Temp. rms: %.3f C" % tsig)
    print("CCD Temp. mean: %.2f C" % tmean)
    print("CCD Temp. median: %.2f C" % tmed)

    # Create master dark by median filter or mean
    if median:
        master_dark = np.median(stack, axis=0)
    else:
        master_dark = np.mean(stack, axis=0)

    # Make a plot
    sig = rb.std(master_dark)
    med = np.median(master_dark)
    vmin = med - 2 * sig
    vmax = med + 2 * sig
    aspect = np.float(xsz) / np.float(ysz)
    plt.figure(37, figsize=(5 * aspect * 1.2, 5))
    plt.clf()
    plt.imshow(master_dark, vmin=vmin, vmax=vmax, cmap='gist_heat', interpolation='nearest', origin='lower')
    plt.colorbar()
    plt.annotate('Dark Current = %.2f cts/sec' % med, [0.72, 0.8], horizontalalignment='right',
                 xycoords='figure fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\sigma$ = %.2f cts/sec' % sig, [0.72, 0.75], horizontalalignment='right',
                 xycoords='figure fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\langle T_{\rm CCD} \rangle$ = %.2f C' % np.median(temps),
                 [0.72, 0.7], horizontalalignment='right',
                 xycoords='figure fraction', fontsize='large')
    #                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.title("Master Dark")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

    # Write out plot and master dark array
    if write:
        name = outdir + 'master_dark' + tag

        plt.savefig(name + '.png', dpi=300)

        hout = fits.Header()
        hout["TEMPMAX"] = (tmax, "Maximum CCD temperature")
        hout["TEMPMIN"] = (tmin, "Minimum CCD temperature")
        hout["TEMPMED"] = (tmed, "Median CCD temperature")
        hout["TEMPMN"] = (tmean, "Mean CCD temperature")
        hout["TEMPSIG"] = (tsig, "CCD temperature RMS")
        hout["EXPMAX"] = (expmax, "Maximum exposure time")
        hout["EXPMIN"] = (expmin, "Minimum exposure time")
        hout["DARKCNT"] = (med, "Median dark current (cts/sec)")
        hout["DARKSIG"] = (sig, "Dark current RMS (cts/sec)")
        if len(glob.glob(name)) == 1:
            os.system('rm ' + name + '.fits')
        if float32:
            fits.writeto(name + '.fits', np.float32(master_dark), hout)
        else:
            fits.writeto(name + '.fits', master_dark, hout)

    return master_dark


# ----------------------------------------------------------------------#
# master_flat:
# ----------------------------------------------------------------------#

def master_flat(files, bias=None, dark=None, write=True, outdir='/', band='V',
                tag='', clobber=False, stretch=3, float32=True, median=True):
    """
    Overview:
    ---------
    Create a master flat using (optionally) a provided bias and dark frame. Output
    is written to "outdir" in FITS format.

    Inputs:
    -------
    files       : List of flat field files from which a master flat will be created.
                  Must be provided, no default.

    Keyword arguments:
    ------------------
    bias        : Master bias frame (default None)
    dark        : Master dark frame calibrated in ADU/sec (default None)
    band        : Band from which flatis being produced
    write       : Toggle to write files to disk (default True)
    outdir      : Directory to which output files are written (default pwd)
    clobber     : Toggle to overwrite files if they already exist in outdir
                  (default False)
    stretch     : Multiple of the noise RMS to stretch image (default 3)


    Calling sequence:
    -----------------
    master_flat = master_flat(flatfiles,bias=master_bias,dark=master_dark,write=True,
                              outdir='/home/users/bob/stuff/')

    """

    # Don't redo master_dark unless clobber keyword set
    name = outdir + 'master_flat' + tag + '_' + band + '.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master flat already exists!")
        master_flat = fits.getdata(name, 0, header=False)
        return master_flat

    # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    filter = header["filter"]

    ysz, xsz = image.shape
    stack = np.zeros((fct, ysz, xsz))

    # Load stack array and get CCD temperatures
    meds = []
    for i in np.arange(fct):
        output = 'Reading {}: frame {} of {} \r'.format(files[i].split('/')[-1], \
                                                        str(i + 1), str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        image = np.float32(image)
        if header["filter"] != filter:
            sys.exit("Filters do not match!")
        if length(bias) > 1:
            image -= bias
        if length(dark) > 1:
            exptime = header['EXPTIME']
            image -= dark * exptime
        meds.append(np.median(image))
        stack[i, :, :] = image / np.median(image)

    # Obtain statistics for the master dark image header
    med = np.median(meds)
    sig = np.std(meds)

    # Create master flat by median filter
    if median:
        master_flat = np.median(stack, axis=0)
    else:
        master_flat = np.mean(stack, axis=0)

    # Make a plot
    sig = rb.std(master_flat)
    med = np.median(master_flat)
    vmin = med - stretch * sig
    vmax = med + stretch * sig
    aspect = np.float(xsz) / np.float(ysz)
    plt.figure(40, figsize=(5 * aspect * 1.2, 5))
    plt.clf()
    plt.imshow(master_flat, vmin=vmin, vmax=vmax, cmap='gist_heat', interpolation='nearest', origin='lower')
    plt.colorbar()
    plt.title("Master Flat")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

    # Write out plot and master flat array
    if write:
        plt.savefig(outdir + 'master_flat' + tag + '_' + band + '.png', dpi=300)
        hout = fits.Header()
        hout["FILTER"] = (filter, "Filter used when taking image")
        hout["MEDCTS"] = (med, "Median counts in individual flat frames")
        hout["MEDSIG"] = (sig, "Median count RMS in individual flat frames")
        if length(bias) > 1:
            hout.add_comment("Bias subtracted")
        if length(dark) > 1:
            hout.add_comment("Dark subtracted")

        if len(glob.glob(outdir + 'master_flat' + tag + '.fits')) == 1:
            os.system('rm ' + outdir + 'master_flat' + tag + '.fits')
        if float32:
            fits.writeto(outdir + 'master_flat' + tag + '_' + band + '.fits', np.float32(master_flat), hout)
        else:
            fits.writeto(outdir + 'master_flat' + tag + '_' + band + '.fits', master_flat, hout)

    return master_flat

########################################################################################################################
def length(something):
    if hasattr(something, "__len__"):
        try:
            return len(something)
        except:
            return len(np.atleast_1d(something))
    else:
        return 1