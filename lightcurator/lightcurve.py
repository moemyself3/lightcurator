#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
lightcurve
==========
Script to plot lightcurves.

Copyright (c) 2018-2019 Moises Castillo

All rights reserved.
"""

import astroalign as aa
from os import listdir, makedirs, remove
from os.path import isfile, join, basename, splitext
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from astropy.table import Table
from astropy.visualization import ZScaleInterval
import multiprocessing as mp
from astropy.time import Time
from datetime import datetime
import matplotlib.pyplot as plt
import subprocess
import ccdproc
import itertools
import time

def plot(object_table):
    """
    Create a time series plot of stars.

    Args:
        candidate_catalog: a catalog of candidate variables with timestamp

    Returns:
        a plot of the candidate variables
    """

    total_frames = len(object_table)

    lens = []
    for source in object_table['sources']:
        lens.append(len(source))

    datelist = []
    for i in range(0,total_frames):
        datelist.append([object_table['date'][i]]*lens[i])

    t = []
    for i in range(0,total_frames):
        t=t+datelist[i]
    t = Time(t)

    flux = []
    for i in range(0,total_frames):
        for j in range(0,lens[i]):
            flux.append(object_table['sources'][i]['flux'][j])

    fig = plt.figure()
    plt.plot(t.plot_date,flux,'o')
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    fig.savefig('plot_'+timestamp+'.png')
    plt.close(fig)
    return object_table

def align(object_table):
    """
    Align images from image_list.

    Args:
        object_table: a list of file names

    Returns:
        a table of objects that have matching field of view
    """

    # use file list to make data objects
    object_table = makedata(object_table)
    data = object_table['data']

    # sort table by date
    object_table.sort('date')

    # pick reference image
    ref_index = len(data)//2
    ref_img = data[ref_index]

    # align images to reference image
    aligned = []
    for item in object_table['data']:
        aligned.append(aa.align_image(ref_img,item))
    object_table['aligned'] = aligned

    return object_table

def extract(object_table):
    """
    Create catalog of sources

    Args:
        aligned_objects: CCDData objects that have matching field of view

    Returns:
        a catalog of sources with their timestamp
    """
    sources = []
    aligned_objects = object_table['aligned']

    for data in aligned_objects:
        mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)

        # Next Subtract background and use DAOStarFinder
        daofind = DAOStarFinder(fwhm = 10.0, threshold=5.*std)
        sources.append(daofind(data - median))

    object_table['sources'] = sources

    return object_table

def matchcat(timeseries_catalog):
    """
    Match catalogs using Astropy matching

    Args:
        timeseries_catalog: catalog of sources with their timestamp
    Returns:
        a catalog of matched sources with their timestamp
    """
    candidate_catalog = 1
    return candidate_catalog

def makelist(mypath):
    """
    Make list of objects

    Args:
        path: path to directory to of objects
    Returns:
        a table of object paths
    """
    # Import file list
    image_list = Table()
    image_list['path'] = [join(mypath,f) for f in listdir(mypath) if isfile(join(mypath, f))]
    return image_list

def paralign(object_table, strict=True):
    """
    Parallel processing: Align images from image_list.

    Args:
        object_table: a list of file names

    Returns:
        a table of objects that have matching field of view
    """

    # use file list to make data objects
    object_table = makedata(object_table)
    data = object_table['filtered']

    # sort table by date
    object_table.sort('date')

    # pick reference image
    ref_index = len(data)//2
    ref_img = data[ref_index]
    object_table['ref'] = [ref_img] * len(data)

    process_limit = mp.cpu_count() - 1
    pool = mp.Pool(processes=process_limit)
    strictlist = [strict]*len(data)
    # align images to reference image
    args_align = zip(object_table['path'], object_table['filtered'], object_table['ref'], object_table['sigclip'], strictlist)
    output = pool.starmap(tryregister, args_align)
    object_table['aligned'], object_table['sigclip'] = zip(*output)
    pool.close()
    pool.join()
    return object_table

def parextract(object_table):
    """
    Create catalog of sources

    Args:
        aligned_objects: CCDData objects that have matching field of view

    Returns:
        a catalog of sources with their timestamp
    """
    sources = []
    aligned_objects = object_table['aligned']
    object_table['kwargs_stat'] = {'sigma':3.0, 'iters':5}

    process_limit = mp.cpu_count() - 1
    pool = mp.Pool(processes=process_limit)

    args_stat = aligned_objects
    kwargs_stat = object_table['kwargs_stat']
    iterable_stat = zip(args_stat,kwargs_stat)

    sigma_clipped_stats = pool.starmap(sigma_clipped_wrapper, iterable_stat)
    mean, median, std = zip(*sigma_clipped_stats)

    # Next Subtract background and use DAOStarFinder
    args_dao = zip(aligned_objects, median)
    kwargs_dao = []
    for item in std:
        kwargs_dao.append({'fwhm':10, 'threshold': 5.*item})
    object_table['kwargs_dao'] = kwargs_dao
    iterable_dao = zip(args_dao, kwargs_dao)

    sources = pool.starmap(DAOStarFinder_wrapper,iterable_dao)

    object_table['sources'] = sources
    pool.close()
    pool.join()

    return object_table

def sigma_clipped_wrapper(data, kwargs):
    return sigma_clipped_stats(data,**kwargs)

def DAOStarFinder_wrapper(data, kwargs):
    daofind = DAOStarFinder(**kwargs)
    return daofind(data[1] - data[0])

def do_lightcurve():
    print('makelist...')
    object_table = makelist('')

    print('makelist complete...')

    #shorttable = object_table #[0:10]
    print('alignment...')
    object_table = paralign(object_table)

    print('alignment complete...')
    print('extraction...')
    object_table = parextract(object_table)

    print('extraction complete...')
    print('plotting')
    plot(object_table)
    return object_table

def do_deepsky_table(object_table):
    deepsky = np.zeros_like(object_table['aligned'][0])
    for item in object_table['aligned']:
        deepsky = deepsky + item
    return deepsky

class mm:
    def minmax(self):
        median = []
        mean = []
        median = np.median(self)*.95
        mean = np.mean(self)
        vmin, vmax = (mean - median, mean + median)
        return vmin, vmax

def makepic(data, filename='frame'):
    """
    Create an image from data frame.

    Args:
        data: CCDData object

    Returns:
        frame as png file
    """
    vmin,vmax = mm.minmax(data)
    fig = plt.figure()
    plt.imshow(data, cmap='Greys', origin='lower', vmin=vmin, vmax=vmax)
    timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')
    fig.savefig(filename+'_'+timestamp+'.png')
    plt.close(fig)

def makepic_zscale(data, filename='frame'):
    zscale = ZScaleInterval()
    fig = plt.figure()
    plt.imshow(zscale(data), cmap='Greys', origin='lower')
    fig.savefig(filename+'.png')
    plt.close(fig)

def makedata(object_table):
    data = []
    date = []
    image_list = object_table['path']
    for item in image_list:
        sci_data = fits.open(item)[0]
        data.append(np.asarray(sci_data.data))
        date.append(sci_data.header['DATE-OBS'])

    object_table['date'] = date
    object_table['data'] = data
    object_table = hotpixfix(object_table)
    return object_table

def hotpixfix_wrapper(sci_data, sigclip=4.5):
    return ccdproc.cosmicray_lacosmic(sci_data, sigclip=sigclip)[0]

def hotpixfix(object_table, sigclip=4.5):
    object_table['sigclip'] = sigclip
    process_limit = mp.cpu_count() - 1
    pool = mp.Pool(processes=process_limit)
    args_hotpixfix = zip(object_table['data'], object_table['sigclip'])
    filtered = pool.starmap(hotpixfix_wrapper, args_hotpixfix)
    object_table['filtered'] = filtered
    pool.close()
    pool.join()
    return object_table

def tryregister(path, source, target, sigclip, strict=True):
    attempts = 0
    while attempts < 3:
        try:
            aligned = aa.register(source, target)
            hdu = fits.PrimaryHDU(aligned)
            hdul = fits.HDUList([hdu])
            filename, _ = splitext(basename(path))
            aligned_fits = 'light_collection/aligned/'+filename+'_aligned.fits'
            hdul.writeto(aligned_fits)
            return aligned, sigclip
        except aa.MaxIterError:
            if strict:
                return np.zeros_like(source), sigclip
            if not strict:
                sigclip+=5
                source = hotpixfix_wrapper(path,sigclip)
                attempts += 1
    return np.zeros_like(source), sigclip

def lightcurator(mypath, parallel=False):
    """
    All in one function to create lightcurves

    Args:
        path: string, path to root directory that contains timeseries data

    Returns:
        a deepsky fits file with WCS
    """
    #Make pathlist
    object_table = makelist(mypath)

    date = []
    image_list = object_table['path']

    # Make a time stamp column
    for item in image_list:
        with fits.open(item) as sci_data:
            date.append(sci_data[0].header['DATE-OBS'])
    object_table['date'] = date

    # find reference image
    object_table.sort('date')

    # sample set for testing
    object_table = object_table[0:100]
    ref_index = len(object_table['date'])//2

    with fits.open(object_table['path'][ref_index]) as sci_data:
         ref_img = hotpixfix_wrapper(sci_data[0].data)

    setup_dirs()

    # read, filter, align images
    sigclip = 4.5
    print('is parallel: ' + str(parallel))
    if not parallel:
        print('deepsky process: serial')
        start = time.time()
        deepsky=np.zeros_like(ref_img)
        for item in object_table['path']:
            with fits.open(item) as sci_data:
                filtered = hotpixfix_wrapper(sci_data[0].data)
            aligned, _ = tryregister(item, filtered, ref_img, sigclip)
            deepsky = deepsky + aligned
        print('deepsky complete!')
        end = time.time()
        print('time: '+str(end - start))
    else:
        print('deepsky process: parallel')
        start = time.time()
        process_limit = mp.cpu_count() - 1
        pool = mp.Pool(processes=process_limit)

        args_do_deepsky = zip(object_table['path'], itertools.repeat(ref_img), itertools.repeat(sigclip))
        deepsky = sum(pool.starmap(do_align, args_do_deepsky))
        pool.close()
        pool.join()
        print('parallel processed deepsky complete!')
        end = time.time()
        print('time: '+str(end - start))

    makepic(deepsky, 'light_collection/deepsky/deepsky_preview')
    hdu = fits.PrimaryHDU(deepsky)
    hdul = fits.HDUList([hdu])
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    deepskyfits = 'light_collection/deepsky/deepsky_'+timestamp+'.fits'
    hdul.writeto(deepskyfits)

    # Do astrometry
    print('starting astrometry...')
    fitsfile_wcs = astrometry(deepskyfits)
    print('astrometry complete!')

    # get WCS from .new file
    with fits.open(fitsfile_wcs) as wcsframe:
        thewcs = WCS(wcsframe[0].header)
        header = thewcs.to_header()

    # attach wcs to aligned frames
    alignpath = 'light_collection/aligned/'
    for aligned_file in listdir(alignpath):
        aligned_filepath = alignpath+aligned_file
        with fits.open(aligned_filepath) as aligned_hdu:
            fits.writeto(aligned_filepath, aligned_hdu[0].data, header, overwrite=True)

    # with open.fits(deepskyfits)

    return deepskyfits

def do_align(path, ref_img, sigclip):
    aligned = np.zeros_like(ref_img)

    with fits.open(path) as sci_data:
        filtered = hotpixfix_wrapper(sci_data[0].data)

    aligned, _ = tryregister(path, filtered, ref_img, sigclip)
    return aligned

def setup_dirs():
    # make directory for output files
    root = 'light_collection'
    output_directories = ['/deepsky', '/aligned', '/cats']
    output_directories = [root+directory for directory in output_directories]
    for directory in output_directories:
        makedirs(directory, exist_ok=True)
    return output_directories

def astrometry(fitsfile):
    # Take aligned image and add wcs
    subprocess.run(['solve-field',fitsfile], check=True)

    fitsfile_wcs = fitsfile.replace('.fits', '.new')

    # check if file exists
    if not isfile(fitsfile_wcs):
        raise FileNotFoundError('The solved fits file does not exist')

    # get RA,DEC coordinates to do a skymatch
    return fitsfile_wcs

def clean():
    dirs_to_clean = setup_dirs()
    for directory in dirs_to_clean:
        for item in listdir(directory):
            remove(directory+'/'+item)

if __name__ == '__main__':
    print('main does nothing')
