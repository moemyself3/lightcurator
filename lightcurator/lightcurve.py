#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
lightcurve
==========
Script to plot lightcurves.

Copyright (c) 2018-2019 Moises Castillo

All rights reserved.
"""

__version__ = '1.0.0a0'

import astroalign as aa
from os import listdir
from os.path import isfile, join
from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from astropy.table import Table
import multiprocessing as mp
from astropy.time import Time
import matplotlib.pyplot as plt
import subprocess

def plot(object_table):
    """
    Create a timeseries plot of stars.

    Args:
        candidate_catalog: a catalog of candidate variables with timestamp

    Returns:
        a plot of the candiate variables
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

    plt.plot(t.plot_date,flux,'o')
    return object_table

def align(object_table):
    """
    Align images from image_list.

    Args:
        object_table: a list of file names

    Returns:
        a table of objects that have matching field of view
    """

    image_list = object_table['path']
    # use file list to make data objects
    data = []
    date = []
    for item in image_list:
        sci_data = fits.open(item)[0]
        data.append(np.asarray(sci_data.data))
        date.append(sci_data.header['DATE-OBS'])
    object_table['date'] = date
    object_table['data'] = data

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
    mypath = '/Users/mcastillo/Desktop/20170420_test'
    image_list = Table()
    image_list['path'] = [join(mypath,f) for f in listdir(mypath) if isfile(join(mypath, f))]
    return image_list

def paralign(object_table):
    """
    Parallel processing: Align images from image_list.

    Args:
        object_table: a list of file names

    Returns:
        a table of objects that have matching field of view
    """

    image_list = object_table['path']
    # use file list to make data objects
    data = []
    date = []
    for item in image_list:
        sci_data = fits.open(item)[0]
        data.append(np.asarray(sci_data.data))
        date.append(sci_data.header['DATE-OBS'])
    object_table['date'] = date
    object_table['data'] = data

    # pick reference image
    ref_index = len(data)//2
    ref_img = data[ref_index]
    object_table['ref'] = [ref_img] * len(data)

    pool = mp.Pool()
    # align images to reference image
    aligned = []
    args_align = zip(object_table['ref'], object_table['data'])
    aligned = pool.starmap(aa.align_image, args_align)
    object_table['aligned'] = aligned

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

    pool = mp.Pool()

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

def do_deepsky(object_table):
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

def makepic(data):
    vmin,vmax = mm.minmax(data)
    plt.imshow(data, cmap='Greys', origin='lower', vmin=vmin, vmax=vmax)

def astrometry(data):
    #Take aligned image and add wcs
    subprocess.run(['solve-field',data])

    #get RA,DEC coordinates to do a skymatch

if __name__ == '__main__':
    print('main does nothing')
