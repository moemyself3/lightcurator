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
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astropy import units as u
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, Column
from astropy.visualization import ZScaleInterval
from astropy.time import Time
from astroquery.vizier import Vizier
from photutils import DAOStarFinder
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
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
        mean, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)

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

def make_catalog(aligned_table):
    # get wcs and RA, DEC from aligned
    for aligned_path, cat_path in zip(aligned_table['aligned_path'], aligned_table['cat_path']):
        with fits.open(aligned_path) as hdu:
            w = WCS(hdu[0].header)
            cat = ascii.read(cat_path)
            pixcrd = list(zip(cat['xcentroid'], cat['ycentroid']))
            world = w.wcs_pix2world(pixcrd, 0)
            ra, dec = zip(*world)
            col_ra = Column(ra, name='ra')
            col_dec = Column(dec, name='dec')
            cat.add_columns([col_ra, col_dec])
            # Update catalog
            ascii.write(cat, cat_path, overwrite=True)

def makelist(mypath, column_title='path'):
    """
    Make list of objects

    Args:
        path: path to directory to of objects
    Returns:
        a table of object paths
    """
    # Import file list
    image_list = Table()
    image_list[column_title] = [join(mypath,f) for f in listdir(mypath) if isfile(join(mypath, f))]
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

    # get date of observation
    date_obs = object_table['date'][0]
    date_obs = datetime.strptime(date_obs,'%Y-%m-%dT%H:%M:%S').strftime('%Y-%m-%dT%H%M%S')

    # make root path
    root = 'light_collection/'+date_obs

    # pick reference image
    ref_index = len(data)//2
    ref_img = data[ref_index]
    object_table['ref'] = [ref_img] * len(data)

    process_limit = mp.cpu_count() - 1
    pool = mp.Pool(processes=process_limit)
    strictlist = [strict]*len(data)
    # align images to reference image
    args_align = zip(object_table['path'], object_table['filtered'], object_table['ref'], object_table['sigclip'], itertools.repeat(root), strictlist)
    output = pool.starmap(tryregister, args_align)
    object_table['aligned'], object_table['sigclip'] = zip(*output)
    pool.close()
    pool.join()
    return object_table

def parextract_table(object_table):
    """
    Create catalog of sources

    Args:
        aligned_objects: CCDData objects that have matching field of view

    Returns:
        a catalog of sources with their timestamp
    """
    sources = []
    aligned_objects = object_table['aligned']
    object_table['kwargs_stat'] = {'sigma':3.0, 'maxiters':5}

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

def parextract_path(root, object_table=None):
    """
    Create catalog of sources

    Args:
        aligned_objects: CCDData objects that have matching field of view

    Returns:
        a catalog of sources with their timestamp

    Todo:
        Need to make a single table with corresponding paths for aligned images
         for time stamps.
    """
    if not object_table:
        object_table = Table()
    alignpath = root+'/aligned/'
    object_table = makelist(alignpath, 'aligned_path')

    aligned_objects = object_table['aligned_path']
    kwargs_stat = itertools.repeat({'sigma':3.0, 'maxiters':5})

    process_limit = mp.cpu_count() - 1
    pool = mp.Pool(processes=process_limit)

    args_stat = aligned_objects
    iterable_stat = zip(args_stat, kwargs_stat)
    sigma_clipped_stat = pool.starmap(sigma_clipped_wrapper, iterable_stat)
    mean, median, std = zip(*sigma_clipped_stat)

    # Next Subtract background and use DAOStarFinder
    args_dao = zip(aligned_objects, median)
    kwargs_dao = []
    for item in std:
        kwargs_dao.append({'fwhm':10, 'threshold': 5.*item})
    iterable_dao = zip(args_dao, kwargs_dao)

    sources = pool.starmap(DAOStarFinder_wrapper,iterable_dao)

    #object_table['sources'] = sources
    pool.close()
    pool.join()

    catpath = []
    name_source_list = zip(aligned_objects, sources)
    for path, source_table in name_source_list:
        filename, _ = splitext(basename(path))
        catfile = root+'/cats/aligned/'+filename+'.cat'
        ascii.write(source_table, catfile, overwrite=True)
        catpath.append(catfile)
    object_table['cat_path'] = catpath
    return object_table

def make_master_cat(path, root):
    """
    Create catalog of sources from deepsky

    Args:
        path: str - path to deepsky

    Returns:
        a catalog of sources with their timestamp
    """
    kwargs_stat = {'sigma':3.0, 'maxiters':5}

    with fits.open(path) as hdu:
        data = hdu[0].data
        mean, median, std = sigma_clipped_stats(data, **kwargs_stat)
        kwargs_dao = {'fwhm':10.0, 'threshold':5*std}
        print('mean: '+str(mean))
        print('median: '+ str(median))
        print('std: '+str(std))
        # Next Subtract background and use DAOStarFinder
        print('kwargs_dao: ' +  str(kwargs_dao))
        daofind = DAOStarFinder(**kwargs_dao)
        sources = daofind(data - median)

    filename, _ = splitext(basename(path))
    catfile = root+'/deepsky/'+filename+'.cat'
    ascii.write(sources, catfile, overwrite=True)

    # add RA and DEC to cat file
    args_makecatalog = Table()
    args_makecatalog['aligned_path'] = [path]
    args_makecatalog['cat_path'] = [catfile]
    make_catalog(args_makecatalog)

    return catfile


def sigma_clipped_wrapper(path, kwargs):
    with fits.open(path) as fitsfile:
        data = fitsfile[0].data
        sigma_clipped_stat = sigma_clipped_stats(data, **kwargs)
    return sigma_clipped_stat

def DAOStarFinder_wrapper(args_dao, kwargs):
    path = args_dao[0]
    median = args_dao[1]
    with fits.open(path) as fitsfile:
        data = fitsfile[0].data
        daofind = DAOStarFinder(**kwargs)
        return daofind(data - median)

def do_lightcurve():
    print('makelist...')
    object_table = makelist('')

    print('makelist complete...')

    #shorttable = object_table #[0:10]
    print('alignment...')
    object_table = paralign(object_table)

    print('alignment complete...')
    print('extraction...')
    object_table = parextract_table(object_table)

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
        sci_data.close()
    object_table['date'] = date
    object_table['data'] = data
    object_table = hotpixfix(object_table)
    return object_table

def hotpixfix_wrapper(sci_data, sigclip=5):
    return ccdproc.cosmicray_lacosmic(sci_data, sigclip=sigclip)[0]

def hotpixfix(object_table, sigclip=5):
    object_table['sigclip'] = sigclip
    process_limit = mp.cpu_count() - 1
    pool = mp.Pool(processes=process_limit)
    args_hotpixfix = zip(object_table['data'], object_table['sigclip'])
    filtered = pool.starmap(hotpixfix_wrapper, args_hotpixfix)
    object_table['filtered'] = filtered
    pool.close()
    pool.join()
    return object_table

def tryregister(path, source, target, sigclip, root, strict=True):
    attempts = 0
    while attempts < 3:
        try:
            aligned = aa.register(source, target)
            hdu = fits.PrimaryHDU(aligned)
            hdul = fits.HDUList([hdu])
            hdr = hdul[0].header

            with fits.open(path) as sci_data:
                hdr['DATE-OBS'] = sci_data[0].header['DATE-OBS']

            filename, _ = splitext(basename(path))
            aligned_fits = root+'/aligned/'+filename+'_aligned.fits'
            hdul.writeto(aligned_fits)
            return aligned, sigclip
        except aa.MaxIterError:
            if strict:
                return np.zeros_like(source), sigclip
            if not strict:
                sigclip+=5
                source = hotpixfix_wrapper(path,sigclip)
                attempts += 1
        except Exception:
                return np.zeros_like(source), sigclip
    return np.zeros_like(source), sigclip

def lightcurator(mypath, parallel=False):
    """
    All in one function to create lightcurves

    Args:
        path: string, path to root directory that contains timeseries data

    Returns:
        a deepsky and aligned data fits files with WCS and catalogs with xy and RA/DEC coordinates
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
    #object_table = object_table[0:100]
    ref_index = len(object_table['date'])//2

    with fits.open(object_table['path'][ref_index]) as sci_data:
         ref_img = hotpixfix_wrapper(sci_data[0].data)

    # setup directories by date and time of beginning of observation
    date_obs = object_table['date'][0]
    date_obs = datetime.strptime(date_obs,'%Y-%m-%dT%H:%M:%S').strftime('%Y-%m-%dT%H%M%S')
    root, _ = setup_dirs(date_obs)

    # read, filter, align images
    sigclip = 5
    print('is parallel: ' + str(parallel))
    if not parallel:
        print('deepsky process: serial')
        start = time.time()
        deepsky=np.zeros_like(ref_img)
        for item in object_table['path']:
            with fits.open(item) as sci_data:
                filtered = hotpixfix_wrapper(sci_data[0].data)
            aligned, _ = tryregister(item, filtered, ref_img, sigclip, root)
            deepsky = deepsky + aligned
        print('deepsky complete!')
        end = time.time()
        print('time: '+str(end - start))
    else:
        print('deepsky process: parallel')
        start = time.time()
        process_limit = mp.cpu_count() - 1
        pool = mp.Pool(processes=process_limit)

        args_do_deepsky = zip(object_table['path'], itertools.repeat(ref_img), itertools.repeat(sigclip), itertools.repeat(root))
        deepsky = sum(pool.starmap(do_align, args_do_deepsky))
        pool.close()
        pool.join()
        print('parallel processed deepsky complete!')
        end = time.time()
        print('time: '+str(end - start))

    makepic(deepsky, root+'/deepsky/deepsky_preview')
    hdu = fits.PrimaryHDU(deepsky)
    hdul = fits.HDUList([hdu])
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    deepskyfits = root+'/deepsky/deepsky_'+timestamp+'.fits'
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
    alignpath = root+'/aligned/'
    for aligned_file in listdir(alignpath):
        aligned_filepath = alignpath+aligned_file
        with fits.open(aligned_filepath) as aligned_hdu:
            aligned_data = aligned_hdu[0].data
            hdr = aligned_hdu[0].header
            header['DATE-OBS'] = hdr['DATE-OBS']
            fits.writeto(aligned_filepath, aligned_data, header, output_verify='fix', overwrite=True)
            #header.tofile(aligned_filepath, endcard=True, overwrite=True)


    # Source extraction
    aligned_table = parextract_path(root)

    # Not necessary to continue
    # Need to sort figure out how to correlate them to the original data
    # object_table = hstack([object_table, aligned_table])

    # Match Cat
    make_catalog(aligned_table)

    master_cat_path = make_master_cat(fitsfile_wcs, root)
    master_cat = ascii.read(master_cat_path)
    master_ra = master_cat['ra']
    master_dec = master_cat['dec']
    master_sky = SkyCoord(ra=master_ra*u.degree, dec=master_dec*u.degree)

    matched_path = []
    for catpath in aligned_table['cat_path']:
        c = ascii.read(catpath)
        c_ra = c['ra']
        c_dec = c['dec']
        c_sky = SkyCoord(ra=c_ra*u.degree, dec=c_dec*u.degree)
        idx, d2d, d3d = c_sky.match_to_catalog_sky(master_sky)
        t = d2d < 1*u.arcsec
        c['sep_flag'] = t
        c['idx'] = idx
        ct = c.dtype
        names, dtypes = zip(*ct.descr)
        matched_cat = Table(names = names, dtype=dtypes)
        for row, sep_flag in enumerate(c['sep_flag']):
            if sep_flag:
                matched_cat.add_row(c[row].as_void())
        filename, _ = splitext(catpath)
        filepath = root+'/cats/matched/'+basename(filename)+'_match.cat'
        ascii.write(matched_cat, filepath, overwrite=True)
        matched_path.append(filepath)
    aligned_table['matched_path'] = matched_path


    names = list(range(len(master_cat)))
    names.append('DATE-OBS')
    dtypes = list(itertools.repeat('float64', len(master_cat)))
    dtypes.append('<U19')

    aligned_table.sort('matched_path')

    row  = [{} for i in range(len(aligned_table['matched_path']))]

    for num, path in enumerate(aligned_table['matched_path']):
        cat = ascii.read(path)
        row[num]=list(itertools.repeat(np.nan, len(master_cat)))
        for ele, idx in enumerate(cat['idx']):
            for source_id in range(len(master_cat)):
                if source_id == idx:
                    row[num][idx]=cat['flux'][ele]
        aligned_image_path, _ = splitext(basename(path))
        aligned_image_path = aligned_image_path.replace('_match','')

        aligned_image_path = root+'/aligned/'+aligned_image_path+'.fits'
        with fits.open(aligned_image_path) as hdu:
            date = hdu[0].header['DATE-OBS']
        row[num].append(date)

    timeseries_catalog = Table(rows = row, names= names, dtype=dtypes)
    index = list(range(0, len(aligned_table['matched_path'])))
    index = Column(index, name='index', dtype='u4')
    timeseries_catalog.add_column(index)
    timeseries_catalog.write(root+'/cats/timeseries_catalog.ecsv', format='ascii.ecsv')

    return object_table, aligned_table, timeseries_catalog

def do_align(path, ref_img, sigclip, root):
    aligned = np.zeros_like(ref_img)

    with fits.open(path) as sci_data:
        filtered = hotpixfix_wrapper(sci_data[0].data)

    aligned, _ = tryregister(path, filtered, ref_img, sigclip, root)
    return aligned

def setup_dirs(date_obs='9999'):
    # make directory for output files
    root = 'light_collection/'+date_obs
    output_directories = ['/deepsky', '/aligned', '/cats/aligned', '/cats/matched']
    output_directories = [root+directory for directory in output_directories]
    for directory in output_directories:
        makedirs(directory, exist_ok=True)
    return root, output_directories

def astrometry(fitsfile):
    # Take aligned image and add wcs
    subprocess.run(['solve-field',fitsfile], check=True)

    fitsfile_wcs = fitsfile.replace('.fits', '.new')

    # check if file exists
    if not isfile(fitsfile_wcs):
        raise FileNotFoundError('The solved fits file does not exist')

    return fitsfile_wcs

def clean(date_obs='9999'):
    _, dirs_to_clean = setup_dirs(date_obs)
    for directory in dirs_to_clean:
        for item in listdir(directory):
            remove(directory+'/'+item)

def query_from_wcs(fits_path, radius=30):
    """
    SIMBAD query of VSX and UCAC catalog from wcs using CRVAL1/CRVAL2 as RA/DEC

    Args:
        path: string, path to wcs file. expects WCS keyword CRVAL1 and CRVAL2
        radius: float, radius of query region in arcmin
    Returns:
        a Table list of VSX, GCVS,  I/340 (UCAC5), I/345 (GAIA DR2) objects in region of query
    """

    with fits.open(fits_path) as hdu:
        header = hdu[0].header
        ra, dec = header['CRVAL1'], header['CRVAL2']

        result = Vizier.query_region(coord.SkyCoord(ra=ra, dec=dec,
                                                unit=(u.deg, u.deg),
                                                frame='fk5'),
                            radius=radius*u.arcmin,
                            catalog=['B/vsx','B/gcvs','I/340','I/345'])
    print(result)
    return result

if __name__ == '__main__':
    print('lightcurve')
