#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
calibration
==========
Simple ccd calibration routine.

Copyright (c) 2018-2019 Moises Castillo

All rights reserved.
"""

import ccdproc
import numpy as np
from astropy.stats import mad_std
from astropy import units as u
import os
import warnings

def add_units(path, units='adu'):
    """
    Check for units in fits header. If none then add Keyword BUNIT with default = adu

    Args:
        path: path to fits file
        units: optional, default set to adu useful options
                adu, photon, and electron
    Returns:
        None
    """
    ic1 = ccdproc.ImageFileCollection(location=path)
    # check for units in header
    if 'bunit' not in ic1.summary.colnames:
        newpath = path+'/unitsadded'
        os.mkdir(newpath)
        warnings.warn('Overwriting Original Header!!!')
        for hdu in ic1.hdus(overwrite=True):
            hdu.header['bunit'] = 'adu'

def validate_units(imagefilecollection):
    if 'bunit' not in imagefilecollection.summary.colnames:
        raise Exception('Missing BUNIT keyword in header. Consider using calibration.add_units()')

# Median Combine
def median_combine(path):
    """
    Median combine with from a list within a directory.

    Args:
        path: path to directory of fits files to be combined

    Returns:
        a median combined CCDData object
    """

    ic1 = ccdproc.ImageFileCollection(location=path)

    validate_units(ic1)

    framelist = ic1.files_filtered(BUNIT='adu', include_path=True)

    median_frame = ccdproc.combine(framelist,
                                    method='median',
                                    sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                    sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                                    mem_limit=350e6
                                    )
    return median_frame

# Make Master Flat
def make_master_flat(path, flat_dark):
    """
    Make master flat.

    Args:
        path: path to directory of fits files to be combined

        flat_dark: CCDData object of dark for flat

    Returns:
        master flat as a median combined CCDData object
    """

    combined_flat = median_combine(path)
    exposure_time = combined_flat.header['EXPTIME']
    dark_exposure = flat_dark.header['EXPTIME']
    master_flat = ccdproc.subtract_dark(
                            combined_flat,
                            flat_dark,
                            dark_exposure=dark_exposure*u.s,
                            data_exposure=exposure_time*u.s,
                            exposure_unit=u.s,
                            scale=True
                            )
    return master_flat

def make_masters(path):
    darkpath = path+'/dark'
    flatpath = path+'/flat'
    master_dark = median_combine(darkpath)
    master_flat = make_master_flat(flatpath, master_dark)
    return master_dark, master_flat

def reduce(ccd, master_dark, master_flat, readnoise, gain):
    """
    reduce raw data frame for science.
        Steps: make dark, make flat, then complete reduction using
        ccdproc.ccd_process which dark subtracts and flat correct

    Args:
        path: path to calibration frames directory to be combined

        ccd: CCDData to be reduced

        readnoise: Quantity, read noise from CCD

        gain: Quantity, gain from CCD

    Returns:
        reduced frame
    """
    # get exposure times for dark and data
    dark_exposure = master_dark.header['EXPTIME']*u.s
    data_exposure = ccd.header['EXPTIME']*u.s

    # assign units to gain and readnoise
    gain = gain*u.electron/u.adu
    readnoise = readnoise*u.electron

    # Gain correct dark and flat before processing
    master_dark = ccdproc.gain_correct(master_dark, gain)
    master_flat = ccdproc.gain_correct(master_flat, gain)

    reduced_ccd=ccdproc.ccd_process(
                                    ccd,
                                    oscan=None,
                                    trim=None,
                                    error=True,
                                    master_bias=None,
                                    dark_frame=master_dark,
                                    master_flat=master_flat,
                                    bad_pixel_mask=None,
                                    gain=gain,
                                    readnoise=readnoise,
                                    oscan_median=True,
                                    oscan_model=None,
                                    min_value=None,
                                    dark_exposure=dark_exposure,
                                    data_exposure=data_exposure,
                                    exposure_key=None,
                                    exposure_unit=u.s,
                                    dark_scale=True
                                    )
    return reduced_ccd

if __name__ == '__main__':
    print('calibration')
