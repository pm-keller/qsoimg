#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

mktb.py

Created on: 2023/06/07
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Make data table

"""

import os
import csv
import pickle
import numpy as np
from astropy.io import fits
from astropy.modeling import models, fitting

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--finaldir', default='/rds/user/pmk46/hpc-work/19A-056/final', type=str, help='Directory where final images are stored.')
parser.add_argument('--imdir', default='/rds/user/pmk46/hpc-work/19A-056/imaging', type=str, help='Directory where intermediate imaging data is stored.')
parser.add_argument('--cell', default=0.3, type=str, help='Cell size of image in arcsec')
parser.add_argument('--cutout', default=2048, type=str, help='Pixel size of cutout images.')

args = parser.parse_args()


def get_cutout(data, cutout):
    """ 
    Get a cutout of the centre of a 2D data array.
    """

    N = data.shape[0]
    pc = N // 2

    nstart = max(pc - cutout // 2, 0)
    nend = min(pc + cutout // 2, N)
    return data[nstart:nend, nstart:nend]


with open('finished.txt', 'r') as file:
    for line in file:
        obsname = line.strip()

        # FITS image files
        imfits = os.path.join(args.finaldir, "images", f"{obsname}-res.im-MFS-image.fits")
        imfitsV = os.path.join(args.finaldir, "stokes-V", f"{obsname}-res-Stokes-V.im-MFS-image.fits")

        # self-calibration solution interval files
        selfcaldir   = os.path.join(args.imdir, f"{obsname}", "selfcal")
        solint_file = os.path.join(selfcaldir, "solint.pickle")

        if os.path.exists(imfitsV) and os.path.exists(solint_file):
            image_data = fits.open(imfits)[0].data[0, 0]
            image_data_V = fits.open(imfitsV)[0].data[0, 0]
            header = fits.open(imfits)[0].header
            beam = header["BMIN"] * header["BMAJ"] / (header["CDELT2"] * 2.355)**2

            # get location and value of peak flux in image
            peak_flux_px = np.unravel_index(np.argmax(image_data), image_data.shape)
            peak_flux = (image_data[peak_flux_px] * 1000).round(1) # mJy
            centre_px = np.array(image_data.shape) // 2
            peak_flux_dist = (args.cell * np.sqrt(np.sum((peak_flux_px - centre_px)**2)) / 60).round(2) # arcmin

            # generate cutout images
            Nsmall = 16
            image_data_small = get_cutout(image_data, Nsmall)
            image_data = get_cutout(image_data, args.cutout)
            image_data_V = get_cutout(image_data_V, args.cutout)

            # cutout statistics
            rms = (np.sqrt(np.mean(image_data**2)) * 1e6).round(1) # uJy
            rmsV = (np.sqrt(np.mean(image_data_V**2)) * 1e6).round(1) # uJy
            mad = np.median(1.4826 * np.abs(image_data - np.median(image_data)) * 1e6).round(1) # uJy
            madV = np.median(1.4826 * np.abs(image_data_V - np.median(image_data_V)) * 1e6).round(1) # uJy
            rms_ratio_IV = (rms / rmsV).round(2)
            mad_ratio_IV = (mad / madV).round(2)
            fluxpeak = np.max(image_data_small * 1000).round(3) # mJy
            snr_peak = (fluxpeak / mad * 1000).round(2)

            # Fit 2D Gaussian to central part of image
            x = np.arange(-Nsmall//2, Nsmall//2)
            y = np.arange(-Nsmall//2, Nsmall//2)
            X, Y = np.meshgrid(x, y)

            # Create a 2D Gaussian model
            gaussian_init = models.Gaussian2D(fluxpeak / 1000, 0, 0, 4, 4, bounds={"x_mean": (-4, 4), "y_mean": (-4, 4), "x_stddev": (2, 8),  "y_stddev": (2, 8), "amplitude": (0, 0.1)})

            # Define the fitter and fit the model to the data
            fitter = fitting.LevMarLSQFitter()
            gaussian_fit = fitter(gaussian_init, X, Y, image_data_small)

            # Extract the fitted parameters
            amplitude = gaussian_fit.amplitude.value * 1000
            x_mean = gaussian_fit.x_mean.value
            y_mean = gaussian_fit.y_mean.value
            x_stddev = gaussian_fit.x_stddev.value
            y_stddev = gaussian_fit.y_stddev.value
            fluxint = (amplitude * x_stddev * y_stddev / beam).round(2)
            snr_amp = (amplitude / mad * 1000).round(2)

            print(fluxpeak, amplitude, mad, snr_peak, snr_amp)

            # load self-calibration solution interval
            with open(solint_file, "rb") as file:
                solint_dict = pickle.load(file)

            # load flagging statistics
            flaggingdir = os.path.join(args.imdir, f"{obsname}", "flagging")
            flagging_before = os.path.join(flaggingdir, "flagging_flags_summary_before.npy")
            flagging_after = os.path.join(flaggingdir, "rflagavg_flags_summary_after.npy")

            flagstats_before = np.load(flagging_before, allow_pickle=True).item()
            flagstats_after = np.load(flagging_after, allow_pickle=True).item()

            flagstat_list = []

            for spw in flagstats_before["spw"]:
                flagstat_list.append((flagstats_after["spw"][spw]["flagged"] - flagstats_before["spw"][spw]["flagged"]) / flagstats_before["spw"][spw]["total"] * 100)

            flagstat_avg = np.mean(flagstat_list).round(2)

            # make table
            table_row = (obsname, rms, mad, rmsV, madV, (rms/rmsV).round(2), (mad/madV).round(2), x_mean.round(2), y_mean.round(2), x_stddev.round(2), y_stddev.round(2), amplitude.round(2), fluxint, fluxpeak, snr_peak, snr_amp, peak_flux, peak_flux_dist, solint_dict["ST"].round(1), flagstat_avg)


            # append table rows
            with open('table.csv', 'a', newline='') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write the new row to the CSV file
                writer.writerow(table_row)