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

from scipy.constants import c
from scipy.special import erf, erfinv
from scipy.optimize import curve_fit

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astroquery.ipac.irsa import Irsa
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--finaldir', default='/rds/user/pmk46/hpc-work/19A-056/final', type=str, help='Directory where final images are stored.')
parser.add_argument('--imdir', default='/rds/user/pmk46/hpc-work/19A-056/imaging', type=str, help='Directory where intermediate imaging data is stored.')
parser.add_argument('--cell', default=0.3, type=str, help='Cell size of image in arcsec')

args = parser.parse_args()

# list of detected sources
det_names = ["QSO-J0100+2802", "QSO-J1034-1425", "QSO-J1427+3312", "QSO-J1429+5447", "QSO-J1558-0724", "QSO-J1602+4228", "QSO-J2318-3113", "QSO-J2348-3054"]

def get_target_pixel(hdu, target_coord_str):
    """ 
    Convert sky coordinates to pixel coordinates.
    """

    # Extract the image data and header
    image_header = hdu.header

    # Create a WCS object from the FITS header
    wcs = WCS(image_header)

    # Convert the target coordinate to SkyCoord object
    target_coord = SkyCoord(target_coord_str, unit=(u.hourangle, u.deg))

    # Convert the target coordinate to pixel coordinates
    target_pixel = wcs.all_world2pix(target_coord.ra.deg, target_coord.dec.deg, 0, 0, 1)
    
    # Round the pixel coordinates to the nearest integer
    target_pixel = np.array(target_pixel).astype(float) - np.array([1, 1, 1, 1])

    return target_pixel

def get_cutout(data, cutout, pc):
    """ 
    Get a cutout of the centre of a 2D data array.
    """

    N = data.shape[0]

    nstart0 = max(pc[0] - cutout // 2, 0)
    nend0 = min(pc[0] + cutout // 2 + 1, N)
    nstart1 = max(pc[1] - cutout // 2, 0)
    nend1 = min(pc[1] + cutout // 2 + 1, N)

    return data[nstart0:nend0, nstart1:nend1]

def gaussian_2d(xy, amplitude, x0, y0):
    """ 
    2D Gaussian function
    """

    x, y = xy
    x0 = float(x0)
    y0 = float(y0)
    a = (np.cos(theta)**2) / (2 * sigmin**2) + (np.sin(theta)**2) / (2 * sigmaj**2)
    b = -(np.sin(2 * theta)) / (4 * sigmin**2) + (np.sin(2 * theta)) / (4 * sigmaj**2)
    c = (np.sin(theta)**2) / (2 * sigmin**2) + (np.cos(theta)**2) / (2 * sigmaj**2)
    g = amplitude * np.exp(- (a * ((x - x0)**2) + 2 * b * (x - x0) * (y - y0) + c * ((y - y0)**2)))

    return g.ravel()

def gauss_fit_2d(XY, data, bounds):
    """ 
    Fit a 2D Gaussian to image cutout data
    """
    # Perform the Gaussian fit
    popt, pcov = curve_fit(gaussian_2d, XY, data.ravel(), p0=(0, 0, 0), bounds=bounds)

    # Extract the fitted parameters
    return popt    


# read Spitzer data
spitzer_name = []
spitzer_S = []
spitzer_errors = []

with open('tables/spitzer.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=" ")
    for row in reader:
        spitzer_name.append(row[0])
        spitzer_S.append(float(row[1]))
        spitzer_errors.append(float(row[2]))

spitzer_name = np.array(spitzer_name)
spitzer_S = np.array(spitzer_S)
spitzer_errors = np.array(spitzer_errors)

# read quasar catalog daa
qso_catalog = QTable.read('tables/Quasar_catalog_Banados+16_Matsuoka+19a_Matsuoka+19b_Wang+18_Wang+19_Reed+19_Yang+20.txt', format='ascii')

# define cosmology 
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)   

with open('qsolist.txt', 'r') as file:
    for line in file:
        obsname = line.strip()

        catalog_name = list(obsname)
        catalog_name[3] = "_"
        catalog_name = "".join(catalog_name)

        # read QSO catalogue
        row = qso_catalog[qso_catalog["QSO-Name"] ==  catalog_name]
        z = row["Redshift"].value[0]
        Mo = row["M_1450"].value[0]
        ra = row["RA"].value[0]
        dec = row["Dec"].value[0]
        dl = cosmo.luminosity_distance(z).to("m").value
        Lsol = 3.828e26
        Lo = 10**((4.83 - Mo) / 2.5) 

        # get Spitzer flux density
        spitzer_idx = np.where(spitzer_name == catalog_name)[0]

        if len(spitzer_idx) > 0:
            spitzer = spitzer_S[spitzer_idx[0]]
            spitzer_err = spitzer_errors[spitzer_idx[0]]
        else: 
            spitzer = np.nan
            spitzer_err = np.nan

        spitzer_AB = (8.9 - 2.5 * np.log10(spitzer / 1e6))
        spitzer_AB_err = 2.5 * spitzer_err / spitzer / np.log(10)

        # get ALLWISE flux density
        irsa_query = Irsa.query_region(f"{ra} {dec}", radius=2*u.arcsec, catalog="allwise_p3as_psd")
        allwise = irsa_query["w1mpro"].value
        allwise_err = irsa_query["w1sigmpro"].value

        # query reject table
        if len(allwise) == 0:
            irsa_query = Irsa.query_region(f"{ra} {dec}", radius=2*u.arcsec, catalog="allwise_p3as_psr")
            allwise = irsa_query["w1mpro"].value
            allwise_err = irsa_query["w1sigmpro"].value

            if len(allwise) == 1:
                print(obsname, " REJECT TABLE")

        if len(allwise) == 1:
            allwise = allwise[0]
            allwise_err = allwise_err[0]
        else:
            allwise = np.nan
            allwise_err = np.nan
        if not isinstance(allwise_err, (int, float, np.number)):
            allwise = np.nan
            allwise_err = np.nan       

        # convert to AB magnitudes
        allwise_AB = allwise + 2.699

        # FITS image files
        imfits = os.path.join(args.finaldir, "images", f"{obsname}-res.im-MFS-image.fits")

        # self-calibration solution interval files
        selfcaldir   = os.path.join(args.imdir, f"{obsname}", "selfcal")
        solint_file = os.path.join(selfcaldir, "solint.pickle")

        if os.path.exists(imfits):
            # read fits
            hdu = fits.open(imfits)[0]
            header = hdu.header
            image_data = hdu.data[0, 0]

            # get synthesised beam parameters
            sigmin = header["BMIN"] / (2.355 * header["CDELT2"])
            sigmaj = header["BMAJ"] / (2.355 * header["CDELT2"])
            beam = sigmin * sigmaj
            theta = np.pi - header["BPA"] * np.pi / 180
            
            # get location and value of peak flux in image
            peak_flux_px = np.unravel_index(np.argmax(image_data), image_data.shape)
            peak_flux = (image_data[peak_flux_px] * 1000).round(1) # mJy
            centre_px = np.array(image_data.shape) // 2
            peak_flux_dist = (args.cell * np.sqrt(np.sum((peak_flux_px - centre_px)**2)) / 60).round(2) # arcmin

            # get target pixel
            pc = get_target_pixel(hdu, f"{ra} {dec}")[[0, 1]]
            pcr = np.round(pc, 0).astype(int)
            pcd = pc - pcr
            fluxpeak = (image_data[pcr[0], pcr[1]] * 1e6).round(3) # uJy
            blank = (image_data[pcr[0], pcr[1] - 200] * 1e6).round(3) # uJy

            # generate cutout images
            Nsmall, Ntiny = 21, 13
            image_data_small = get_cutout(image_data, Nsmall, pcr)
            image_data_tiny = get_cutout(image_data, Ntiny, pcr)
            image_data_blank = get_cutout(image_data, Nsmall, [pcr[0], pcr[1] - 200])
            image_data = get_cutout(image_data, 41, pcr)

            # cutout statistics
            rms = (np.sqrt(np.mean(image_data**2)) * 1e6).round(1) # uJy
            mad = np.median(1.4826 * np.abs(image_data - np.median(image_data)) * 1e6).round(1) # uJy

            # frequency 
            #freq = header["CRVAL3"]
            freq = 1.6e9
            f0 = (1 + z) * freq

            # spectral index
            alpha = -0.75


            # Fit 2D Gaussian to central part of the image
            x = np.arange(-Nsmall//2, Nsmall//2) + 1.0 + pcd[0]
            y = np.arange(-Nsmall//2, Nsmall//2) + 1.0 + pcd[1]
            X, Y = np.meshgrid(x, y)
            
            amplitude, x_mean, y_mean = gauss_fit_2d((X, Y), image_data_small, bounds=([-1e5, -3, -3], [1e5, 3, 3]))            
            amplitude_fixed, _, _ = gauss_fit_2d((X, Y), image_data_small, bounds=([-1e5, -1e-16, -1e-16], [1e5, 1e-16, 1e-16]))
            amplitude_blanks, _, _ = gauss_fit_2d((X, Y), image_data_blank, bounds=([-1e5, -1e-16, -1e-16], [1e5, 1e-16, 1e-16]))

            # fitte in-band spectrum
            if obsname == "QSO-J1427+3312":
                amplitude = 1281 * 1e-6
                fluxpeak = 1281
                alpha = -1.43
            elif obsname == "QSO-J1429+5447":
                amplitude = 3003 * 1e-6
                fluxpeak = 3003
                alpha = -0.73
            elif obsname == "QSO-J2318-3113":
                amplitude = 637 * 1e-6
                fluxpeak = 637
                alpha = -1.48

            # get integrated flux densities
            fluxint = (amplitude * sigmin * sigmaj / beam)
            fluxsum = np.sum(image_data_tiny) / beam / (2 * np.pi)

            # estimate S/N
            snr_amp = (amplitude / mad * 1e6).round(2)

            # integrated radio luminosity
            log_L_int = np.log10(4 * np.pi * dl**2 * 1e-26 * (5e9 / f0)**alpha * fluxint / (1+z) / 1e6 / Lsol * 5e9)
            log_L_int_err = mad / fluxint / np.log(10)


            # compute radio luminosity
            log_L_nu = np.log10(4 * np.pi * dl**2 * 1e-26 * (5e9 / f0)**alpha * np.abs(fluxpeak) / (1+z) / 1e6)
            log_L_peak = log_L_nu + np.log10(5e9 / Lsol)
            log_L_peak_err = mad / fluxpeak / np.log(10)
            snr_peak = (fluxpeak / mad).round(2)
            Lsgn = np.sign(fluxpeak)

            # compute upper limit
            E = erf(fluxpeak / mad / np.sqrt(2))
            #ul = mad * np.sqrt(2) * (erfinv(1 - (1 + E) * 0.002699796) + fluxpeak / mad / np.sqrt(2))
            ul = 3 * mad
            log_L_ul = np.log10(4 * np.pi * dl**2 * 1e-26 * (5e9 / f0)**alpha * ul  / 1e6 / (1+z) / Lsol * 5e9)

            # compute the luminosity at 4400A
            if not np.isnan(spitzer):
                log_Lo = np.log10(4 * np.pi * dl**2 * 1e-26 * (36000 / 4400 / (1+z))**-0.5 * spitzer / 1e6 / (1+z) / Lsol * c / 4.4e-7)
                log_Lo_err = spitzer_err / spitzer / np.log(10)
            elif not np.isnan(allwise):
                fw1 = 10**((8.926 - allwise_AB) / 2.5)
                log_Lo = np.log10(4 * np.pi * dl**2 * 1e-26 * (33677 / 4400 / (1+z))**-0.5 * fw1 / (1+z) / Lsol * c / 4.4e-7)
                log_Lo_err = allwise_err / 2.5
            else:
                log_Lo = (34.1 - Mo) / 2.5 + np.log10((1450 / 4400)**-0.5 / Lsol * c / 4.4e-7)
                log_Lo_err = 0.3 / np.log(10)
            
            
            # compute radio loudness parameter
            R = 10 ** (log_L_peak - log_Lo) * (c / 4.4e-7 / 5e9)
            Rul = 10 ** (log_L_ul - log_Lo) * (c / 4.4e-7 / 5e9)
            Rerr = np.log(10) * R * np.sqrt(log_L_peak_err**2 + log_Lo_err**2)
                
            # load self-calibration solution interval
            if os.path.exists(solint_file):
                with open(solint_file, "rb") as file:
                    solint = pickle.load(file)["ST"].round(1)
            else:
                solint = "-"

            # load flagging statistics
            flaggingdir = os.path.join(args.imdir, f"{obsname}", "flagging")
            flagging_before = os.path.join(flaggingdir, "flagging_flags_summary_before.npy")
            flagging_after = os.path.join(flaggingdir, "rflagavg_flags_summary_after.npy")

            try:
                flagstats_before = np.load(flagging_before, allow_pickle=True).item()
                flagstats_after = np.load(flagging_after, allow_pickle=True).item()
                flagstat_list = []

                for spw in flagstats_before["spw"]:
                    flagstat_list.append((flagstats_after["spw"][spw]["flagged"] - flagstats_before["spw"][spw]["flagged"]) / flagstats_before["spw"][spw]["total"] * 100)

                flagstat_avg = np.mean(flagstat_list).round(2)
            except:
                flagstat_avg = 0

            # make tables
            table1_row = (obsname, rms, mad, x_mean.round(2), y_mean.round(2), np.round(sigmin, 2), np.round(sigmaj, 2), np.round(amplitude * 1e6, 0), amplitude_fixed * 1e6, np.round(amplitude_blanks * 1e6, 0), np.round(fluxint * 1e6, 0), np.round(fluxsum * 1e6, 0), fluxpeak, blank, snr_peak, snr_amp, peak_flux, peak_flux_dist, solint, flagstat_avg)
            table2_row = (obsname, z, freq, amplitude * 1e6, mad, ul, Mo, allwise_AB, allwise_err, spitzer_AB, spitzer_AB_err, Lsgn, log_L_nu, log_L_peak, log_L_peak_err, log_L_ul, log_Lo, log_Lo_err, R, Rerr, Rul)

            # append table rows
            with open('tables/table1.csv', 'a', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(table1_row)  

            with open('tables/table2.csv', 'a', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(table2_row)  