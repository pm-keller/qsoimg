#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

sensplots.py

Created on: 2023/06/07
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Plot image noise RMS as a function of self-calibration round.

"""

import os
import numpy as np
import argparse
import matplotlib.pyplot as plt

from astropy.io import fits


parser = argparse.ArgumentParser()
parser.add_argument('--imdir', default='/rds/user/pmk46/hpc-work/19A-056/imaging', type=str, help='Directory where intermediate imaging data is stored.')
parser.add_argument('--finaldir', default='/rds/user/pmk46/hpc-work/19A-056/final/images', type=str, help='Directory where final images are stored.')
parser.add_argument('--rounds', default=3, type=int, help='Number of rounds of selfcal.')
parser.add_argument('--cutout', default=2048, type=int, help='Image cutout size.')
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

        cleandir = os.path.join(args.imdir, obsname, "clean")
        plotfile = os.path.join(cleandir, f"{obsname}-sensplot.png")
        mad_list = []

        if os.path.exists(os.path.join(cleandir, f"{obsname}-selfcal-0.im-MFS-image.fits")) and os.path.exists(os.path.join(args.finaldir, f"{obsname}-res.im-MFS-image.fits")) and not os.path.exists(plotfile):
            imfits = os.path.join(cleandir, f"{obsname}.im-MFS-image.fits")
            image_data = fits.open(imfits)[0].data[0, 0]
            image_data = get_cutout(image_data, args.cutout)
            mad_list.append(np.median(1.4826 * np.abs(image_data - np.median(image_data)) * 1e6).round(1)) # uJy
            
            for rnd in range(args.rounds):
                imfits = os.path.join(cleandir, f"{obsname}-selfcal-{rnd}.im-MFS-image.fits")
                image_data = fits.open(imfits)[0].data[0, 0]
                image_data = get_cutout(image_data, args.cutout)
                mad_list.append(np.median(1.4826 * np.abs(image_data - np.median(image_data)) * 1e6).round(1)) # uJy

            imfits = os.path.join(args.finaldir, f"{obsname}-res.im-MFS-image.fits")
            image_data = fits.open(imfits)[0].data[0, 0]
            image_data = get_cutout(image_data, args.cutout)
            mad_list.append(np.median(1.4826 * np.abs(image_data - np.median(image_data)) * 1e6).round(1)) # uJy

            rnds = np.arange(args.rounds+2).astype(int)

            print(plotfile)

            plt.plot(rnds, mad_list)
            plt.xlabel("Self-Calibration Round", color="k")
            plt.ylabel(r"Noise RMS ($\mu$Jy)")
            plt.xlim([0, args.rounds+1])
            plt.savefig(plotfile)
            plt.close()