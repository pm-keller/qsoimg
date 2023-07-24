#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

preppeel.py

Created on: 2023/06/13
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Prepare data and source lists for source peeling.

"""

import os
import re
import csv
import argparse
from astropy.coordinates import SkyCoord
import numpy as np
import casatools
import casatasks

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, help='Path to Measurement Set.')
parser.add_argument('--imfits', type=str, help='Path to latest CLEAN image.')
parser.add_argument('--srcs', type=str, help='List of CLEAN components.')
parser.add_argument('--first', type=str, help='FIRST survey regions.')
parser.add_argument('--selfcal', type=str, help='Directory where self-calibration tables are stored.')
parser.add_argument('--rnds', type=int, help='Number of selfcal-rounds to apply.')
args = parser.parse_args()


coords = []
fluxds = []
distances = []

def get_target():
    """ 
    Read source list and return component with the highest flux
    """
    
    coords = []
    fluxds = []

    with open(args.srcs, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row
        
        for row in reader:
            ra = row[2].replace(":", "h", 1).replace(":", "m", 1) + "s"
            dec = row[3].replace(".", "d", 1).replace(".", "m", 1) + "s"

            coordstr = ra + ' ' + dec
            coords.append(coordstr)
            fluxds.append(float(row[4]))

    coords = np.array(coords)

    idx = np.argmax(fluxds)
    return coords[idx], fluxds[idx]

def get_first_region(coord):
    """ 
    Read FIRST regions list and return corresponding region coordinates and size.
    """

    # read CRTF regions files
    with open(args.first) as file:
        regions = [line.rstrip() for line in file]

    distances = []
    coords = []
    sizes = []

    for region in regions[1:]:
        # get RA/DEC
        coordstr = re.findall("[+-]?[0-9]+.[0-9]+deg", region)
        ra = coordstr[0]
        dec = coordstr[1]

        # get region size
        sizes.append(float(re.findall("[0-9]+.[0-9]+arcsec", region)[0][:-6]))

        firstcoord = ra + ' ' + dec
        skycoord0 = SkyCoord(coord, frame='icrs')
        skycoord1 = SkyCoord(firstcoord, frame='icrs')
        coords.append(skycoord1.to_string('hmsdms'))
        distances.append(skycoord0.separation(skycoord1).arcmin)

    idx = np.argmin(distances)

    return coords[idx], sizes[idx]


def make_source_lists(coord, size):
    """ 
    Create two source lists: 1) sources in target region 2) sources outside target region
    """
    
    with open(args.srcs, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row

        rows_in = "Format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='1692896269.520322', MajorAxis, MinorAxis, Orientation\n"
        rows_out = "Format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='1692896269.520322', MajorAxis, MinorAxis, Orientation\n"
        
        for row in reader:
            ra = row[2].replace(":", "h", 1).replace(":", "m", 1) + "s"
            dec = row[3].replace(".", "d", 1).replace(".", "m", 1) + "s"

            skycoord0 = SkyCoord(ra + ' ' + dec, frame='icrs')
            skycoord1 = SkyCoord(coord, frame='icrs')
            distance = skycoord0.separation(skycoord1).arcsec

            if distance < size:
                rows_in += ','.join(row) + '\n'
            else:
                rows_out += ','.join(row) + '\n'

        with open('target.txt', 'w') as f:
            f.write('%s' % rows_in)
        with open('srcsout.txt', 'w') as f:
            f.write('%s' % rows_out)


def get_solint(fluxp, snr=10):
    """
    Return self-calibration solution interval for a given peak flux, total observation time and SNR.
    """

    # get observation time    
    tb = casatools.table()
    tb.open(args.ms)
    time_array = tb.getcol("TIME")
    tb.close()
    ttot = max(time_array) - min(time_array)

    # get image statistics
    stat = casatasks.imstat(args.imfits)

    # per-antenna RMS. 5~sqrt(N_antennas - 3)
    rms_ant = stat['medabsdevmed'] * 5.0

    # self-calibration solution interval
    tself = ttot / (fluxp / rms_ant / snr)**2 

    # double up tself, because it will be used to calibrate on polarisations separately
    tself = max(2, 2 * tself[0])

    with open('tself.txt', 'w') as f:
        f.write('%d' % tself)
    
    return tself

def apply_selfcal_tables():
    """
    Apply self-calibration tables to measurement set.
    """
    gaintables = []
    spwmaps = []

    for rnd in range(args.rnds):
        for mode in ["G", "T", "S", "ST"]:
            gaintable = f"{args.selfcal}/selfcal-{rnd}-{mode}.tb"

            if os.path.exists(gaintable):
                gaintables.append(gaintable)
                if "S" in mode:
                    spwmaps.append([0,0,0,0,0,0,0,0,0])
                else:
                    spwmaps.append([0,1,2,3,4,5,6,7,8])

    casatasks.applycal(args.ms, gaintable=gaintables, spwmap=spwmaps)


if __name__ == "__main__":
    coord, fluxd = get_target()
    coord, size = get_first_region(coord)
    
    with open('region.txt', 'w') as f:
        f.write('%s\n%d' % (coord, size))

    make_source_lists(coord, size)
    get_solint(fluxd)
    apply_selfcal_tables()
    casatasks.setjy(args.ms, fluxdensity=[0, 0, 0, 0], usescratch=True)