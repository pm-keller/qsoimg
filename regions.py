#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

regions.py

Created on: 2023/04/27
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Read a CRTF regions file, select sources that need to be treated as outliers and make CASA mask.

"""

import os
import re
import shutil
import casatasks
import casatools
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, help='Path to measurement set.')
parser.add_argument('--image', type=str, help='FITS image for creating CASA mask.')
parser.add_argument('--regions', type=str, help='Path to regions file. Must be in CRTF format.')
parser.add_argument('--first', type=str, default="", help='Path to FIRST regions file. Must be in CRTF format.')
parser.add_argument('--cell', type=float, default=0.3, help='Cell size in arcsec, e.g. 0.3.')
parser.add_argument('--imsize', type=int, default=12144, help='Image size in pixels, e.g. 8192.')
parser.add_argument('--radius', type=float, help='Radius in arcmin above which sources will be treated as outlier sources.')
args = parser.parse_args()


# get phase centre of measurement set
md = casatools.msmetadata()
md.open(args.ms)
pc = md.phasecenter()
c0 = SkyCoord(ra=pc["m0"]["value"], dec=pc["m1"]["value"], unit=pc["m0"]["unit"])

# read CRTF regions files
with open(args.regions) as file:
    regions = [line.rstrip() for line in file]

clist = []
ralist, declist, malist, sizelist = [], [], [], []
sizes = np.array([32, 64, 128, 256, 512, 1024])


for region in regions[1:]:
    # get RA/DEC
    ra = re.search("[0-9]+:[0-9]+:[0-9]+.[0-9]+", region)[0]
    dec = re.search("[+-][0-9]+.[0-9]+.[0-9]+.[0-9]+", region)[0]
    orientation = re.search("[0-9]+.[0-9]+deg", region)[0]

    masize = re.findall("[0-9]+.[0-9]+arcsec", region)
    axlen = (float(masize[0].replace("arcsec", "")) / 0.3, float(masize[1].replace("arcsec", "")) / 0.3)
    imsize = 2 * max(axlen)
    dsize = sizes - imsize
    pidx = np.where(dsize > 0)
    imsize = np.min(sizes[pidx])

    ma = f"ellipse [[{int(imsize) // 2}pix, {int(imsize) // 2}pix], [{axlen[0]}pix, {axlen[1]}pix], {orientation}]"

    # add RA/DEC to lists
    ralist.append(ra)
    declist.append(dec)
    malist.append(ma)
    sizelist.append(imsize)

    # replace the first two '.' with ':' in dec string
    dec = list(dec)
    dec[4]=":"
    dec[7]=":"

    # add sky coordinates to list
    clist.append(SkyCoord(ra+" "+"".join(dec), unit=(u.hourangle, u.deg)))

# separate outlier and field sources
out_regions = [regions[0],]
ma_regions = [regions[0],]

# outlier file template string
outstr1 = """imagename=src{:d}
imsize=[{:d},{:d}]
cell=[0.3arcsec,0.3arcsec]
phasecenter=J2000 {:s} {:s}"""

outstr2 = """{:s} {:s} {:d}"""

outstr1_list = []
outstr2_list = []
cnt = 0

# remove existing coordinate file
if os.path.exists(args.regions + "-out.txt"):
    os.remove(args.regions + "-out.txt")

for i, ci in enumerate(clist):
    dra, ddec = c0.spherical_offsets_to(ci)
    dc = max((np.abs(dra.to(u.arcmin)), np.abs(ddec.to(u.arcmin))))

    if dc > (args.radius - 1) * u.arcmin:
        out_regions.append(regions[i+1])
        outstr1_list.append(outstr1.format(cnt, sizelist[i], sizelist[i], ralist[i], declist[i]))
        outstr2_list.append(outstr2.format(ralist[i], declist[i], sizelist[i]))
        cnt += 1

    else:
        ma_regions.append(regions[i+1])

# use FIRST regions to create mask, if they are provided
if args.first != "" and os.path.exists(args.first):
    with open(args.first) as file:
        ma_regions = [regions[0],] + [line.rstrip() for line in file][1:]

print(ma_regions)

# save new CRTF region files for the corresponding lists
with open(args.regions + "-out.txt", mode='wt') as wfile:
    wfile.write('\n'.join(out_regions))
with open(args.regions + "-mask.txt", mode='wt') as wfile:
    wfile.write('\n'.join(ma_regions))
with open(args.regions + "-outlierfile.txt", mode='wt') as wfile:
    wfile.write('\n'.join(outstr1_list))
with open(args.regions + "-outliers.txt", mode='wt') as wfile:
    wfile.write('\n'.join(outstr2_list))

# make casa mask
if os.path.exists(args.regions+"-tmp.im"):
    shutil.rmtree(args.regions+"-tmp.im")
if os.path.exists(args.regions+"-mask.im"):
    shutil.rmtree(args.regions+"-mask.im")
    
casatasks.importfits(args.image, imagename=args.regions+"-tmp.im")
casatasks.makemask(inpimage=args.regions+"-tmp.im", inpmask=args.regions+"-mask.txt",  mode="copy", output=args.regions+"-mask.im")
shutil.rmtree(args.regions+"-tmp.im")