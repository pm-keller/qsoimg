#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

makeim.py

Created on: 2023/04/16
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Make a CASA image cube from a list of FITS files

"""

import casatools
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fits', type=str, nargs="+", help='List of FITS images which will be concatenated and converted to CASA images.')
parser.add_argument('-o', '--out', type=str, help='Name of output image.')
args = parser.parse_args()

N = len(args.fits)
imlist = []

print("convert FITS to CASA images:")

for i, im in enumerate(args.fits):
   print(im)
   ia = casatools.image()
   ia.fromfits(f"part-{i}.im", im, overwrite=True)
   imlist.append(f"part-{i}.im")

# concatenate image
print(f"concatenate images: {args.out}")
ia = casatools.image()
ia.imageconcat(args.out, imlist, axis=-1, relax=True, reorder=True, overwrite=True)

