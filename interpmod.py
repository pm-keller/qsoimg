#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

fitmod.py

Created on: 2023/05/01
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Fit a polynomial model to model FITS images, resample at new frequencies and write to disk.

"""

import re
import glob
import copy 
import time
import argparse
import numpy as np
from scipy.interpolate import interp1d
from multiprocessing import Pool, cpu_count

from astropy.io import fits
import casatools

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, help='Path to measurement set.')
args = parser.parse_args()

# Get a list of files that need to be interpolated and read data and frequency arrays
name = re.search("QSO-J[0-9]+[+-][0-9]+", args.ms)[0]
filelist = glob.glob(f"{name}-*0*-model.fits")

# relax glob expression if filelist is empty
if len(filelist) == 0:
    filelist = glob.glob(f"*-0*-model.fits")

print(f"reading files: {filelist}")

hdulist = [fits.open(filename)[0] for filename in filelist]
data_array = np.array([hdu.data for hdu in hdulist])
freq_array = np.array([hdu.header['CRVAL3'] for hdu in hdulist])

# sort by frequency
sortidx = np.argsort(freq_array)
freq_array = freq_array[sortidx]
data_array = data_array[sortidx]

# get frequencies from measurement set
print(f"reading frequencies from {args.ms}")
md = casatools.msmetadata()
md.open(args.ms)
new_freq_range = np.concatenate([md.chanfreqs(spw=i) for i in range(md.nspw())])
md.close()

# create a list of headers for the new data
headerlist = [copy.deepcopy(hdulist[0].header) for i in range(len(new_freq_range))]

# add a new frequency to each header
for i in range(len(headerlist)):
    headerlist[i]["CRVAL3"] = new_freq_range[i]

# prepare the data for the polynomial fit
shape = data_array.shape
data_array = data_array.reshape(shape[0], -1)
plist = []

# fit a polynomial along frequency at every image pixel and store coefficients in list
print(f"\ninterpolate data")
t1 = time.time()

datamax = np.max(data_array, axis=0)
idx_list = np.where(datamax > 0)[0]

print(f"interpolating {len(idx_list)} pixels\n")

# function for parallelising over pixels
def interp(idx):
    data = data_array.T[idx]
    return interp1d(freq_array, data, kind='cubic', fill_value="extrapolate")

with Pool(cpu_count()) as p:
    interplist = p.map(interp, idx_list)

t2 = time.time()

# function for parallelising over frequencies
def write_im(i):
    freq = new_freq_range[i]
    new_data = np.zeros((data_array.shape[-1],))

    for j, idx in enumerate(idx_list):
        new_data[idx] = interplist[j](freq)

    # put the data back into its proper shape
    new_data = new_data.reshape(shape[1:])

    # write to file
    hdu = fits.PrimaryHDU(data=new_data, header=headerlist[i])
    hdu.writeto('out-{:04d}-model.fits'.format(i), overwrite=True)

# write the new fits images to the disk
print("writing data to FITS")

t3 = time.time()

with Pool(cpu_count()) as p:
    p.map(write_im, range(len(new_freq_range)))

t4 = time.time()

print("\ninterpolation: {:.2f}s, evaluation and writing: {:.2f}s".format(t2-t1, t4-t3))


"""
### An alternative which might be faster, but with higher memory requirements ###

new_data = []

for data in data_array.T:
    p = np.polyfit(freq_array, data, 4)
    pdata = np.polyval(p, new_freq_range)
    new_data.append(pdata)

new_data = np.swapaxes(new_data, 0, 1).reshape((len(new_freq_range),) + shape[1:])

def write_im(i):
    hdu = fits.PrimaryHDU(data=new_data[i].data, header=headerlist[i])
    hdu.writeto('out-{:04d}-model.fits'.format(i), overwrite=True)
    
with Pool(cpu_count()) as p:
    p.map(write_im, range(len(new_freq_range)))
"""
