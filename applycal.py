#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

applycal.py

Created on: 2023/05/31
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Apply self-calibration tables to original measurement set. 
Self-calibration tables must be in format selfcal-<round>-<mode>.tb,
where mode is one of "G", "T", "S" or "ST".

"""

import os
import shutil
import casatasks
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, help='Path to measurement set.')
parser.add_argument('--nself', type=int, default=3, help='Number of self-calibration tables to apply. Default is 3.')
parser.add_argument('--dir', type=str, help='Directory where self-calibration tables are stored.')
parser.add_argument('--nspw', type=int, default=9, help='Number of spectral windows. Default is 9.')
args = parser.parse_args()

gaintables = []
spwmaps = []

# create list of gaintables and spwmaps
for round in range(args.nself):
    for mode in ["G", "T", "S", "ST"]:
        gaintbpath = os.path.join(args.dir, f"selfcal-{round}-{mode}.tb")

        if os.path.exists(gaintbpath):
            gaintables.append(gaintbpath)

            if mode in ["G", "T"]:
                spwmaps.append(list(range(0, args.nspw)))
            elif mode in ["S", "ST"]:
                spwmaps.append(args.nspw * [0,])

# copy measurement set to self-calibration directory
mstmp = os.path.join(args.dir, f"{args.ms.split('/')[-1][:-3]}-selfcal-{args.nself-1}-tmp.ms")
ms = os.path.join(args.dir, f"{args.ms.split('/')[-1][:-3]}-selfcal-{args.nself-1}.ms")

if os.path.exists(mstmp):
    shutil.rmtree(mstmp)

shutil.copytree(args.ms, mstmp)

# apply gaintables
print(gaintables)
print(spwmaps)
casatasks.applycal(mstmp, gaintable=gaintables, spwmap=spwmaps)

# freeze calibrated measurement set
if os.path.exists(ms):
    shutil.rmtree(ms)
if os.path.exists(f"{ms}.flagversions"):
    shutil.rmtree(f"{ms}.flagversions")

casatasks.split(mstmp, ms, datacolumn="corrected")
shutil.rmtree(mstmp)

