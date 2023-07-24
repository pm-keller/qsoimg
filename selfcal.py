#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""
  
selfcal.py
  
Created on: 2023/06/14
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk
  
Description: Self-Calibrtaion
  
"""

import argparse
import casatasks
import casatools
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, help='Path to Measurement Set.')
parser.add_argument('--solintfile', type=str, help='Path solution interval file.')
parser.add_argument('--calmode', type=str, default='p', help='Calibration mode. One of "p" (phase) or "ap" (amplitude and phase)')
parser.add_argument('--snr', type=str, default=10, help='Min. SNR for self-calibration.')
parser.add_argument('--rnd', type=str, default=0, help='Self-calibration round.')
args = parser.parse_args()

# make solution intervals that are a fraction of the total observation time
with open(args.solintfile, "rb") as file:
  solint_dict = pickle.load(file)

# get observation time
tb = casatools.table()
tb.open(args.ms)
time_array = tb.getcol("TIME")
tb.close()
ttot = max(time_array) - min(time_array)

print(solint_dict, ttot)

# G: gaincal without combining data, T: combine polarisations, S: combine spectral windows
# The solution interval needs to be adjusted accordingly and is doubled for calmode = "ap"
N = 1
if args.calmode == "ap":
  N = 2

solintG = N * max(solint_dict["G"], 36 * solint_dict["ST"])
solintT = N * max(solint_dict["T"], 18 * solint_dict["ST"])
solintS = N * max(solint_dict["S"], 2 * solint_dict["ST"])
solintST = min(N * solint_dict["ST"], ttot)

gaintables, spwmaps = [], []

if (ttot // solintG):
    solintG = (ttot + 2) / (ttot // solintG)

    # SPW self-calibrationm combine polarisations
    casatasks.gaincal(
      args.ms,
      caltable=f"selfcal-{args.rnd}-G.tb",
      solint=f"{solintG}s",
      refant="ea19",
      calmode=args.calmode,
      solnorm=True,
      solmode="R",
      gaintype="G",
      minsnr=args.snr
    )
    gaintables.append(f"selfcal-{args.rnd}-G.tb")
    spwmaps.append([0,1,2,3,4,5,6,7,8])

if (ttot // solintT):
    solintT = (ttot + 2) / (ttot // solintT)
    
    # SPW self-calibrationm combine polarisations
    casatasks.gaincal(
      args.ms,
      caltable=f"selfcal-{args.rnd}-T.tb",
      solint=f"{solintT}s",
      refant="ea19",
      calmode=args.calmode,
      solnorm=True,
      solmode="R",
      gaintype="T",
      minsnr=args.snr,
      gaintable=gaintables,
      spwmap=spwmaps
    )
    gaintables.append(f"selfcal-{args.rnd}-T.tb")
    spwmaps.append([0,1,2,3,4,5,6,7,8])

if (ttot // solintS):
    solintS = (ttot + 2) / (ttot // solintS)

    # SPW self-calibrationm combine polarisations
    casatasks.gaincal(
      args.ms,
      caltable=f"selfcal-{args.rnd}-S.tb",
      solint=f"{solintS}s",
      refant="ea19",
      calmode=args.calmode,
      solnorm=True,
      solmode="R",
      gaintype="G",
      combine="spw",
      minsnr=args.snr,
      gaintable=gaintables,
      spwmap=spwmaps
    )
    gaintables.append(f"selfcal-{args.rnd}-S.tb")
    spwmaps.append([0,0,0,0,0,0,0,0,0])

if (ttot // solintST):
    # SPW self-calibrationm combine polarisations
    casatasks.gaincal(
      args.ms,
      caltable=f"selfcal-{args.rnd}-ST.tb",
      solint=f"{solintST}s",
      refant="ea19",
      calmode=args.calmode,
      solnorm=True,
      solmode="R",
      gaintype="T",
      combine="spw",
      minsnr=args.snr,
      gaintable=gaintables,
      spwmap=spwmaps
    )
    gaintables.append(f"selfcal-{args.rnd}-ST.tb")
    spwmaps.append([0,0,0,0,0,0,0,0,0])