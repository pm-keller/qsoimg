#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

selfcalplots.py

Created on: 2023/06/07
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Make self-calibration plots

"""

import os
import argparse
import casaplotms
from pyvirtualdisplay import Display

display = Display(visible=0,size=(2048,2048))
display.start( )

parser = argparse.ArgumentParser()
parser.add_argument('--imdir', default='/rds/user/pmk46/hpc-work/19A-056/imaging', type=str, help='Directory where intermediate imaging data is stored.')
parser.add_argument('--rounds', default=4, type=int, help='Number of rounds of selfcal.')
args = parser.parse_args()

with open('finished.txt', 'r') as file:
    for line in file:
        obsname = line.strip()

        selfcaldir = os.path.join(args.imdir, obsname, "selfcal")

        for rnd in range(args.rounds):
            for mode in ["G", "T", "S", "ST"]:
                selfcalfile = os.path.join(selfcaldir, f"selfcal-{rnd}-{mode}.tb")
                plotfile = os.path.join(selfcaldir, f"selfcal-{rnd}-{mode}.png")

                if os.path.exists(selfcalfile) and not os.path.exists(plotfile):
                    print(plotfile)

                    if "S" in mode:
                        cax = "corr"
                    else:
                        cax = "spw"
                    
                    for antenna in ["1~9", "10~19", "19~26"]:
                        casaplotms.plotms(selfcalfile, xaxis="time", yaxis="phase", xconnector="line", coloraxis=cax, iteraxis="antenna", antenna=antenna, gridrows=3, gridcols=3, plotfile=plotfile, overwrite=True, highres=True, showgui=False)

                        if rnd > 1:
                            plotfileamp = os.path.join(selfcaldir, f"selfcal-{rnd}-{mode}-amp.png")
                            casaplotms.plotms(selfcalfile, xaxis="time", yaxis="amplitude", xconnector="line", coloraxis=cax, iteraxis="antenna", antenna=antenna, gridrows=3, gridcols=3, plotfile=plotfileamp, overwrite=True, highres=True, showgui=False)
