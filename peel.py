#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

peel.py

Created on: 2023/06/12
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Peel a source from the visibility data.

"""

import os
import shutil
import numpy as np 
import casatools
import casatasks
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, default='', help='Path to measurement set from which to peel the source. The model column should contain a calibrated model of the source to be peeled.')
parser.add_argument('--selfcal', default='', type=str, help='Directory containing self-calibration tables for source to be peeled.')
args = parser.parse_args()


def invert_gains(gaintable):
    """
    Compute inverse of gaintable.
    """

    if os.path.exists(gaintable + ".rev"):
        shutil.rmtree(gaintable + ".rev")

    shutil.copytree(gaintable, gaintable + ".rev")

    tb = casatools.table()
    tb.open(gaintable)
    gains = tb.getcol("CPARAM")
    tb.close()

    tb.open(gaintable + ".rev", nomodify=False)
    tb.putcol("CPARAM", 1.0 / gains)
    tb.close()

    return gaintable + ".rev"

def apply_selfcal_tables(reverse=False):
    """
    Apply self-calibration tables to measurement set.
    """
    gaintables = []
    spwmaps = []

    for mode in ["G", "T", "S", "ST"]:
        gaintable = os.path.join(args.selfcal, f"selfcal-0-{mode}.tb")
        
        if os.path.exists(gaintable):
            if reverse:
                gaintable = invert_gains(gaintable)
        
            gaintables.append(gaintable)
            if "S" in mode:
                spwmaps.append([0,0,0,0,0,0,0,0,0])
            else:
                spwmaps.append([0,1,2,3,4,5,6,7,8])

    casatasks.applycal(args.ms, gaintable=gaintables, spwmap=spwmaps)

if __name__ == "__main__":
    print("APPLY CALIBRATION")
    apply_selfcal_tables()

    print("SUBTRACT SOURCE")
    casatasks.uvsub(args.ms)

    #print("APPLY INVERSE CALIBRATION")
    #apply_selfcal_tables(reverse=True)