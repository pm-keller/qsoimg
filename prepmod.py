#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

prepmod.py

Created on: 2023/04/16
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Prepare a measurement set for DP3 model prediction

"""

import shutil
import casatools
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, help='measurement set to take model from')
args = parser.parse_args()

shutil.copytree(args.ms, args.ms + ".tmp")

tb = casatools.table()

# fill datacolumn with zeros
tb.open(args.ms + ".tmp", nomodify=False)
data = tb.getcol("DATA") * 0
tb.putcol("DATA", data)
tb.close()