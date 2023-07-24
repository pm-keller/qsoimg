#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

putmod.py

Created on: 2023/04/16
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Put data from one measurement set to the model of another

"""

import casatasks
import casatools
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--take', type=str, help='measurement set to take model from')
parser.add_argument('--add', type=str, help='measurement set to add model to')
args = parser.parse_args()

tb = casatools.table()

# get model data to be added
tb.open(args.take)
model = tb.getcol("MODEL_DATA")
tb.close()

# put data to model column
casatasks.setjy(args.add, fluxdensity=[0, 0, 0, 0], usescratch=True)
tb.open(args.add, nomodify=False)
tb.putcol("MODEL_DATA", model)
tb.close()
