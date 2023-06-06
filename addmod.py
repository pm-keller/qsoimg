#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

addmod.py

Created on: 2023/04/16
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Add model from one measurement set to that of another

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
tb.unlock()
tb.close()

# replace data with model and add to existing model column
tb.open(args.add, nomodify=False)
casatasks.uvsub(args.add) # This creates a corrected data column if it does not exist
tb.putcol("CORRECTED_DATA", model)
casatasks.uvsub(args.add, reverse=True)

# transfer updated model from data column to model column
model = tb.getcol("CORRECTED_DATA")
tb.putcol("MODEL_DATA", model)

# unlock and close table
tb.unlock()
tb.close()
