#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

tclean.py

Created on: 2023/05/08
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Execute CASA's tclean on the data.

"""

# execute tclean
tclean(ms, 
  imagename=imname,
  imsize=1024, 
  cell="0.3arcsec", 
  pblimit=-0.1,
  niter=10000, 
  threshold=threshold,
  nsigma=10,
  gain=0.1, 
  weighting="briggs",
  robust=-0.5, 
  uvtaper="120klambda", 
  uvrange=["3~200 klambda"],
  gridder="wproject", 
  wprojplanes=2,
  specmode="mfs",
  deconvolver="mtmfs",
  nterms=2,
  outlierfile=outlierfile,
  mask=mask,
  parallel=True,
)