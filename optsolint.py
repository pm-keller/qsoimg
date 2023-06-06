#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

optsolint.py

Created on: 2023/05/23
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Find optimal self-calibration solution interval for a measurement set.

"""

import argparse
import casatools
import casatasks
import numpy as np
import pickle

from scipy.interpolate import interp1d

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, help='Measurement set.')
parser.add_argument('--solint_min', type=float, help='Minimum solution interval.')
parser.add_argument('--refant', type=str, help='Reference antenna.', default="ea19")
parser.add_argument('--calmode', type=str, help='Calibration mode "p" or "ap" (phase or phase + amplitude).', default="p")
parser.add_argument('--solmode', type=str, help='Solve mode, e.g. "R", "L1" or "L1R" for more robustness', default="R")
parser.add_argument('--nspw', type=int, help='Number of spectral windows in data.', default=9)
args, unknown = parser.parse_known_args()

def get_obs_time(ms):
    """Get the duration of an observation, i.e. max(t)-min(t) in seconds.

    Parameters
    ----------
    ms : string
        path to measurement set

    Returns
    -------
    float
        observation duration (+ buffer = integration time) (s)
    """
    
    tb = casatools.table()
    tb.open(ms)
    time_array = tb.getcol("TIME")
    exposure = np.nanmean(tb.getcol("EXPOSURE"))
    tb.close()

    return max(time_array) - min(time_array) + exposure


def get_solint_list(tmin, tmax, dt):
    """Get list of possible solution intervals between tmin and tmax seconds.
    A solution interval must be a proper fraction of tmax.

    Parameters
    ----------
    tmin : float
        minimum solution interval (s)
    tmax : float
        maximum solution interval (s)
    dt : float
        integration time (s)

    Returns
    -------
    list of strings
        the solution interval strings that can be used for self-calibration
    """

    solint_list = []
    t = tmax
    c = 0

    while t > tmin:
        c += 1
        t = dt * ((tmax / c) // dt)
        t = tmax / (tmax // t)
        solint_list.append(t)

    solint_list = [f"{solint}s" for solint in np.unique(solint_list)]

    return solint_list
        

def gain_cal_inf(ms, refant="ea19", calmode="p", solmode="R"):
    """Make gaintable with gaintype "G" and solint "inf". 
    The resulting gaintable is used to calibrate an offset between polarisation products.

    Parameters
    ----------
    ms : string
        path to measurement set
    refant : str, optional
        reference antenna, by default "ea19"
    calmode : str, optional
        calibration mode "p" (phase) or "ap" (amplitude & phase), by default "p"
    solmode : str, optional
        solve mode, e.g. "R", "L1" or "L1R" for more robustness, by default "R"

    Returns
    -------
    string
        path to calibration table
    """

    casatasks.gaincal(
        ms,
        caltable="selfcal-Ginf.tb",
        solint="inf",
        refant=refant,
        calmode=calmode,
        solnorm=True,
        solmode=solmode, 
        gaintype="G", 
        combine="spw",
    )

    return "selfcal-Ginf.tb"


def gain_cal_short(ms, gaintable=[], spwmap=[], solint="int", refant="ea19", calmode="p", solmode="R"):
    """Make gaintable with gaintype "G" and a short solint. 
    The resulting gaintable used as a proxy of the true gain as a function of time.

    Parameters
    ----------
    ms : string
        path to measurement set
    gaintable : list
        list of caltables to apply before solving for gains, by default []
    spwmap : list,
        list of spwmaps for the corresponding caltables
    refant : str, optional
        reference antenna, by default "ea19"
    calmode : str, optional
        calibration mode "p" (phase) or "ap" (amplitude & phase), by default "p"
    solmode : str, optional
        solve mode, e.g. "R", "L1" or "L1R" for more robustness, by default "R"

    Returns
    -------
    string
        path to calibration table
    """

    casatasks.gaincal(
        ms,
        caltable="selfcal-Gshort.tb",
        solint=solint,
        refant=refant, 
        calmode=calmode, 
        solnorm=True,
        solmode=solmode, 
        gaintype="G", 
        combine="spw",
        gaintable=gaintable,
        spwmap=spwmap,
    )

    return "selfcal-Gshort.tb"


def get_gain_data(caltable):
    """Get data and time array from a calibration table.

    Parameters
    ----------
    caltable : string
        path to calibration tables

    Returns
    -------
    tuple of numpy masked arrays
        data and time arrays
    """

    tb = casatools.table()
    tb.open(caltable)
    
    ants = tb.getcol("ANTENNA1")
    data, flags, time = [], [], []

    for ant in np.unique(ants):
        t1 = tb.query(f"ANTENNA1=={ant}")
        data.append(t1.getcol("CPARAM"))
        flags.append(t1.getcol("FLAG"))
        time.append(t1.getcol("TIME"))     
        t1.close() 

    data = np.ma.array(data, mask=flags)[:, :, 0].T 
    time = np.mean(time, axis=0)

    tb.close()

    return (data, time)

def get_window_sizes(N):
    """Get window sizes for data smoothing.


    Parameters
    ----------
    Nmax : int
        Number of data points

    Returns
    -------
    list of int
        list of window sizes
    """

    dt_list = []

    for nt in range(1, N+1):
        dt_list.append(N // nt)

    return np.unique(dt_list)[:-1]


def moving_average(x, y, dt):
    """Compute moving average

    Parameters
    ----------
    x : array
        x-data
    y : array
        y-data
    dt : int
        window size

    Returns
    -------
    tuple of arrays
        smoothed x and y data
    """

    xavg, yavg = [], []
    N = len(x)
    nt = N // dt

    for i in range(nt):
        xavg.append(np.nanmean(x[i*dt:(i+1)*dt], axis=0))
        yavg.append(np.nanmean(y[i*dt:(i+1)*dt], axis=0))

    return xavg, yavg


def estimate_variances_from_mavg(time, data):
    """Estimate the variance of different solution intervals using moving averages.

    Parameters
    ----------
    time : array
        time array
    data : array
        data array

    Returns
    -------
    tuple of arrays
        the solution intervals and the variances due to interpolation error, noise, and both effects combined.
    """
    varI, varN = [], []
    var0 = np.nanvar((data[:, 0] - data[:, 1]) / np.sqrt(2))

    N = len(time)
    dt_list = get_window_sizes(N)

    for dt in dt_list:
        data_nt, time_nt = [], []
        nt = N // dt

        time_nt, data_nt = moving_average(time, data, dt)
        data_nt = interp1d(time_nt, data_nt, fill_value="extrapolate", axis=0)(time)[:, [1, 0]]

        varI.append(np.nanvar(data - data_nt) - var0 * (1 + nt / N / dt))
        varN.append(var0 * nt / N / dt)

    varI.append(np.nanvar(data) - var0 * (1 + 1 / N))
    varN.append(var0 / N)

    varTot = np.array(varI) + np.array(varN)
    solint0 = np.mean(time[1:] - time[:-1])
    solints = solint0 * np.array(np.append(dt_list, [N,]))

    return (solints, np.array(varI), np.array(varN), varTot)


def get_optimal_solints(solints, varI, varN, ttot, nspw=9):
    """Get optimal solution intervals for different calibration modes
    and estimate the possible dynamic range after self-calibration
    with these solution intervals.

    Parameters
    ----------
    solints : list
        list of solution intervals
    varI : list
        list of variances due to interpolation errors
    varN : list
        list of variances due to noise
    ttot : float
        total observation time in seconds
    nspw : int, optional
        number of spectral windows in data, by default 9

    Returns
    -------
    dictionary
        solution intervals for calibration modes "G", "T", "S" and "ST".
    """

    # make solution intervals a proper fraction of the observation time
    solints = ttot / (ttot // solints)

    # calmode "G"
    varTotG = varI + varN * 2 * nspw
    iminG = np.argmin(varTotG)
    solintG = solints[iminG]

    # calmode "T"
    varTotT = varI + varN * nspw
    iminT = np.argmin(varTotT)
    solintT = solints[iminT]

    # calmode "G", combine="SPW"
    varTotS = varI + varN
    iminS = np.argmin(varTotS)
    solintS = solints[iminS]

    # calmode "T", combine="SPW"
    varTotST = varI + varN / 2
    iminST = np.argmin(varTotST)
    solintST = solints[iminST]

    return {"G": solintG, "T": solintT, "S": solintS, "ST": solintST}


if __name__ == "__main__":
    ttot = get_obs_time(args.ms)
    solint0 = ttot / (ttot // args.solint_min)
    ginf = gain_cal_inf(args.ms, refant=args.refant, calmode=args.calmode, solmode=args.solmode)
    gshort = gain_cal_short(args.ms, ginf, args.nspw*[0,], solint0, refant=args.refant, calmode=args.calmode, solmode=args.solmode)
    data, time = get_gain_data(gshort)
    res = estimate_variances_from_mavg(time, data)
    solint_dict = get_optimal_solints(*res[:3], ttot, args.nspw)

    with open('solint.pickle', 'wb') as file:
        pickle.dump(solint_dict, file)

    print(solint_dict)
