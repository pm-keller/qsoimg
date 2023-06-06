#!/usr/bin/env python
# -*-coding:utf-8 -*-
"""

saveflags.py

Created on: 2023/05/31
Author: Pascal M. Keller
Affiliation: Cavendish Astrophysics, University of Cambridge, UK
Email: pmk46@cam.ac.uk

Description: Post-imaging flagging on residual

"""

import os
import casatasks
import casatools
import numpy as np
import argparse

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, help='Path to measurement set.')
parser.add_argument('--name', type=str, help='Name of flag version.')
parser.add_argument('--when', type=str, help='Before/After flagging.')
args = parser.parse_args()

def get_data_arrays(ms, datacol="DATA", nchan=64):
    """
    Get data arrays from measurement set.
    """

    # get array from table
    tb = casatools.table()
    tb.open(ms)
    time = tb.getcol("TIME")
    spw = tb.getcol("DATA_DESC_ID")
    uvw = tb.getcol("UVW")
    data = tb.getcol(datacol)
    flags = tb.getcol("FLAG")
    data = np.ma.array(data, mask=flags)
    tb.close()

    ntimes = np.unique(time).shape[0]
    nspw = np.unique(spw).shape[0]

    # reshape
    data = data.reshape(2, 64, nspw, ntimes, -1)
    data = np.ma.array([data[:, i % nchan, i // nchan] for i in range(64 * nspw)])
    data = np.moveaxis(data, 0, -1)

    # compute uv-distance averaged over times and frequencies
    uvdist = np.ma.mean(np.sqrt(np.sum(np.abs(uvw)**2, axis=0)).reshape(nspw, ntimes, -1), axis=(0, 1))

    # reorder data
    data = data[:, :, np.argsort(uvdist)]

    return time, spw, uvw, data


def get_versionnames(ms):
    """Get a list of flagversions

    Parameters
    ----------
    ms : str
        path to measurement set

    Returns
    -------
    list of str
        flag version names
    """

    if os.path.exists(ms + ".flagversions"):
        return [vname.split(".")[-1] for vname in os.listdir(ms + ".flagversions")]
    else:
        return []
    

def save(ms, name, when, field=""):
    """save flagversion and flagging summary

    Parameters
    ----------
    ms : str
        path to measurement set
    name : str
        name of processing step status
    when : str
        'before' or 'after' processing step
    field : str
        field to summarise
    """
    versionnames = get_versionnames(ms)

    print(f"\nsave flag version: {when}_{name}_flags")
    if f"{when}_{name}_flags" in versionnames:
        casatasks.flagmanager(
            ms,
            mode="delete",
            versionname=f"{when}_{name}_flags",
        )
    casatasks.flagmanager(ms, mode="save", versionname=f"{when}_{name}_flags")

    print(f"\nwriting summary to: {name}_flags_summary_{when}.npy")
    summary_1 = casatasks.flagdata(
        ms,
        mode="summary",
        field=field,
        name=f"{when}_{name}_flagging",
    )
    np.save(f"{name}_flags_summary_{when}.npy", summary_1)


def diagnostic_plots_2D(data, when="before"):
    """
    Plot Amplitude vs Baseline and Frequency, and vs Time and Frequency.
    """

    data_avg_1 = np.abs(np.ma.mean(data, axis=(0, 1)))
    data_avg_2 = np.abs(np.ma.mean(data, axis=(0, 2)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    ax1.imshow(data_avg_1, aspect="auto")
    ax1.set_title("Averaged Visibilities")
    ax1.set_xlabel("Channel")
    ax1.set_ylabel("Baseline Number")

    im = ax2.imshow(data_avg_2, aspect="auto")
    ax2.set_title("Averaged Visibilities")
    ax2.set_xlabel("Channel")
    ax2.set_ylabel("Time Integration")

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.02, 0.8])
    plt.colorbar(im, cax=cax, label="Amplitude (Jy)")
    
    plt.savefig(f"{when}_flagging_2D.png")
    plt.close()


def diagnostic_plots_1D(data, when="before"):
    """ 
    Plot Amplitude vs Frequency.
    """
    
    fig = plt.figure()

    data_avg_1 = np.abs(np.ma.mean(data, axis=(0, 1)))

    plt.plot(np.arange(data_avg_1.shape[-1]), data_avg_1.T)
    plt.xlim([0, data_avg_1.shape[-1]-1])
    plt.xlabel("Channel")
    plt.ylabel("Amplitude (Jy)")
    plt.title("Polarisation and Time-Averaged Visibilities")

    plt.savefig(f"{when}_flagging_1D.png")
    plt.close()

if __name__=="__main__":
    # plot diagnostic plots
    _, _, _, data = get_data_arrays(args.ms)
    diagnostic_plots_1D(data, f"{args.when}_{args.name}")
    diagnostic_plots_2D(data, f"{args.when}_{args.name}")

    save(args.ms, name=args.name, when=args.when)