# -*- coding: utf-8 -*-
'''
JNW's mystery assignment
'''
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
#import matplotlib as mpl
#mpl.use('TkAgg')

import matplotlib.pyplot as plt
from astrobase import astrotess as at
from lightcurve_utils import median_filter_lightcurve
from lightcurve_utils import run_two_periodograms

from astrobase.periodbase import kbls
from astrobase.plotbase import plot_phased_mag_series
from astrobase import periodbase, checkplot
from astrobase.lcmath import phase_magseries, sigclip_magseries
from astrobase.varbase import lcfit
from astrobase.periodbase import get_snr_of_dip

from astrobase.varbase.transits import get_transit_times, \
        given_lc_get_transit_tmids_tstarts_tends



def plot_phased_lcs(times, fluxs, period, epoch,
                    savdir='../results/tic_238176110_winn/',
                    substr='discovery', pb=0.01, phaselim=[-.6,.6]):

    outfile = (
        savdir+'phased_on_{:s}.png'.format(substr)
    )
    plot_phased_mag_series(times, fluxs, period, magsarefluxes=True,
                           errs=None, normto=False, epoch=epoch,
                           outfile=outfile, sigclip=False, phasebin=pb,
                           plotphaselim=phaselim, plotdpi=400)

    outfile = (
        savdir+'phased_on_{:s}_asy_sclip_20-3.png'.format(substr)
    )
    plot_phased_mag_series(times, fluxs, period, magsarefluxes=True,
                           errs=None, normto=False, epoch=epoch,
                           outfile=outfile, sigclip=[20,3], phasebin=pb,
                           plotphaselim=phaselim, plotdpi=400)

    outfile = (
        savdir+'phased_on_{:s}_asy_sclip_3-3.png'.format(substr)
    )
    plot_phased_mag_series(times, fluxs, period, magsarefluxes=True,
                           errs=None, normto=False, epoch=epoch,
                           outfile=outfile, sigclip=[3,3], phasebin=pb,
                           plotphaselim=phaselim, plotdpi=400)




def visually_look_thru(time, flux):

    f, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(12,6))
    ax.scatter(time, flux, c='k', alpha=0.5, label='PDCSAP', zorder=1,
                   s=1.5, rasterized=True, linewidths=0)

    ax.set(xlabel='time [days]', ylabel='relative flux')
    f.tight_layout(h_pad=0, w_pad=0)

def initial_wrangle():
    savdir = '../results/tic_238176110_winn/'
    fname = 'tess2018206045859-s0001-0000000238176110-111-s_llc.fits.gz'
    infile = savdir+fname
    lctype = 'PDCSAP'

    filter_savpath = (savdir+
                      'tess2018206045859-s0001-0000000238176110-111-s'
                      '_llc_PDCSAP_blsfit_phased.png'
    )


    time, flux, err_flux = at.get_time_flux_errs_from_Ames_lightcurve(
        infile, lctype)

    sel = ~( (time > 1348) & (time < 1350) )
    time, flux, err_flux = time[sel], flux[sel], err_flux[sel]

    flux = median_filter_lightcurve(time, flux, err_flux, filter_savpath)

    return time, flux, err_flux


def in_out_transit_plot(time, flux, intransit, ootransit, savpath):

    f, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(8,4))

    ax.scatter(time[ootransit], flux[ootransit], c='k', s=1.5, rasterized=True,
              linewidths=0)
    ax.scatter(time[intransit], flux[intransit], c='r', s=1.5, rasterized=True,
              linewidths=0)

    ax.set_ylabel('relative flux')
    ax.set_xlabel('time [days]')
    f.tight_layout(h_pad=0, w_pad=0)
    f.savefig(savpath, dpi=400, bbox_inches='tight')


if __name__ == "__main__":
    pb = 0.01
    plot_phased = False

    time, flux, err_flux = initial_wrangle()

    if plot_phased:

        period = 2.798581
        epoch = 2456297.719 - 2457000 # BJD --> BTJD 
        plot_phased_lcs(time, flux, period, epoch, substr='discovery', pb=pb)

        blsd, _ = run_two_periodograms(time, flux, err_flux, period-0.02,
                                       period+0.02,
                                       outdir='../results/tic_238176110_winn/',
                                       outname='two_periodograms.png',
                                       magsarefluxes=True, sigclip=None)

        blsd = kbls.bls_stats_singleperiod(time, flux, err_flux,
                                           blsd['bestperiod'],
                                           magsarefluxes=True, sigclip=None,
                                           perioddeltapercent=5)

        plot_phased_lcs(time, flux, blsd['period'], blsd['epoch'],
                        substr='blsparams', pb=pb)

        plot_phased_lcs(time, flux, blsd['period'], blsd['epoch'],
                        substr='blsparams_smaller', pb=0.005, phaselim=[-0.25,0.25])

    import pandas as pd
    df = pd.DataFrame({'time':time, 'flux':flux, 'err_flux':err_flux})
    df.to_csv('../results/tic_238176110_winn/lightcurve.csv', index=False)

    # do masking just to have it implemented
    tmids_obsd, t_starts, t_ends = given_lc_get_transit_tmids_tstarts_tends(
        time, flux, err_flux,
        blsfit_savpath='../results/tic_238176110_winn/blsfit_0.png',
        trapfit_savpath='../results/tic_238176110_winn/trapfit_0.png',
        magsarefluxes=True, nworkers=8, sigclip=None, N=0.03)

    in_transit = np.zeros_like(time).astype(bool)

    for t_start, t_end in zip(t_starts, t_ends):

        this_transit = ( (time > t_start) & (time < t_end) )

        in_transit |= this_transit

    out_of_transit = ~in_transit

    in_out_transit_plot(time, flux, in_transit, out_of_transit,
                        '../results/tic_238176110_winn/in_out_transit.png')

    # rerun BLS and GLS
    time, flux, err_flux = (time[out_of_transit], flux[out_of_transit],
                            err_flux[out_of_transit])

    _, _ = run_two_periodograms(time, flux, err_flux, 0.1,
                                (np.max(time)-np.min(time))/2,
                                outdir='../results/tic_238176110_winn/',
                                outname='bls_gls_with_hj_masked.png',
                                magsarefluxes=True, sigclip=None, autofreq=True)

    import IPython; IPython.embed()
