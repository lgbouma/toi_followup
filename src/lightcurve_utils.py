# -*- coding: utf-8 -*-
from __future__ import division, print_function

import matplotlib as mpl
mpl.use('Agg')

import numpy as np, matplotlib.pyplot as plt
from astropy.io import fits

from astrobase import periodbase, checkplot
from astrobase import astrotess as at
from astrobase.periodbase import kbls
from astrobase.varbase import lcfit
from astrobase.varbase.trends import smooth_magseries_ndimage_medfilt
from astrobase import lcmath

from glob import glob
from parse import search
import os

from astrobase.services.tic import tic_single_object_crossmatch

def get_a_over_Rstar_guess(lcfile, period):
    # xmatch TIC. get Mstar, and Rstar.
    # with period and Mstar, you can guess the semimajor axis.
    # then calculate a/Rstar.

    # get object RA/dec, so that you can get the Teff/logg by x-matching TIC,
    # so that you can get the theoretical limb-darkening coefficients.
    hdulist = fits.open(lcfile)
    main_hdr = hdulist[0].header
    lc_hdr = hdulist[1].header
    lc = hdulist[1].data
    ra, dec = lc_hdr['RA_OBJ'], lc_hdr['DEC_OBJ']
    sep = 0.1*u.arcsec
    obj = tic_single_object_crossmatch(ra,dec,sep.to(u.deg).value)
    if len(obj['data'])==1:
        rad = obj['data'][0]['rad']
        mass = obj['data'][0]['mass']
    else:
        raise NotImplementedError
    if not isinstance(rad ,float) and not isinstance(mass, float):
        raise NotImplementedError

    # P^2 / a^3   = 4pi^2 / (GM)
    # a^3 = P^2 * GM / (4pi^2)
    # a = ( P^2 * GM / (4pi^2) )^(1/3)

    P = period*u.day
    M = mass*u.Msun
    a = ( P**2 * const.G*M / (4*np.pi**2) )**(1/3)
    R = rad*u.Rsun
    a_by_Rstar = (a.cgs/R.cgs).value

    return a_by_Rstar


def get_limb_darkening_initial_guesses(lcfile):
    '''
    CITE: Claret 2017, whose coefficients we're parsing
    '''

    # get object RA/dec, so that you can get the Teff/logg by x-matching TIC,
    # so that you can get the theoretical limb-darkening coefficients.
    hdulist = fits.open(lcfile)
    main_hdr = hdulist[0].header
    lc_hdr = hdulist[1].header
    lc = hdulist[1].data
    ra, dec = lc_hdr['RA_OBJ'], lc_hdr['DEC_OBJ']
    sep = 0.1*u.arcsec
    obj = tic_single_object_crossmatch(ra,dec,sep.to(u.deg).value)
    if len(obj['data'])==1:
        teff = obj['data'][0]['Teff']
        logg = obj['data'][0]['logg']
        metallicity = obj['data'][0]['MH'] # often None
        if not isinstance(metallicity,float):
            metallicity = 0 # solar
    else:
        raise NotImplementedError

    # get the Claret quadratic priors for TESS bandpass
    # the selected table below is good from Teff = 1500 - 12000K, logg = 2.5 to
    # 6. We choose values computed with the "r method", see
    # http://vizier.u-strasbg.fr/viz-bin/VizieR-n?-source=METAnot&catid=36000030&notid=1&-out=text
    if not 2300 < teff < 12000:
        raise AssertionError
    if not 2.5 < logg < 6:
        raise AssertionError

    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    catalog_list = Vizier.find_catalogs('J/A+A/600/A30')
    catalogs = Vizier.get_catalogs(catalog_list.keys())
    t = catalogs[1]
    sel = (t['Type'] == 'r')
    df = t[sel].to_pandas()

    # since we're using these as starting guesses, not even worth
    # interpolating. just use the closest match!
    # each Teff gets 8 logg values. first, find the best teff match.
    foo = df.iloc[(df['Teff']-teff).abs().argsort()[:8]]
    # then, among those best 8, get the best logg match.
    bar = foo.iloc[(foo['logg']-logg).abs().argsort()].iloc[0]

    u_linear = bar['aLSM']
    u_quad = bar['bLSM']

    return float(u_linear), float(u_quad)

def make_LS_ACF_periodograms(infile):

    lctype = 'PDCSAP'
    time, flux, err_flux = at.get_time_flux_errs_from_Ames_lightcurve(infile, lctype)

    # determine number of cadences over which to smooth ACF.
    cadence_day = 2/60/24
    smooth_over = 1 # day
    smoothacf = np.round(smooth_over/cadence_day,0)
    if smoothacf % 2 == 0:
        smoothacf += 1
    smoothacf = int(smoothacf)

    # make the path
    outdir = '../results/sector_1_tois/'
    outpath = outdir+infile.split('/')[-1].replace(
        '.fits.gz','_{:s}_ACF_and_LSP.png'.format(lctype))

    gls = periodbase.pgen_lsp(time, flux, err_flux, magsarefluxes=True)
    acf = periodbase.macf_period_find(time, flux, err_flux,
                                      smoothacf=smoothacf,
                                      fillgaps='noiselevel',
                                      smoothfunckwargs={'polyorder':3},
                                      magsarefluxes=True, maxacfpeaks=5)

    cpf = checkplot.twolsp_checkplot_png(gls, acf, time, flux, err_flux,
                                         outfile=outpath, objectinfo=infodict)


def _calculate_panel_nrows_ncols(n_transits):
    if n_transits <= 4:
        nrows, ncols = 2, 2
    elif 4 < n_transits <= 8:
        nrows, ncols = 2, 4
    elif 8 < n_transits <= 12:
        nrows, ncols = 3, 4
    elif 12 < n_transits <= 16:
        nrows, ncols = 4, 4
    elif 16 < n_transits <= 20:
        nrows, ncols = 4, 5
    elif 20 < n_transits <= 25:
        nrows, ncols = 5, 5
    elif 25 < n_transits <= 30:
        nrows, ncols = 5, 6
    elif 30 < n_transits <= 36:
        nrows, ncols = 6, 6
    elif 36 < n_transits <= 42:
        nrows, ncols = 6, 7
    elif 42 < n_transits <= 48:
        nrows, ncols = 6, 8
    elif n_transits > 48:
        raise AssertionError('consider a different visualization approach')

    return nrows, ncols


def plot_individual_transit_inspection_panel(time, flux, err_flux,
                                             indivtransit_savpath=None,
                                             blsfit_savpath=None,
                                             trapfit_savpath=None,
                                             sma_guess=10,
                                             u_linear=0.2,
                                             u_quad=0.2
                                             ):
    '''
    Given a lightcurve, make a plot of the individual transits. You need
    individual transits to have SNR>~5 for this to be worth your time. E.g.,

    transit_001         transit_002             transit_003
    transit_004         transit_005             transit_006
    '''

    if not isinstance(indivtransit_savpath,str):
        raise AssertionError

    # first, run BLS to get an initial epoch and period.
    endp = 1.05*(np.nanmax(time) - np.nanmin(time))/2
    blsdict = kbls.bls_parallel_pfind(time, flux, err_flux, magsarefluxes=True,
                                      startp=0.1, endp=endp,
                                      maxtransitduration=0.3, nworkers=8,
                                      sigclip=15.)
    blsd = kbls.bls_stats_singleperiod(time, flux, err_flux,
                                       blsdict['bestperiod'],
                                       magsarefluxes=True, sigclip=15.,
                                       perioddeltapercent=5)
    #  plot the BLS model.
    lcfit._make_fit_plot(blsd['phases'], blsd['phasedmags'], None,
                         blsd['blsmodel'], blsd['period'], blsd['epoch'],
                         blsd['epoch'], blsfit_savpath, magsarefluxes=True)

    ingduration_guess = blsd['transitduration']*0.2
    transitparams = [blsd['period'], blsd['epoch'], blsd['transitdepth'],
                     blsd['transitduration'], ingduration_guess
                    ]

    # fit a trapezoidal transit model; plot the resulting phased LC.
    trapd = lcfit.traptransit_fit_magseries(time, flux, err_flux,
                                            transitparams, magsarefluxes=True,
                                            sigclip=15.,
                                            plotfit=trapfit_savpath)

    # use the trapezoidal model's epoch as the guess to identify (roughly) in
    # and out of transit points
    tmids, t_starts, t_ends = get_transit_times(blsd, time, 2, trapd=trapd)

    n_transits = len(tmids)
    nrows, ncols = _calculate_panel_nrows_ncols(n_transits)

    # make the panel plot
    plt.close('all')
    f, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8,6), sharex=False,
                         sharey=True)

    ymin, ymax = ( 0.998*np.nanmin(flux),
                   np.minimum(
                       np.nanmedian(flux)+2*np.nanstd(flux),
                       1.001 * np.nanmax(flux)
                   )
                 )

    for transit_ix, t_start, t_end, ax in list(
        zip(range(len(t_starts)), t_starts, t_ends, axs.flatten())):

        try:
            sel = (time < t_end) & (time > t_start)
            sel_time = time[sel]
            sel_flux = flux[sel]
            sel_err_flux = err_flux[sel]

            pctl = 0.2 # desired percentile
            ix_near = abs(sel_flux -
                          np.percentile(sel_flux,pctl,interpolation='nearest')
                         ).argmin()
            t_mid_guess = sel_time[ix_near]

            init_rp = np.sqrt(trapd['fitinfo']['finalparams'][2])
            fitparams = {'t0':t_mid_guess, 'rp':init_rp, 'sma':sma_guess,
                         'incl':85, 'u':[u_linear,u_quad]}
            fixedparams = {'ecc':0., 'omega':90., 'limb_dark':'quadratic',
                           'period':blsd['period'] }
            priorbounds = {'rp':(init_rp-0.05, init_rp+0.05),
                           'u_linear':(u_linear-1, u_linear+1),
                           'u_quad':(u_quad-1, u_quad+1),
                           't0':(t_mid_guess-0.1, t_mid_guess+0.1),
                           'sma':(0.7*sma_guess,1.3*sma_guess),
                           'incl':(75,90) }

            maf = lcfit.mandelagol_fit_magseries(
                sel_time, sel_flux, sel_err_flux, fitparams, priorbounds,
                fixedparams, sigclip=30, n_mcmc_steps=500,
                mcmcprogressbar=True, samplesavpath='../data/temp.h5',
                nworkers=8, overwriteexistingsamples=True
            )

            t_mid = maf['fitinfo']['fitepoch']
            fit_flux = maf['fitinfo']['fitmags']
            fit_time = maf['magseries']['times']

            ax.scatter(24*(sel_time-t_mid), sel_flux, c='k', alpha=0.8, s=4,
                       zorder=1, rasterized=True, linewidths=0)
            # don't show the fit; we just want a good midtime quickly.
            # ax.plot(24*(fit_time-t_mid), fit_flux, alpha=0.8, zorder=2, lw=1,
            #         rasterized=True)

            ax.vlines(0, ymin, ymax, colors='k', linestyles='solid', alpha=0.3,
                     zorder=-1, lw=1)

            ax.text(0.01,0.99,'t'+str(transit_ix).zfill(3),
                    ha='left',va='top',zorder=2,fontsize='xx-small',
                    transform=ax.transAxes)

            if transit_ix >= (nrows-1)*ncols:
                ax.set_xlabel('time from tmid [hr]',fontsize='xx-small')
            if transit_ix % ncols == 0:
                ax.set_ylabel('relative flux', fontsize='xx-small')

            ax.set_ylim([ymin,ymax])
            twidth = 24*min(
                abs(np.nanmin(sel_time)-t_mid),abs(np.nanmax(sel_time)-t_mid)
            )
            ax.set_xlim([-twidth,twidth])

            del maf

        except Exception as e:
            print(e)
            print('transit {:d} failed, continue'.format(transit_ix))
            continue

    for ax in axs.flatten():
        ax.get_yaxis().set_tick_params(which='both', direction='in')
        ax.get_xaxis().set_tick_params(which='both', direction='in')
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize('xx-small')
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize('xx-small')

    f.tight_layout(h_pad=0, w_pad=0)
    f.savefig(indivtransit_savpath, dpi=350, bbox_inches='tight')
    print('saved {:s}'.format(indivtransit_savpath))
    f.savefig(indivtransit_savpath.replace('.png','.pdf'), bbox_inches='tight')
    print('saved {:s}'.format(indivtransit_savpath.replace('.png','.pdf')))



def single_whitening_plot(time, flux, smooth_flux, whitened_flux, savpath):
    f, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8,6))
    axs[0].scatter(time, flux, c='k', alpha=0.5, label='PDCSAP', zorder=1,
                   s=1.5, rasterized=True, linewidths=0)
    axs[0].plot(time, smooth_flux, 'b-', alpha=0.9, label='median filter',
                zorder=2)
    axs[1].scatter(time, whitened_flux, c='k', alpha=0.5,
                   label='PDCSAP/median filtered',
                   zorder=1, s=1.5, rasterized=True, linewidths=0)

    for ax in axs:
        ax.legend(loc='best')

    axs[0].set(ylabel='relative flux')
    axs[1].set(xlabel='time [days]', ylabel='relative flux')
    f.tight_layout(h_pad=0, w_pad=0)
    f.savefig(savpath, dpi=400, bbox_inches='tight')


def median_filter_lightcurve(time, flux, err_flux, savpath, mingap_min=240.,
                             smooth_window_day=2., cadence_min=2,
                             make_diagnostic_plots=True):
    '''
    detrending parameters.
    mingap_min: minimum gap to determine time group size. (in minutes).
    smooth_window_day: window for median filtering.
    '''

    mingap = mingap_min/60/24
    cadence_day = cadence_min/60/24

    windowsize = int(smooth_window_day/cadence_day)
    if windowsize % 2 == 0:
        windowsize += 1

    # get time groups, and median filter each one
    ngroups, groups = lcmath.find_lc_timegroups(time, mingap=mingap)

    tg_smooth_flux = []
    for group in groups:
        tg_flux = flux[group]
        tg_smooth_flux.append(
            smooth_magseries_ndimage_medfilt(tg_flux, windowsize)
        )

    smooth_flux = np.concatenate(tg_smooth_flux)
    filtered_flux = flux/smooth_flux

    if make_diagnostic_plots:
        single_whitening_plot(time, flux, smooth_flux, filtered_flux, savpath)

    return filtered_flux


def make_individual_transit_inspection_panel(infile):
    '''
    wrapper to plot_individual_transit_inspection_panel that interfaces with
    TOI alert lightcurves
    '''

    lctype = 'PDCSAP'

    # make the paths
    basedir = '../results/individual_transit_inspection/'
    outdir = basedir+'sector_1_alerts/'
    for dirname in [basedir, outdir]:
        if not os.path.exists(dirname):
            os.mkdir(dirname)
    outpath = outdir+infile.split('/')[-1].replace(
        '.fits.gz','_{:s}_individual_transits.png'.format(lctype))
    blsfit_savpath = outdir+infile.split('/')[-1].replace(
        '.fits.gz','_{:s}_blsfit_phased.png'.format(lctype))
    trapfit_savpath = outdir+infile.split('/')[-1].replace(
        '.fits.gz','_{:s}_trapfit_phased.png'.format(lctype))
    filter_savpath = outdir+infile.split('/')[-1].replace(
        '.fits.gz','_{:s}_medianfilter.png'.format(lctype))

    time, flux, err_flux = at.get_time_flux_errs_from_Ames_lightcurve(
        infile, lctype)

    flux = median_filter_lightcurve(time, flux, err_flux, filter_savpath)

    plot_individual_transit_inspection_panel(time, flux, err_flux,
                                             indivtransit_savpath=outpath,
                                             blsfit_savpath=blsfit_savpath,
                                             trapfit_savpath=trapfit_savpath)

def run_two_periodograms(times, mags, errs, startp, endp,
                         outdir='../results/',
                         outname='foobar.png', magsarefluxes=False,
                         sigclip=None, autofreq=False):
    '''
    return blsdict, gls
    '''

    if not autofreq:
        blsdict = kbls.bls_parallel_pfind(times, mags, errs,
                                          magsarefluxes=magsarefluxes,
                                          startp=startp, endp=endp,
                                          maxtransitduration=0.3, nworkers=8,
                                          sigclip=sigclip, stepsize=1e-7,
                                          autofreq=False)
    elif autofreq:
        blsdict = kbls.bls_parallel_pfind(times, mags, errs,
                                          magsarefluxes=magsarefluxes,
                                          startp=startp, endp=endp,
                                          maxtransitduration=0.3, nworkers=8,
                                          sigclip=sigclip, autofreq=True)

    gls = periodbase.pgen_lsp(times, mags, errs, magsarefluxes=magsarefluxes,
                              startp=startp, endp=endp, nworkers=8,
                              sigclip=sigclip)
    outpath = outdir+outname
    cpf = checkplot.twolsp_checkplot_png(blsdict, gls, times, mags, errs,
                                         outfile=outpath, objectinfo=None)

    return blsdict, gls
