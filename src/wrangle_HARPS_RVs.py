# -*- coding: utf-8 -*-
'''
There is a TOI you want ESO/HARPS archival data for. Go to

    http://archive.eso.org/wdb/wdb/adp/phase3_spectral/form

select all the observations
download the shell script to e.g.,

    /home/luke/local/HARPS_data/tic_62483237_harps

execute it. The code in this file:

    * extracts the HARPS tars
    * moves things into a better directory structure
    * gets the raw RVs
    * plots the RVs
    * performs a maximum likelihood keplerian fit, given the TESS epoch and
      period
'''
from __future__ import division, print_function

import numpy as np, matplotlib.pyplot as plt, pandas as pd
from glob import glob
import os, tarfile

import radvel
from radvel.plot import orbit_plots, mcmc_plots
from scipy import optimize

import astropy.units as u
from astropy.io import ascii, fits

from numpy import array as nparr
from scipy.interpolate import interp1d

from astrobase.services.tic import tic_single_object_crossmatch

from read_harps import read_spec

############
# OUTDATED #
############
def _read_raw_harps():

    pattern = '*.fits'

    files = np.sort(glob(pattern))

    hdulist = fits.open(files[0])

    spec_hdr = hdulist[1].header
    spec_data = hdulist[1].data

    pipe_hdr = hdulist[0].header

def _search_headers_for_RVs():

    searchstrs = ['DVRMS', 'RADVEL', 'VELOCITY', 'RADIAL', 'CCF RVC', 'RVC']

    pattern = 'data/reduced/*e2ds*.fits'
    files = np.sort(glob(pattern))
    hdulist = fits.open(files[0]) #FIXME
    e2ds_hdr = hdulist[0].header
    print(hdulist.info())
    for searchstr in searchstrs:
        print([l for l in repr(e2ds_hdr).split('\n') if searchstr in l])

    pattern = 'data/reduced/*bis*.fits'
    files = np.sort(glob(pattern))
    hdulist = fits.open(files[0]) #FIXME
    bis_hdr = hdulist[0].header
    print(hdulist.info())
    for searchstr in searchstrs:
        print([l for l in repr(bis_hdr).split('\n') if searchstr in l])

    pattern = 'data/reduced/*s1d*.fits'
    files = np.sort(glob(pattern))
    hdulist = fits.open(files[0]) #FIXME
    s1d_hdr = hdulist[0].header
    print(hdulist.info())
    for searchstr in searchstrs:
        print([l for l in repr(bis_hdr).split('\n') if searchstr in l])


##########
# USEFUL #
##########
def extract_tars(tardir='/home/luke/local/HARPS_data/tic_62483237_harps/'):

    tarfiles = np.sort(glob(tardir+'*.tar'))

    for tf in tarfiles:
        tar = tarfile.open(tf)
        tar.extractall(path=tardir)
        tar.close()
        print('extracted %s to %s' % (tf, tardir))


def fix_file_structure(tardir, extracted_dir):
    # HARPS data from different times goes to different places. After ~2012 it
    # automatically extracts to tardir+'/data/reduced/20??-??-??' 
    # Before ~2012 it's automatically placed into tardir.

    if not os.path.exists(extracted_dir):
        os.mkdir(extracted_dir)

    newpattern = tardir+'data/reduced/20??-??-??/HARPS*.fits'
    new_harps_files = np.sort(glob(newpattern))
    oldpattern = tardir+'HARPS*.fits'
    old_harps_files = np.sort(glob(oldpattern))

    # Move the post-2012 data into tar_extracted
    for srcpath in np.concatenate((new_harps_files,old_harps_files)):
        dstpath = extracted_dir+os.path.basename(srcpath)
        os.rename(srcpath, dstpath)
        print('moved %s\n\tto %s' % (srcpath, dstpath))


def get_rvs(extracted_dir):
    '''
    get RVs, errors, times, exposure times, SNR, bisector width, and FWHM of
    CCF RVs from archival HARPS measurements.
    '''

    ccfpattern = extracted_dir+'*_ccf_*.fits'
    ccffiles = np.sort(glob(ccfpattern))

    v_p, sig_ccf_rv, v_p_sig, bjd, exptime, snr = [],[],[],[],[],[]

    for ccffile in ccffiles:
        hdulist = fits.open(ccffile)
        ccf_hdr = hdulist[0].header

        # barycentric RV, drift corrected [km/s]
        v_p.append(ccf_hdr['HIERARCH ESO DRS CCF RVC'])
        # photon noise in CCF RV [km/s]
        sig_ccf_rv.append(ccf_hdr['HIERARCH ESO DRS CCF NOISE'])
        # estimated RV uncertainties [m/s]
        v_p_sig.append(ccf_hdr['HIERARCH ESO DRS DVRMS'])
        # barycentric julian date
        bjd.append(ccf_hdr['HIERARCH ESO DRS BJD'])
        # total integration time
        exptime.append(ccf_hdr['EXPTIME'])
        # S_N order center48. there are 72 of these, so I'm not sure why one
        # would select "48".
        snr.append(ccf_hdr['HIERARCH ESO DRS SPE EXT SN48'])

        hdulist.close()

    v_p, sig_ccf_rv, v_p_sig, bjd, exptime, snr = (
        nparr(v_p), nparr(sig_ccf_rv), nparr(v_p_sig), nparr(bjd),
        nparr(exptime), nparr(snr)
    )

    v_p *= 1e3 # get v_p in meters per second
    sig_ccf_rv *= 1e3 # ditto for the photon noise in CCF RV

    # it's unclear whether the "estimated RV uncertainties" or the "photon
    # noise in the CCF RV" is the relevant error to use. Take the maximum of
    # the two (generally, this seems to imply the photon noise in the CCF RV).
    err_v_p = np.maximum(v_p_sig, sig_ccf_rv)

    if not np.all(sig_ccf_rv > v_p_sig):
        raise AssertionError('expected errors to be from photon noise in the'
                             'CCF RV. keep this as the general case.')

    # retrieve the activity indicators: the bisector span and the FWHM.
    bispattern = extracted_dir+'*_bis_*.fits'
    bisfiles = np.sort(glob(bispattern))

    biswidth, fwhm, ccf_mask = [],[],[]
    for bisfile in bisfiles:
        hdulist = fits.open(bisfile)
        bis_hdr = hdulist[0].header

        # bisector velocity span [km/s]
        biswidth.append(bis_hdr['HIERARCH ESO DRS BIS SPAN'])
        # FWHM of the CCF [km/s]
        fwhm.append(bis_hdr['HIERARCH ESO DRS CCF FWHM'])
        # mask filename to CCF mask
        ccf_mask.append(bis_hdr['HIERARCH ESO DRS CCF MASK'])

        hdulist.close()

    biswidth, fwhm, ccf_mask = (
        nparr(biswidth), nparr(fwhm), nparr(ccf_mask)
    )

    biswidth *= 1e3 # convert to m/s
    fwhm *= 1e3 # convert to m/s

    return v_p, err_v_p, bjd, exptime, snr, biswidth, fwhm, ccf_mask


def plot_rvs(rv, err_rv, bjd, savname):

    f,ax = plt.subplots(figsize=(8,6))

    t0 = int(np.round( np.min(bjd), 0 ) )

    ax.errorbar(bjd-t0, rv, yerr=err_rv,
                elinewidth=0.3, capsize=1, capthick=1, linewidth=0, fmt='s')

    ax.set_xlabel('BJD - {:d}'.format(t0))
    ax.set_ylabel('barycentric RV, drift corrected [m/s]')

    f.tight_layout()
    f.savefig(savname, dpi=350)
    print('made {:s}'.format(savname))


def plot_rvs_and_model(rv, err_rv, bjd, savname, likelihood, model_times,
                       xlim=None):

    f,ax = plt.subplots(figsize=(8,6))

    t0 = int(np.round( np.min(bjd), 0 ) )

    ax.errorbar(bjd-t0, rv, yerr=err_rv,
                elinewidth=0.3, capsize=1, capthick=1, linewidth=0, fmt='s')

    ax.plot(model_times-t0, likelihood.model(model_times))

    ax.set_xlabel('BJD - {:d}'.format(t0))
    ax.set_ylabel('barycentric RV, drift corrected [m/s]')

    if isinstance(xlim,list):
        ax.set_xlim(xlim)

    f.tight_layout()
    f.savefig(savname, dpi=350)
    print('made {:s}'.format(savname))


def initialize_rv_model(n_planets, time_base, per1, tc1, secosw1, sesinw1,
                        logk1, dvdt, curv):

    params = radvel.Parameters(n_planets ,basis='per tc secosw sesinw logk')

    params['per1'] = radvel.Parameter(value=per1)
    params['tc1'] = radvel.Parameter(value=tc1)
    params['secosw1'] = radvel.Parameter(value=secosw1)
    params['sesinw1'] = radvel.Parameter(value=sesinw1)
    params['logk1'] = radvel.Parameter(value=logk1)
    mod = radvel.RVModel(params, time_base=time_base)
    mod.params['dvdt'] = radvel.Parameter(value=dvdt)
    mod.params['curv'] = radvel.Parameter(value=curv)

    return mod


def _get_WM14_mass(Rp):
    # Rp given in earth radii. (float)
    # Mp returned in earth masses (float)
    # Weiss & Marcy 2014, Weiss+ 2013.
    # see Eqs 5,6,7,8 from Weiss+ 2018, CKS VI

    R_p = Rp*u.Rearth

    if R_p < 1.5*u.Rearth:
        ρ_p = (2.43 + 3.39*(R_p.to(u.Rearth).value))*(u.g/(u.cm**3))
        M_p = (ρ_p/(5.51*(u.g/(u.cm**3)))) * \
                (R_p.to(u.Rearth).value)**3 * u.Mearth

    elif R_p >= 1.5*u.Rearth and R_p <= 4.*u.Rearth:
        M_p = 2.69*(R_p.to(u.Rearth).value)**(0.93) * u.Mearth

    elif R_p > 4.*u.Rearth and R_p < 9.*u.Rearth:
        M_p = 0.86*(R_p.to(u.Rearth).value)**(1.89) * u.Mearth

    elif R_p > 9.*u.Rearth:
        M_p = 318*u.Mearth

    return float(M_p.to(u.Mearth).value)


def radvel_max_likelihood(rv, err_rv, bjd):

    # parameters from the TESS alert
    period = 11.05842 # days
    epoch = 1334.8964 # BTJD
    epoch += 2457000  # now BJD
    logk1 = 1 # 2.71 m/s as initial guess for semiamplitude.
    dvdt = 0. # slope
    curv = 0. # curvature
    secosw, sesinw = 0., 0.

    # for the prior on semiamplitude: assume the TESS planet is real. convert
    # from radius to mass using e.g., the Weiss relations.
    rp_by_rstar = 0.0333
    rp = 2.69 # r_earth
    mp_guess = _get_WM14_mass(rp)*u.Mearth
    sini = 1

    # guess the expected semiamplitude by using stellar parameters from TIC and
    # planet parameters from the alert. (and a planet mass-radius relation).
    result = tic_single_object_crossmatch(336.40231, -34.90962,
                                          (1*u.arcsec).to(u.deg).value)
    if not len(result['data'])==1:
        raise NotImplementedError('assumes you match only a single object')
    mstar = result['data'][0]['mass']*u.Msun
    rstar = result['data'][0]['rad']
    k1_guess = (
        (28.4329*u.m/u.s) * mp_guess * sini / (1*u.Mjup)
        * ( (mp_guess + mstar) / (1*u.Msun) )**(-2/3)
        * ( (period*u.day)/(1*u.year) )**(-1/3)
    )
    print('k1_guess [m/s]: {:.2f}'.format(k1_guess.to(u.m/u.s)))
    k1_guess = k1_guess.to(u.m/u.s).value

    # add priors
    min_logk1, max_logk1 = np.log(0.1*k1_guess), np.log(10*k1_guess)

    # white noise jitter [m/s]
    jitter = 1.0
    # offset between absolute velocity of the star and 0 km/s. AKA the "mean
    # center of mass velocity", from the relative velocity vectors in the
    # galaxy of the sun and the star.
    gamma = np.nanmean(rv)

    t_model = np.linspace( np.min(bjd)-5, np.max(bjd)+5, 10000)

    # abscissa for slope and curvature terms. should be near the midpoint of
    # the time baseline.
    time_base = int(np.nanmedian(t_model))

    # assume circular orbit, fixed period and time of transit
    mod = initialize_rv_model(1, time_base, period, epoch, secosw,
                              sesinw, logk1, dvdt, curv)

    likelihood = radvel.likelihood.RVLikelihood(mod, bjd, rv, err_rv)

    likelihood.params['gamma'] = radvel.Parameter(value=gamma)
    likelihood.params['jit'] = radvel.Parameter(value=jitter)

    # choose which parameters to vary/fix. by default, they all vary.
    # options are:
    # per1, tc1, secosw1, sesinw1, logk1, dvdt (linear term), curv (quadratic)
    # gamma, jit

    likelihood.params['per1'].vary = False
    likelihood.params['tc1'].vary = False
    likelihood.params['secosw1'].vary = False
    likelihood.params['sesinw1'].vary = False
    likelihood.params['logk1'].vary = True
    likelihood.params['dvdt'].vary = False
    likelihood.params['curv'].vary = False

    print(likelihood)

    # plot initial model
    initmodel_name = ('../results/sector_1_tois/'
                      'tic_62483237_harps_rv_initial_model.png')
    plot_rvs_and_model(rv-gamma, err_rv, bjd, initmodel_name, likelihood, t_model)

    # initialize posterior object; add priors
    posterior = radvel.posterior.Posterior(likelihood)

    posterior.priors += [radvel.prior.Gaussian('jit', np.log(3), 0.5)]
    posterior.priors += [radvel.prior.HardBounds('logk1', min_logk1, max_logk1)]

    # maximize the likelihood, print the updated posterior

    results = optimize.minimize(
        posterior.neglogprob_array,
        posterior.get_vary_params(),
        method='Powell'
    )

    # plot maximum likelihood model
    mlmodel_name = ('../results/sector_1_tois/'
                      'tic_62483237_harps_rv_maxlikelihood_model_full.png')
    plot_rvs_and_model(rv-gamma, err_rv, bjd, mlmodel_name, likelihood, t_model)
    # plot maximum likelihood model with narrow time limits
    mlmodel_name = ('../results/sector_1_tois/'
                      'tic_62483237_harps_rv_maxlikelihood_model_xlimcut1.png')
    plot_rvs_and_model(rv-gamma, err_rv, bjd, mlmodel_name, likelihood,
                       t_model, xlim=[980,1220])
    mlmodel_name = ('../results/sector_1_tois/'
                      'tic_62483237_harps_rv_maxlikelihood_model_xlimcut2.png')
    plot_rvs_and_model(rv-gamma, err_rv, bjd, mlmodel_name, likelihood,
                       t_model, xlim=[600,800])
    mlmodel_name = ('../results/sector_1_tois/'
                      'tic_62483237_harps_rv_maxlikelihood_model_xlimcut3.png')
    plot_rvs_and_model(rv-gamma, err_rv, bjd, mlmodel_name, likelihood,
                       t_model, xlim=[-20,100])

    print(posterior)

def _K_channel_response(wvlen):
    # arg: wavelength vector
    # return response function to multiply the flux

    response = np.zeros_like(wvlen)

    Ca_K_line = 3933.66
    K_channel_min = Ca_K_line - 2*1.09
    K_channel_max = Ca_K_line + 2*1.09

    x = np.array([ K_channel_min, Ca_K_line, K_channel_max ])
    y = np.array([ 0, 1, 0 ])

    fn = interp1d(x, y, kind='linear', bounds_error=False, fill_value=0)

    response = fn(wvlen)

    return response

def _H_channel_response(wvlen):
    # arg: wavelength vector
    # return response function to multiply the flux

    response = np.zeros_like(wvlen)

    Ca_H_line = 3968.47
    H_channel_min = Ca_H_line - 2*1.09
    H_channel_max = Ca_H_line + 2*1.09

    x = np.array([ H_channel_min, Ca_H_line, H_channel_max ])
    y = np.array([ 0, 1, 0 ])

    fn = interp1d(x, y, kind='linear', bounds_error=False, fill_value=0)

    response = fn(wvlen)

    return response

def C_cf(T_eff):
    return 10**( (-1.7e-7)*T_eff**2 + (2.25e-3)*T_eff - 7.31 )

def R_phot(T_eff):
    return 10**(-4.78845 - 3.707 / (1 + (T_eff/4587.82)**17.5272))

def calculate_R_HK_prime(extracted_dir, T_eff=4454.4):
    '''
    R_HK := flux in Ca H & and K emission lines / bolometric flux.

    R_HK_prime = R_HK - correction for star's chromospheric flux.

    For details, see ../doc/181001_R_HK_prime_details.txt

    Teff = 4454.4 K is TIC 62483237's RAVE Teff.
    '''

    # build the response function for the triangle notches
    specpattern = extracted_dir+'*_s1d_*.fits'
    specfiles = np.sort(glob(specpattern))

    for specfile in specfiles:

        wvlen, flux = read_spec(specfile)

        R_channel = (4011.07 > wvlen) & (wvlen > 3991.07)
        V_channel = (3911.07 > wvlen) & (wvlen > 3891.07)

        R = np.sum(flux[R_channel])
        V = np.sum(flux[V_channel])

        H_response = _H_channel_response(wvlen)
        K_response = _K_channel_response(wvlen)

        H = np.sum(flux*H_response)
        K = np.sum(flux*K_response)

        #NOTE are you SURE this shouldn't be actually integrated, with a d-lambda?

        # follow Lorenzo-Oliviera et al (2018) for all the response
        S_HARPS = 18.349 * (H+K)/(R+V)
        S_MW = 0.9444 * S_HARPS + 0.0475
        R_HK = 1.34e-4 * C_cf(T_eff) * S_MW
        R_HK_prime = R_HK - R_phot(T_eff)

        print('S_HARPS: {:.2e}'.format(S_HARPS))
        print('S_MW: {:.2e}'.format(S_MW))
        print('R_HK_prime: {:.2e}'.format(R_HK_prime))
        print('log10(R_HK_prime): {:.2f}'.format(np.log10(R_HK_prime)))


def plot_spectra(extracted_dir):
    '''
    The "s1d" file is a 1-dimensional extracted spectrum which has been
    corrected for barycentric motion and rebinned to 0.1 Ang steps.

    Plots:
        * the full spectrum.
        * Ca K 3934
        * Ca H 3969
        * Lithium I resonance doublet: 6707.76 and 6707.91 A.
        * a blended Fe I line to Li I at 6707.44 A.

    '''

    specpattern = extracted_dir+'*_s1d_*.fits'
    specfiles = np.sort(glob(specpattern))

    for specfile in specfiles:

        wvlen, flux = read_spec(specfile)

        savdir = '../results/harps_1d_spectra/tic_62483237_toi_139/'

        plt.close('all')
        f,ax = plt.subplots(figsize=(12,4))
        ax.plot(wvlen, flux)
        ax.set_xlabel('wavelength [Angstr]')
        ax.set_ylabel('relative flux')
        savname = os.path.basename(specfile).replace('.fits','_full.png')
        f.savefig(savdir+savname, dpi=350)
        print('saved %s' % savname)

        # Ca K: 3934 Angstrom
        # Ca H: 3969 Angstrom
        plt.close('all')
        f,axs = plt.subplots(nrows=1, ncols=2, figsize=(12,4))
        for ax in axs:
            ax.plot(wvlen, flux)
            ax.set_xlabel('wavelength [Angstr]')
            ax.set_ylabel('relative flux')

        CaK_wvlen = 3934
        CaH_wvlen = 3969
        delta_wvlen = 50

        axs[0].set_xlim([CaK_wvlen-delta_wvlen, CaK_wvlen+delta_wvlen])
        axs[1].set_xlim([CaH_wvlen-delta_wvlen, CaH_wvlen+delta_wvlen])

        CaK_sel = (
            (wvlen > CaK_wvlen-delta_wvlen) & (wvlen < CaK_wvlen+delta_wvlen)
        )
        CaH_sel = (
            (wvlen > CaH_wvlen-delta_wvlen) & (wvlen < CaH_wvlen+delta_wvlen)
        )

        axs[0].set_title('Ca K 3934 angstr')
        axs[1].set_title('Ca H 3969 angstr')

        axs[0].vlines(CaK_wvlen, np.min(flux[CaK_sel])-10,
                      1.02*np.max(flux[CaK_sel]), color='r',zorder=-1)
        axs[1].vlines(CaH_wvlen, np.min(flux[CaH_sel])-10,
                      1.02*np.max(flux[CaH_sel]), color='r',zorder=-1)

        axs[0].set_ylim([np.min(flux[CaK_sel])-10, 1.02*np.max(flux[CaK_sel])])
        axs[1].set_ylim([np.min(flux[CaH_sel])-10, 1.02*np.max(flux[CaH_sel])])

        savname = os.path.basename(specfile).replace('.fits','_CaHK.png')
        f.savefig(savdir+savname, dpi=350)
        print('saved %s' % savname)

        # Lithium I resonance doublet:
        # one transition at 6707.76 and other at 6707.91 A.
        # note there is a Fe I line at 6707.44 A.

        plt.close('all')
        f,axs = plt.subplots(nrows=1, ncols=2, figsize=(12,4))
        for ax in axs:
            ax.plot(wvlen, flux)
            ax.set_xlabel('wavelength [Angstr]')
            ax.set_ylabel('relative flux')

        Li_doublet_wvlen = (6707.76 + 6707.91)/2
        FeI_line = 6707.44
        delta_wvlen = 10

        axs[0].set_xlim([Li_doublet_wvlen-delta_wvlen,
                         Li_doublet_wvlen+delta_wvlen])
        axs[1].set_xlim([FeI_line-delta_wvlen,
                         FeI_line+delta_wvlen])

        Li_doublet_sel = (
            (wvlen > Li_doublet_wvlen-delta_wvlen)
            &
            (wvlen < Li_doublet_wvlen+delta_wvlen)
        )
        FeI_sel = (
            (wvlen > FeI_line-delta_wvlen) & (wvlen < FeI_line+delta_wvlen)
        )

        axs[0].set_title('Li doublet (6707.76, 6707.91)')
        axs[1].set_title('FeI line (6707.44)')

        axs[0].vlines([(6707.76, 6707.91)], np.min(flux[Li_doublet_sel])-10,
                      1.02*np.max(flux[Li_doublet_sel]), color='r',zorder=-1)
        axs[1].vlines(FeI_line, np.min(flux[FeI_sel])-10,
                      1.02*np.max(flux[FeI_sel]), color='r',zorder=-1)

        axs[0].set_ylim([np.min(flux[Li_doublet_sel])-10,
                         1.02*np.max(flux[Li_doublet_sel])])
        axs[1].set_ylim([np.min(flux[FeI_sel])-10, 1.02*np.max(flux[FeI_sel])])

        savname = os.path.basename(specfile).replace('.fits','_lithium.png')
        f.savefig(savdir+savname, dpi=350)
        print('saved %s' % savname)



if __name__=="__main__":

    # functionality
    run_plot_spectra = False

    # define paths
    tardir='/Users/luke/local/HARPS_data/tic_62483237_harps/'
    extracted_dir=tardir+'tar_extracted/'

    resultdir = '../results/sector_1_tois/'
    rv_rawplot_savname=resultdir+'tic_62483237_harps_rv.png'
    rv_k5_mask_savname=resultdir+'tic_62483237_harps_rv_k5_mask_only.png'
    rv_df_savname=resultdir+'tic_62483237_harps_rv.csv'
    rv_df_k5mask_savname=resultdir+'tic_62483237_harps_rv_k5_mask_only.csv'

    # get the data
    extract_tars(tardir)
    fix_file_structure(tardir, extracted_dir)

    if run_plot_spectra:
        plot_spectra(extracted_dir)

    calculate_R_HK_prime(extracted_dir)

    rv, err_rv, bjd, exptime, snr, biswidth, fwhm, ccf_mask = (
        get_rvs(extracted_dir)
    )

    # plot and save the rvs, depending on the CCF mask used.
    plot_rvs(rv, err_rv, bjd, rv_rawplot_savname)
    df = pd.DataFrame({'bjd':bjd, 'rv':rv, 'err_rv':err_rv,
                       'biswidth':biswidth, 'fwhm':fwhm, 'ccf_mask':ccf_mask})
    df.to_csv(rv_df_savname,index=False)
    print('saved %s' % rv_df_savname)
    del df

    k5mask = (ccf_mask == 'K5')
    plot_rvs(rv[k5mask], err_rv[k5mask], bjd[k5mask], rv_k5_mask_savname)
    df = pd.DataFrame({'bjd':bjd[k5mask], 'rv':rv[k5mask],
                       'err_rv':err_rv[k5mask], 'biswidth':biswidth[k5mask],
                       'fwhm':fwhm[k5mask], 'ccf_mask':ccf_mask[k5mask]})
    df.to_csv(rv_df_k5mask_savname,index=False)
    print('saved %s' % rv_df_k5mask_savname)
    del df

    # now perform a maximum likelihood keplerian fit, assuming circular orbit
    # and a fixed period and time of transit.
    radvel_max_likelihood(rv[k5mask], err_rv[k5mask], bjd[k5mask])
