import numpy as np
from astropy.io import fits

from astrobase import periodbase, checkplot

from glob import glob
from parse import search

def make_LS_ACF_periodograms(infile):

    hdulist = fits.open(infile)

    main_hdr = hdulist[0].header
    lc_hdr = hdulist[1].header
    lc = hdulist[1].data

    time = ['TIME']

    assert 'Ames' in main_hdr['ORIGIN']
    assert 'LIGHTCURVE' in lc_hdr['EXTNAME']

    infodict = {'objectid':main_hdr['OBJECT'],
                'ra':main_hdr['RA_OBJ'],
                'decl':main_hdr['DEC_OBJ']}

    # run Lomb-Scargle on SAP lightcurve
    outdir = '../results/sector_1_tois/'
    outpath = outdir+infile.split('/')[-1].replace(
        '.fits.gz','_SAP_ACF_and_LSP.png')
    time, flux, err_flux = lc['TIME'], lc['SAP_FLUX'], lc['SAP_FLUX_ERR']

    # ensure 2 minute cadence
    assert np.abs(np.nanmedian(np.diff(time))*24*60 - 2) < 1e-2
    cadence_day = 2/60./24.
    smooth_over = 1 # day
    smoothacf = np.round(smooth_over/cadence_day,0)
    #smoothacf /= 5
    if smoothacf % 2 == 0:
        smoothacf += 1
    smoothacf = int(smoothacf)

    # gls = periodbase.pgen_lsp(time, flux, err_flux,
    #                           magsarefluxes=True)
    # acf = periodbase.macf_period_find(time, flux, err_flux,
    #                                   smoothacf=smoothacf,
    #                                   fillgaps='noiselevel',
    #                                   smoothfunckwargs={'polyorder':3},
    #                                   magsarefluxes=True, maxacfpeaks=5)
    # cpf = checkplot.twolsp_checkplot_png(gls, acf, time, flux, err_flux,
    #                                      outfile=outpath, objectinfo=infodict)
    # run Lomb-Scargle on PDC-SAP lightcurve
    outpath = outdir+infile.split('/')[-1].replace(
        '.fits.gz','_PDCSAP_ACF_and_LSP.png')
    time, flux, err_flux = lc['TIME'], lc['PDCSAP_FLUX'], lc['PDCSAP_FLUX_ERR']

    gls = periodbase.pgen_lsp(time, flux, err_flux, magsarefluxes=True)
    acf = periodbase.macf_period_find(time, flux, err_flux,
                                      smoothacf=smoothacf,
                                      fillgaps='noiselevel',
                                      smoothfunckwargs={'polyorder':3},
                                      magsarefluxes=True, maxacfpeaks=5)
    cpf = checkplot.twolsp_checkplot_png(gls, acf, time, flux, err_flux,
                                         outfile=outpath, objectinfo=infodict)

def check_62483237():

    lcdir = '../data/alert_lightcurves/'
    fname = glob(lcdir+'*62483237*')
    assert len(fname) == 1
    fname = fname[0]

    make_LS_ACF_periodograms(fname)

if __name__ == "__main__":
    check_62483237()
