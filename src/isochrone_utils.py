# -*- coding: utf-8 -*-
'''
given spectroscopic (or photometric) parameters, see whether isochrones tell
you anything useful
'''
from __future__ import division, print_function

import numpy as np, pandas as pd, matplotlib.pyplot as plt

from astrobase.services.tic import tic_single_object_crossmatch
import astropy.units as u

def get_star_params_tic():

    try:
        result = tic_single_object_crossmatch(336.40231, -34.90962,
                                              (1*u.arcsec).to(u.deg).value)
        df = pd.DataFrame(result['data'][0], index=[0])
        outname = '../results/tic_62483237_maybe_young/62483237_tic_properties.csv'
        df.to_csv(outname, index=False)
        print('saved {:s}'.format(outname))

    except:
        try:
            outname = '../results/tic_62483237_maybe_young/62483237_tic_properties.csv'
            df = pd.read_csv(outname)
        except Exception as e:
            print(e)

    _teff = np.array(df['Teff'])
    e_teff = np.array(df['e_Teff'])
    logg = np.array(df['logg'])
    e_logg = np.array(df['e_logg'])
    mh = np.array(df['MH'])
    e_mh = np.array(df['e_MH'])
    if mh==None:
        mh = np.array(0)
    if e_mh==None:
        e_mh = np.array(0.1)

    Vmag = np.array(df['Vmag'])
    e_Vmag = np.array(df['e_Vmag'])
    Bmag = np.array(df['Bmag'])
    e_Bmag = np.array(df['e_Bmag'])
    Bmag = np.array(df['Bmag'])
    e_Bmag = np.array(df['e_Bmag'])

    #spectroscopic properties (value, uncertainty)

    Teff = (float(_teff), float(e_teff))
    logg = (float(logg), float(e_logg))
    feh = (float(mh), float(e_mh))
    parallax = (23.5527, 0.0421) # in mas -- from Gaia DR2

    print(
        'from TIC,'+
        '\n\tTeff = {:s}'.format(repr(Teff))+
        '\n\tlogg = {:s}'.format(repr(logg))+
        '\n\tfeh = {:s}'.format(repr(feh))+
        '\n\tparallax = {:s}'.format(repr(parallax))
    )

    return Teff, logg, feh, parallax, df


def get_star_params_RAVE():
    Teff = (4454.4, 67)
    logg = (4.89, 0.11)
    feh = (0.04, 0.08)
    parallax = (23.5527, 0.0421) # in mas -- from Gaia DR2

    print(
        'from RAVE,'+
        '\n\tTeff = {:s}'.format(repr(Teff))+
        '\n\tlogg = {:s}'.format(repr(logg))+
        '\n\tfeh = {:s}'.format(repr(feh))+
        '\n\tparallax = {:s}'.format(repr(parallax))
    )

    return Teff, logg, feh, parallax

