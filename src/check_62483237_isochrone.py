# -*- coding: utf-8 -*-
'''
given spectroscopic (or photometric) parameters, see whether isochrones tell
you anything useful
'''
from __future__ import division, print_function

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from glob import glob
import os, tarfile

from isochrones import StarModel
from isochrones.mist import MIST_Isochrone
from isochrones.priors import age_prior
print(age_prior.bounds)

import isochrone_utils as iu

if __name__ == '__main__':

    outcornerfile = '../results/tic_62483237_maybe_young/isochrones_corner.png'
    outsampleh5 = '../results/tic_62483237_maybe_young/starmodel.h5'

    prevchains = glob('chains/*')
    for prevchain in prevchains:
        print('rm %s' % prevchain)
        os.remove(prevchain)

    teff, logg, feh, parallax, tic_df = iu.get_star_params_tic()
    teff, logg, feh, parallax = iu.get_star_params_RAVE()
    rho = (1.693, 0.05) # from transit model of SPOC

    mist = MIST_Isochrone()

    model = StarModel(mist, Teff=teff, logg=logg, feh=feh, parallax=parallax,
                      density=rho)

    model.fit()

    model.corner_physical()
    plt.savefig(outcornerfile, dpi=400)

    model.save_hdf(outsampleh5, overwrite=True)
