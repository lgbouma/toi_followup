# -*- coding: utf-8 -*-
'''
crossmatch lists of known planets that i expect to be observed with what was
actually alerted.
'''

from __future__ import division, print_function

import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use('Agg')

import os
from glob import glob

def main(
    sector_number,
    mastcsvpath='../data/Kane_MAST_Crossmatch_CTL.csv',
    alertcsvpath='../data/toi-2018-12-15.csv',
    knownplanetdir='/home/luke/Dropbox/proj/tessmaps/results'
):

    exo_df = pd.read_csv(mastcsvpath)
    alert_df = pd.read_csv(alertcsvpath)

    mdf = alert_df.merge(exo_df[['MatchID','pl_hostname']], how='left',
                         left_on='tic_id', right_on='MatchID')

    # known planets if they had a name in exoplanet archive
    alertedknown_df = mdf[['pl_hostname','tic_id']].dropna(subset=['pl_hostname'])
    alertedknown_df = alertedknown_df.drop_duplicates()

    print('known planets that have been alerted\n')
    print(alertedknown_df)

    # sector_number: 1-based. (!!!)


    kpglob_sectornumber = sector_number - 1
    knownplanetglob = (
        'kane_knownplanets_sector{:d}.csv'.format(kpglob_sectornumber)
    )

    knownplanetfiles = glob(os.path.join(knownplanetdir, knownplanetglob))
    if len(knownplanetfiles) != 1:
        raise AssertionError('something wrong in glob')
    knownplanetfile = knownplanetfiles[0]
    knownthissector_df = pd.read_csv(knownplanetfile)

    print('\noverlap with planets expected to be observed in sector {:d}\n'.
         format(sector_number))

    knownthissector_df = knownthissector_df.merge(alertedknown_df,
                                                  on='pl_hostname', how='left')

    knownthissectorandalerted = knownthissector_df.dropna(subset=['tic_id'])
    knownthissectorandalerted['tic_id'] = (
        np.array(knownthissectorandalerted['tic_id']).astype(int)
    )

    knownthissector_transiting_notalerted = (
        knownthissector_df[(knownthissector_df['is_transiting']==1) &
                           (np.isnan(knownthissector_df['tic_id']))]
    )

    print('known in sector {:d} and alerted:'.format(sector_number))
    print(knownthissectorandalerted)

    print('known in sector {:d}, transiting, and not alerted:'.format(sector_number))
    print(knownthissector_transiting_notalerted)

    print('the corresponding TICIDs for known in sector {:d}, transiting, and not alerted:'.
          format(sector_number)
         )

    _ktna = pd.DataFrame(
        {'name':np.array(knownthissector_transiting_notalerted['pl_hostname'])}
    )

    ktna_df = exo_df[['pl_hostname','MatchID']].merge(
        _ktna, how='right', left_on='pl_hostname', right_on='name')
    print(ktna_df)




if __name__=="__main__":
    sector_number = 2 # 1-based !!

    print('\nWRN!: be sure you passed 1-based sector_number')
    main(sector_number)
