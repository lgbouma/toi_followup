from glob import glob
from lightcurve_utils import make_LS_ACF_periodograms

def check_62483237():

    lcdir = '../data/alert_lightcurves/'
    fname = glob(lcdir+'*62483237*')
    assert len(fname) == 1
    fname = fname[0]

    make_LS_ACF_periodograms(fname)

if __name__ == "__main__":
    check_62483237()
