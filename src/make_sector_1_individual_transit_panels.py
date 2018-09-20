from glob import glob

from lightcurve_utils import make_individual_transit_inspection_panel

def make_sector_1_individual_transit_panels():

    lcdir = '/Users/luke/local/tess_alert_lightcurves/'
    pattern = 'tess2018*-s0001-*-111-s_llc.fits.gz'
    fnames = glob(lcdir+pattern)

    # WASP-62, HATS-3, WASP-73, WASP-95, WASP-100, WASP-94A (planet host), HATS-13, ?,
    # HATS-30, ?, WASP-119, WASP-124, WASP-126, ?, WASP-91
    desired = [149603524, 336732616, 231670397, 144065872, 38846515, 92352620,
               289793076, 29344935, 281459670, 355703913, 388104525, 97409519,
               25155310, 281541555, 238176110]
    okids = [str(d) for d in desired]

    thesenames = []
    for f in fnames:
        for okid in okids:
            if okid in f:
                thesenames.append(f)

    for fname in thesenames:
        make_individual_transit_inspection_panel(fname)

if __name__ == "__main__":
    make_sector_1_individual_transit_panels()
