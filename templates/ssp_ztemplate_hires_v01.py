import numpy as n
import matplotlib as m
m.interactive(True)
from matplotlib import pyplot as p
from astropy.io import fits
#from astropy import wcs
import os
#import pixelsplines as pxs
from redmonster.datamgr import io
from redmonster.physics import pixelsplines as pxs
import copy
from redmonster.physics import misc

# ssp_ztemplate_hires_v01.py
# Script to generate redshift templates at higher-than-SDSS sampling
# from Charlie Conroy's theoretical SSPs.
# bolton@utah@iac 2014junio


# Set this:
# export SSP_DATA_DIR=/data/SSP

# First read in the templates:
sspf = os.getenv('SSP_DATA_DIR') + '/SSP_Padova_RRLIB_Kroupa_Z0.0190.fits'
ssp_data = fits.getdata(sspf, 1)
ssp_meta = fits.getdata(sspf, 2)
ssp_wave = ssp_data['LAMBDA'][0].copy()
ssp_flux = ssp_data['SPEC'][0].copy()
ssp_agegyr = ssp_meta['AGEGYR'].copy()
n_age, n_pix = ssp_flux.shape

# Need to sort out conversion from Lnu to Llambda..
# These come in Lsun per Hz.
# So, we need to first multiply by...
Lsun_in_cgs = 3.839e33
# which will get us to ergs per second per Hz.
# Then we need to convert from angstroms to cm
# for purposes of flux-unit conversion:
ang_per_cm = 1.e8
#ssp_wave_cm = ssp_wave / ang_per_cm
# We also need speed of light in cgs units:
c_cgs = 2.99792458e10
# So, the conversion vector is:
#(Lsun_in_cgs * c_cgs / (ssp_wave_cm**2)) * ang_per_cm
# where the last ang_per_cm puts it in per-angstrom.
# This is all equivalent to:
fluxconv = Lsun_in_cgs * c_cgs / (ssp_wave**2 * ang_per_cm)
# by cancelling one of the powers of ang_per_cm,
# since (again) ssp_wave is in Angstroms.

# the actual conversion:
ssp_flam = ssp_flux * fluxconv.reshape((1,-1))

# Exponent to divide out:
exp_div = 12
ssp_flam /= 10.**exp_div

# Get a subset of the templates:
izero = n.argmin(n.abs(n.log10(ssp_agegyr) - (-2.5)))
idx = 10 * n.arange(15) + izero
n.log10(ssp_agegyr[idx])
nsub_age = len(idx)

# Work out the pixel boundaries and widths
# (we go to log-lambda because the
# spacing is uniform in that):
ssp_loglam = n.log10(ssp_wave)
ssp_logbound = pxs.cen2bound(ssp_loglam)
ssp_wavebound = 10.**ssp_logbound
#ssp_dwave = ssp_wavebound[1:] - ssp_wavebound[:-1]
# Note that there is some slight unsmoothness in the
# width of the sampling intervals.  I think this is fine,
# but it will imply a slightly variable resolution
# implicit in the data as it currently stands if we dial
# it in straight from the pixel widths.  So instead, let's
# try setting it based on the average loglam spacing,
# since that is essentially uniform.
# p.plot(ssp_logbound[1:] - ssp_logbound[:-1], hold=False)
dloglam = n.mean(ssp_logbound[1:] - ssp_logbound[:-1])

# This will work:
#ssp_dwave = 10.**(ssp_loglam + 0.5 * dloglam) - \
#            10.**(ssp_loglam - 0.5 * dloglam)
#ssp_dwave = ssp_wave * (10**(dloglam/2.) - 10.**(-dloglam/2.))

# THIS is the best way to do it!!
ssp_dwave = dloglam * ssp_wave * n.log(10.)

# What is the resolution that we want to add in order
# to make the SDSS pixel width?
rebin_dloglam = 0.000025
rebin_dwave = rebin_dloglam * ssp_wave * n.log(10.)

# This constructs the blurring matrix to output resolution:
rebin_blur = pxs.gauss_blur_matrix(ssp_wavebound, n.sqrt(rebin_dwave**2 - ssp_dwave**2))

# wavelength range for rebinned version to be bounded by:
wave_rebin_lo = 1525.
wave_rebin_hi = 10850.
coeff1 = rebin_dloglam
coeff0 = coeff1 * n.floor(n.log10(wave_rebin_lo) / coeff1)
naxis1 = int(n.ceil(1. + (n.log10(wave_rebin_hi) - coeff0) / coeff1))

loglam_rebin = coeff0 + coeff1 * n.arange(naxis1)
logbound_rebin = pxs.cen2bound(loglam_rebin)
wavebound_rebin = 10.**logbound_rebin

# Initialize array for rebinned blurred SSPs:
ssp_rebin = n.zeros((nsub_age, naxis1), dtype=float)

# Do the rebinning:
for i_age in xrange(nsub_age):
    print(i_age)
    pxspline = pxs.PixelSpline(ssp_wavebound, rebin_blur * ssp_flam[idx[i_age]])
    ssp_rebin[i_age] = pxspline.resample(wavebound_rebin)

baselines = [n.log10(ssp_agegyr[idx])]

this_class = 'ssp_hires_galaxy'
this_version = 'v001'

infodict = {'par_names': n.asarray(['log10-age']),
            'par_units': n.asarray(['log10-Gyr']),
            'par_axistype': n.asarray(['regular']),
            'coeff0': coeff0, 'coeff1': coeff1,
            'fluxunit': '10^'+str(exp_div)+' erg/s/Ang/M_sun_init',
            'filename': 'ndArch-'+this_class+'-'+this_version+'.fits'}

io.write_ndArch(n.float32(ssp_rebin), baselines, infodict)

#junk, bjunk, ijunk = io.read_ndArch(infodict['filename'])

