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
from redmonster.physics import airtovac as a2v

# cap_ztemplate_v01.py
# Script to generate redshift templates at SDSS sampling
# from Carlos Allende-Prieto's theoretical stellar spectra.
# bolton@utah@iac 2014mayo
# updated 2014junio for TiO inclusion


# Set this:
# export CAP_DATA_DIR=/data/CAP/coarse

# Set the filename:
# Need to loop over 1--5 somehow...
fname = os.getenv('CAP_DATA_DIR') + '/ndArch-nsc5-coarse00.fits'
# And then these:
#fname = os.getenv('CAP_DATA_DIR') + '/ndArch-DATNLTE-coarse00.fits'
#fname = os.getenv('CAP_DATA_DIR') + '/ndArch-DATLTE-coarse00.fits'

# Get the data:
data, baselines, infodict = io.read_ndArch(fname)

# Build the log-lambda baseline and pixel boundaries:
hires_loglam_air = infodict['coeff0'] + infodict['coeff1'] * n.arange(infodict['nwave'])
hires_logbound_air = pxs.cen2bound(hires_loglam_air)
hires_wave_air = 10.**hires_loglam_air
hires_wavebound_air = 10.**hires_logbound_air
hires_dwave_air = hires_wavebound_air[1:] - hires_wavebound_air[:-1]

# Transform the baselines to vacuum wavelengths:
hires_wave = a2v.a2v(hires_wave_air)
hires_wavebound = a2v.a2v(hires_wavebound_air)
hires_dwave = hires_wavebound[1:] - hires_wavebound[:-1]

# Not sure what the interpretation is at wavelengths below 2000ang...
# We will probably cut that out for now anyway, since these are the
# first-pass redshift templates that we are building.

# Now, we want a flattened view (not copy) of the hi-res data:
npars = data.size // infodict['nwave']
data_flat = data.reshape((npars,infodict['nwave']))

# This next bit converts to flambda in vacuum, via the flattened copy:
for i in xrange(npars):
    data_flat[i] *= (hires_dwave_air / hires_dwave)

# What is the resolution that we want to add in order
# to make the SDSS pixel width?
sdss_dloglam = 0.0001
sdss_dwave = sdss_dloglam * hires_wave * n.log(10.)

# This constructs the blurring matrix to SDSS resolution.
# We use hires_dwave_air to avoid the 2000Ang airtovac discontinuity
# without being significantly incorrect about wavelength.
sdss_blur = pxs.gauss_blur_matrix(hires_wavebound, n.sqrt(sdss_dwave**2 - hires_dwave_air**2))

# Wavelength range for rebinned SDSS version to be bounded by
# the following, which at least covers ugriz:
wave_sdss_lo = 2900.
wave_sdss_hi = 12000.
coeff1 = sdss_dloglam
coeff0 = coeff1 * n.floor(n.log10(wave_sdss_lo) / coeff1)
naxis1 = int(n.ceil(1. + (n.log10(wave_sdss_hi) - coeff0) / coeff1))

loglam_sdss = coeff0 + coeff1 * n.arange(naxis1)
logbound_sdss = pxs.cen2bound(loglam_sdss)
wavebound_sdss = 10.**logbound_sdss

# Initialize array for rebinned blurred stars:
cap_sdssbin = n.zeros((npars, naxis1), dtype=float)

# Do the blurring and rebinning:
for i in xrange(npars):
    print i
    blurspec = sdss_blur * data_flat[i]
    pxspline = pxs.PixelSpline(hires_wavebound, sdss_blur * data_flat[i])
    cap_sdssbin[i] = pxspline.resample(wavebound_sdss)


#i = -1L
#
#i += 1
#p.plot(10.**loglam_sdss, cap_sdssbin[i], drawstyle='steps-mid', hold=False)

# reshape the array:
out_shape = data[...,0].shape + (naxis1,)
cap_sdssbin.resize(out_shape)

out_infodict = copy.deepcopy(infodict)
out_infodict['version'] = infodict['version'] + '_lr'
out_infodict['filename'] = 'ndArch-' + out_infodict['class'] + '-' + out_infodict['version'] + '.fits'
out_infodict['coeff0'] = coeff0
out_infodict['coeff1'] = coeff1

io.write_ndArch(n.float32(cap_sdssbin), baselines, out_infodict)

#junk, bjunk, ijunk = io.read_ndArch(infodict['filename'])

