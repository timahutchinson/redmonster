#
# pixelsplines_demo.py
#
# Script to demonstrate Bolton's "pixelsplines" code.
#
# Meant to be executed incrementally/interactively
# from a Python command line.
#
# bolton@utah2013may
#


# Import statements.
import os
import numpy as n
import matplotlib as m
m.rc('text', usetex=True)
m.use('TkAgg')
m.interactive(True)
m.rc('font', size=16)
from matplotlib import pyplot as p
#import pyfits as pf
from astropy.io import fits as pf
import pixelsplines as pxs


# We'll assume that you have an up-to-date idlspec2d
# from which we'll get the toy data to play with:
tfile = os.getenv('IDLSPEC2D_DIR') + '/templates/spEigenStar-55734.fits'
hdulist = pf.open(tfile)
flux = hdulist[0].data
npix = hdulist[0].header['NAXIS1']
nstar = hdulist[0].header['NAXIS2']
coeff0 = hdulist[0].header['COEFF0']
coeff1 = hdulist[0].header['COEFF1']
loglam = coeff0 + coeff1 * n.arange(npix)
wave = 10.**loglam


#i = -1L
#i+=1
#p.plot(wave, flux[i], hold=False)

# Choose a particular spectrum that you like, and have a look at it:
istar = 40
p.plot(wave, flux[istar], hold=False, drawstyle='steps-mid')

# Now initialize a "pixelspline" object for this spectrum.
# If the units are f-lambda, then we really should use wavelength
# and not log-wavelength as the baseline to initialize.
# First convert the log-lambda pixel-center baseline
# to a log-lambda pixel-boundary baseline:
logbound = pxs.cen2bound(loglam)
wavebound = 10.**logbound
PS = pxs.PixelSpline(wavebound, flux[istar])


# Now let's try out a few of the methods...


# First, let's do a point-evaluation of the pixel spline
# over a finely gridded baseline:
nsub = 10
nfine = npix * nsub
# (The following line ensures that the first "nsub" points in logfine
# have a mean value equal to the first single point in loglam, and so on.)
logfine = coeff0 + coeff1 * (n.arange(nfine) - 0.5 * (nsub-1)) / float(nsub)
wavefine = 10.**logfine
fluxfine = PS.point_evaluate(wavefine)

# Plot point-evaluated and pixel-integrated together:
p.plot(wave, flux[istar], hold=False, drawstyle='steps-mid')
p.plot(wavefine, fluxfine, hold=True)

# Find the analytic minima of the underlying pixel spline,
# and overplot them:
wmin = PS.find_extrema(minima=True)
fmin = PS.point_evaluate(wmin)
p.plot(wmin, fmin, 'o', hold=True)

# Same, but for the maxima -- whatever that's worth in a stellar spectrum:
wmax = PS.find_extrema()
fmax = PS.point_evaluate(wmax)
p.plot(wmax, fmax, 'd', hold=True)

# Average up the finely gridded spectrum manually, and verify
# that it matches the original pixelized spectrum, within the
# expected limits of precision:
flux_reavg = fluxfine.reshape((npix,nsub)).sum(axis=1) / float(nsub)

# Plot both:
p.plot(wave, flux[istar], hold=False)
p.plot(wave, flux_reavg, hold=True)

# Plot fractional difference:
p.plot((flux[istar] - flux_reavg) / flux[istar], hold=False)

# Now do an analytic rebinning with, say, two new pixels
# for every three old pixels:
nper_new = 2
nper_old = 3
# Integer arithmetic to make the number of pixels work out right:
npix_old = nper_old * (npix / nper_old)
npix_new = nper_new * npix_old / nper_old
# Build new wavelength baseline and boundaries:
coeff1_new = nper_old * coeff1 / float(nper_new)
logbound_new = logbound[0] + coeff1_new * n.arange(npix_new+1)
loglam_new = logbound_new[:-1] + 0.5 * coeff1_new
wavebound_new = 10.**logbound_new
wave_new = 10.**loglam_new

flux_new = PS.resample(wavebound_new)

# Look at the original and the rebinned, overplotted against each other:
p.plot(wave[:npix_old], flux[istar][:npix_old], drawstyle='steps-mid', hold=False)
p.plot(wave_new, flux_new, drawstyle='steps-mid', hold=True)


# Now if we manually bin 3 against 2, we should get very exactly
# the same, to within numerical precision (not just sampling precision.)

check_old = flux[istar][0:npix_old].reshape((npix_old/nper_old,nper_old)).sum(axis=1) / float(nper_old)
check_new = flux_new.reshape((npix_new/nper_new,nper_new)).sum(axis=1) / float(nper_new)

p.plot((check_old - check_new) / check_old, hold=False)

# Looks pretty good.  There's probably some ringing off
# the end conditions that we could control to get the residuals
# down even further, but it's in pretty good shape.



# Last thing, out of the pixelsplines module, but not
# a method of the PixelSpline class, is the convolution
# by variable resolution.  This works analytically,
# but on pixelized data rather than on any spline
# form of data.

# Let's say we want to convolve our spectrum by something
# with a sigma that would be 1 Angstrom at 3000 Ang,
# and 1 Ang at 1.2 micron, and 2 Ang at 7500 Ang, varying
# quadratically with wavelength.

sig_conv = 1. - (wave - 3000.) * (wave - 12000.) / (4500.**2)
blurmatrix = pxs.gauss_blur_matrix(wavebound, sig_conv)

# Apply the blur matrix to our test spectrum to see what the outcome looks like:
blurflux = blurmatrix * flux[istar]
p.plot(wave, flux[istar], drawstyle='steps-mid', hold=False)
p.plot(wave, blurflux, drawstyle='steps-mid', hold=True)

# Let's make a second blur matrix, and a third that is purportedly
# the quadrature-sum of the two, and see whether the composition
# of the first two is equal to the third in its action:

sig_conv_2 = 3. - sig_conv
sig_conv_3 = n.sqrt(sig_conv**2 + sig_conv_2**2)

blurmatrix_2 = pxs.gauss_blur_matrix(wavebound, sig_conv_2)
blurmatrix_3 = pxs.gauss_blur_matrix(wavebound, sig_conv_3)

blurflux_12 = blurmatrix_2 * (blurmatrix * flux[istar])
blurflux_3 = blurmatrix_3 * flux[istar]

# Plot them together:
p.plot(wave, blurflux_12, drawstyle='steps-mid', hold=False)
p.plot(wave, blurflux_3, drawstyle='steps-mid', hold=True)

# What's the fractional difference?
p.plot(wave, (blurflux_12 - blurflux_3)/blurflux_3, drawstyle='steps-mid', hold=False)
# Pretty good, all things considered -- except at the edges.

