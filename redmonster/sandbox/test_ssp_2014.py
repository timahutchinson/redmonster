import numpy as n
import matplotlib as m
m.interactive(True)
from matplotlib import pyplot as p
from astropy.io import fits
import os
import pixelsplines as pxs

# Script to test new SSPs from Charlie Conroy, obtained 2014A
#
# Written by A. Bolton, U. of Utah

# Set these environment variables:
# export BOSS_SPECTRO_REDUX=/data/BOSS/redux/dr10mini
# export RUN2D=v5_5_12
# export RUN1D=v5_5_12
# export SSP_DATA_DIR=/data/SSP

# Get a spectrum to work with:
plate = 4388
mjd = 55536
spf = os.getenv('BOSS_SPECTRO_REDUX') + '/' + os.getenv('RUN2D') + \
      '/' + str(plate) + '/spPlate-' + str(plate) + '-' + str(mjd) + '.fits'

flux = fits.getdata(spf, 0)
invvar = fits.getdata(spf, 1)
hdr = fits.getheader(spf)
loglam = hdr['COEFF0'] + hdr['COEFF1'] * n.arange(hdr['NAXIS1'])
n_fib, n_data = flux.shape
ifiber = 100
p.plot(10**loglam, flux[ifiber], hold=False)


# Get the solar-metallicity SSP:
sspf = os.getenv('SSP_DATA_DIR') + '/SSP_Padova_CKC14_new_Kroupa_Z0.0190.out.fits'
ssp_data = fits.getdata(sspf, 1)
ssp_meta = fits.getdata(sspf, 2)
ssp_wave = ssp_data['LAMBDA'][0].copy()
ssp_flux = ssp_data['SPEC'][0].copy()
n_age, n_pix = ssp_flux.shape

# Convert from fnu to flambda with arbitrary normalization:
for i in xrange(n_age):
    ssp_flux[i] /= ssp_wave**2
    ssp_flux[i] /= n.median(ssp_flux[i])

# Check it out, to see if it looks like we expect:
i_age = 170
p.plot(ssp_wave, ssp_flux[i_age], hold=False)
p.xlim(1800.,11000.)


# Pick out some wavelength range of interest:
# Nah, not worth it...
#lam_lo = 1750.
#lam_hi = 11000.
#ilam_lo = min(n.where(ssp_wave >= lam_lo)[0])
#ilam_hi = max(n.where(ssp_wave <= lam_hi)[0])

# Does a blur matrix at this stage give us something tractable?
vdisp = 225. # (in km/s)
c_in_km_per_s = 299792.458
sigconv = ssp_wave * vdisp / c_in_km_per_s
ssp_bound = pxs.cen2bound(ssp_wave)
gblur = pxs.gauss_blur_matrix(ssp_bound, sigconv)

ssp_blur = gblur * ssp_flux[i_age]
p.plot(ssp_wave, ssp_flux[i_age], hold=False)
p.plot(ssp_wave, ssp_blur, hold=True)
p.xlim(1800.,11000.)

# Looks good.
# Now let's bin it down to the BOSS coadd resoultion
# and dial in the wavelength range we want.
sspSpline = pxs.PixelSpline(ssp_bound, ssp_blur)

ssp_coeff0 = 3.225
ssp_coeff1 = 0.0001
ssp_naxis1 = 2**13
ssp_loglam = ssp_coeff0 + ssp_coeff1 * n.arange(ssp_naxis1)
print 10.**ssp_loglam.min(), 10.**ssp_loglam.max()
ssp_logbound = pxs.cen2bound(ssp_loglam)
ssp_wavebound = 10.**ssp_logbound

ssp_boss = sspSpline.resample(ssp_wavebound)

p.plot(10.**ssp_loglam, ssp_boss, hold=False)

# Test speed of FFTs from Numpy:
junk1 = n.random.uniform(size=2**17)
junk2 = n.random.uniform(size=(2**17-1))
f_junk1 = n.fft.fft(junk1)
f_junk2 = n.fft.fft(junk2)
# It definitely matters to have power-of-two!!

# Another test of the FFT convolution convention:
junk1 = n.zeros(128, dtype=float)
junk2 = n.zeros(128, dtype=float)
junk1[0] = 1.
junk2[0] = 2.
f_junk1 = n.fft.fft(junk1)
f_junk2 = n.fft.fft(junk2)
junk3 = n.fft.ifft(f_junk1 * f_junk2)

# Yet another test of the FFT convolution convention:
junk1 = n.random.uniform(size=128)
junk2 = n.random.uniform(size=128)
f_junk1 = n.fft.fft(junk1)
f_junk2 = n.fft.fft(junk2)
junk3 = n.fft.ifft(f_junk1.conj() * f_junk2).real
junk4 = n.convolve(junk1, junk2)

p.plot(junk3, hold=False)
p.plot(junk4, hold=True)


# And another test of the FFT convolution convention:
junk1 = n.random.uniform(size=128)
junk2 = 0.1*n.random.uniform(size=128)
#f_junk1 = n.fft.fft(junk1)
#f_junk2 = n.fft.fft(junk2)
#junk3 = n.fft.ifft(f_junk1 * f_junk2).real
f_junk1 = n.fft.fft(junk1)
f_junk2 = n.fft.fft(junk2)
junk3 = n.fft.ifft(f_junk1 * f_junk2.conj()).real
junk4 = n.zeros(128,dtype=float)
for i in xrange(128):
    junk4[i] = n.sum(junk1 * n.roll(junk2,i))

p.plot(junk3, hold=False)
p.plot(junk4, hold=True)

# OK, that does what we want, I think!!



# Make a padded version of the data and invvar:
data_pad = n.zeros(ssp_naxis1, dtype=float)
data_pad[0:len(flux[ifiber])] = flux[ifiber]
invvar_pad = n.zeros(ssp_naxis1, dtype=float)
invvar_pad[0:len(invvar[ifiber])] = invvar[ifiber]

# Generate simple polynomial basis:
npoly = 3
poly_base = n.arange(ssp_naxis1, dtype=float) / float(ssp_naxis1-1)
poly_set = n.zeros((npoly,ssp_naxis1), dtype=float)
for i in xrange(npoly):
    poly_set[i] = poly_base**i

# The FFTs we need:
ivar_fft = n.fft.fft(invvar_pad)
t_fft = n.fft.fft(ssp_boss)
t2_fft = n.fft.fft(ssp_boss**2)
poly_fft = n.zeros((npoly, ssp_naxis1), dtype=complex)
for i in xrange(npoly):
    poly_fft[i] = n.fft.fft(poly_set[i] * invvar_pad)

data_fft = n.fft.fft(data_pad * invvar_pad)

# Build and populate the array of inversion matrices and right-hand sides:
alpha_big = n.zeros((npoly+1,npoly+1,ssp_naxis1), dtype=float)
rhs_big = n.zeros((npoly+1,ssp_naxis1), dtype=float)
alpha_big[0,0] = n.fft.ifft(t2_fft * ivar_fft.conj()).real
rhs_big[0] = n.fft.ifft(t_fft * data_fft.conj()).real
for i in xrange(npoly):
    alpha_big[i+1,0] = alpha_big[0,i+1] = n.fft.ifft(t_fft * poly_fft[i].conj()).real

ipoly_set = poly_set * invvar_pad.reshape((1,-1))
alpha_poly = n.tensordot(poly_set, ipoly_set, (1,1))
rhs_poly = n.sum(poly_set * data_pad.reshape((1,-1)) * invvar_pad.reshape((1,-1)), axis=1)

alpha_big[1:,1:] = alpha_poly.reshape((npoly,npoly,1)) * n.ones((1,1,ssp_naxis1))
rhs_big[1:] = rhs_poly.reshape((npoly,1)) * n.ones((1,ssp_naxis1))

# Test these to see if they give the same numbers as one expects
# from straight calculation, for the zero-lag case:
print n.sum(ssp_boss**2 * invvar_pad), alpha_big[0,0,0]
print n.sum(ssp_boss * poly_set[0] * invvar_pad), alpha_big[0,1,0], alpha_big[1,0,0]
print n.sum(ssp_boss * poly_set[1] * invvar_pad), alpha_big[0,2,0], alpha_big[2,0,0]
print n.sum(ssp_boss * poly_set[2] * invvar_pad), alpha_big[0,3,0], alpha_big[3,0,0]
print n.sum(ssp_boss * data_pad * invvar_pad), rhs_big[0,0]
print n.sum(poly_set[1] * poly_set[2] * invvar_pad), alpha_big[2,3,0], alpha_big[3,2,0]
print n.sum(poly_set[1] * data_pad * invvar_pad), rhs_big[2,0]

# Number of redshifts to consider:
n_z = ssp_naxis1 - n_data + 1

# Here's where we find out if we're making sense...
sn_squared = n.zeros(n_z, dtype=float)

for i in xrange(n_z):
    coeffs = n.linalg.solve(alpha_big[:,:,i], rhs_big[:,i])
    sn_squared[i] = n.dot(coeffs, n.dot(alpha_big[:,:,i], coeffs))

chisq = n.sum(data_pad**2 * invvar_pad) - sn_squared

bestlag = n.argmin(chisq)

zbase = 10.**loglam[0] / 10.**ssp_loglam[0:n_z] - 1

zbase[bestlag]

# Revisit at the best lag and look explicitly at model:
best_coeffs = n.linalg.solve(alpha_big[:,:,bestlag], rhs_big[:,bestlag])
best_basis = n.zeros((npoly+1,n_data), dtype=float)
best_basis[0] = ssp_boss[bestlag:bestlag+n_data]
best_basis[1:] = poly_set[:,0:n_data]
best_model = n.dot(best_coeffs, best_basis)
p.plot(10.**loglam, flux[ifiber], hold=False)
p.plot(10.**loglam, best_model, hold=True)

print n.sum((flux[ifiber]-best_model)**2 * invvar[ifiber]), chisq[bestlag]

