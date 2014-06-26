#
# multi_exposure_dev.scr.py
#
# script for developping methods that handle multiple
# exposures of a single fiber.
#
import numpy as n
import matplotlib as m
m.interactive(True)
from matplotlib import pyplot as p
from astropy.io import fits
import os
from redmonster.math import misc
#from scipy import signal as sig
from redmonster.datamgr import sdss
from redmonster.datamgr import io
from scipy import optimize as opt
import gc
import multifit as mf

#import pixelsplines as pxs

# Set the following:
# export BOSS_SPECTRO_REDUX='/data/BOSS/redux/dr10mini'
# export RUN2D=v5_5_12
# export RUN1D=v5_5_12

plate = 3686
mjd = 55268
fiberid = 265

# Get the data:
SpC = sdss.SpCFrameAll(plate, mjd)
SpC.set_fiber(fiberid)

spZbest_file = os.getenv('BOSS_SPECTRO_REDUX') + '/' + os.getenv('RUN2D') + '/' + \
               str(plate).rjust(4,'0') + '/' + os.getenv('RUN1D') + '/spZbest-' + \
               str(plate).rjust(4,'0') + '-' + str(mjd) + '.fits'

spz = fits.getdata(spZbest_file, 1)

# Get the models:
data, baselines, infodict = io.read_ndArch('../templates/ndArch-ssp_hires_galaxy-v001.fits')
loglam_ssp = infodict['coeff0'] + infodict['coeff1'] * n.arange(infodict['nwave'])
logbound_ssp = misc.cen2bound(loglam_ssp)
wave_ssp = 10.**loglam_ssp
wavebound_ssp = 10.**logbound_ssp
n_age = data.shape[0]

# What we probably want to do is build the broadening matrix from
# the higher resolution at more or less the same rest-frame coverage,
# then slide that within the template set.
# We should in the end also just do our velocity broadening
# on the models and have a precomputed grid on hand for the data.

# What is the inherent velocity width of the models as they
# have currently been imported?
c_kms = 2.99792458e5
vdisp_ssp = c_kms * n.log(10.) * infodict['coeff1']
# A little over 17 km/s in this case.
# So, let's specify what vdisp baseline we want...
# I think we should increment logarithmically.
#dlog10_vdisp = 0.075
#n_vdisp = 25
#vdisp_base = vdisp_ssp * 10.**(n.arange(n_vdisp) * dlog10_vdisp)#
#
# Or maybe not...
# That seems to stack up oddly in terms of what we actually
# want to sample as the more common values in the universe.
# Let's go linearly from 25 km/s up to something very large...

dvdisp = 25.
n_vdisp = 32
vdisp_base = dvdisp * (1. + n.arange(n_vdisp))

# OK, now we need to sweat the issue of a velocity-convolution
# buffer in constructing our velocity-broadened grid of SSPs.
# If we "say" we're going to 1000 km/s, then 6 sigma is 6000 km/s.
# So how many pixels is that in the current input sampling?
pixbuff_ssp = int(n.ceil(6000. / vdisp_ssp))

# So here is our "target" baseline for the SSPs:
loglam_vmodel = loglam_ssp[pixbuff_ssp:-pixbuff_ssp].copy()
logbound_vmodel = misc.cen2bound(loglam_vmodel)
wave_vmodel = 10.**loglam_vmodel
wavebound_vmodel = 10.**logbound_vmodel

# How much broadening do we need to impart to the SSPs in
# current form, in order to get the desired vdisp?
vdisp_net = n.sqrt(vdisp_base**2 - vdisp_ssp**2)

# Here's the bit where we do the velocity broadening
# to make the new grid that includes multiple vdisps.
# This kind of takes a while, so maybe this should
# be precomputed and stuffed into a file.
big_grid = n.zeros((n_vdisp, n_age, len(loglam_vmodel)), dtype=float)

for i_vdisp in xrange(n_vdisp):
    print i_vdisp
    sigma = wave_ssp * vdisp_net[i_vdisp] / c_kms
    v_matrix = misc.gaussproj(wavebound_ssp, sigma, wavebound_vmodel)
    for j_age in xrange(n_age):
        big_grid[i_vdisp, j_age] = v_matrix * data[j_age]
    gc.collect()

# Now we need to sort out the sigma broadening arrays for the individual exposures.

i=-1
i+=1
p.plot(10.**SpC.loglam_fib[i], SpC.disp_fib[i], hold=False)
# OK, looks sensible, although we've got the blue side padded,
# which probably wasn't necessary, but maybe made the downstream
# IDL work better with the spCFrames...

# First, we need to sort out observed-frame wavelength buffers
# around each of the individual exposures, because that will
# define the input space that we index within the models.

# I will assume at most a 3-pixel sigma, which means an 18-pixel
# buffer on either end.  Some computations:

logbound_fib = [misc.cen2bound(this_loglam) for this_loglam in SpC.loglam_fib]
wave_fib = [10.**this_loglam for this_loglam in SpC.loglam_fib]
wavebound_fib = [10.**this_logbound for this_logbound in logbound_fib]

minwave_fib = [min(this_wave) for this_wave in wave_fib]
maxwave_fib = [max(this_wave) for this_wave in wave_fib]

# Convert LSF SDSS-pixel sigmas into wavelength sigmas:
sigwave_fib = [0.0001 * n.log(10.) * wave_fib[ispec] * SpC.disp_fib[ispec]
               for ispec in xrange(SpC.nspec_fib)]

# To be safe, let's go ten sigmas beyond in either direction:
wavelim_lo = [wave_fib[k][0] - 10. * sigwave_fib[k][0] for k in xrange(SpC.nspec_fib)]
wavelim_hi = [wave_fib[k][-1] + 10. * sigwave_fib[k][-1] for k in xrange(SpC.nspec_fib)]

# Now, what are the indices of the bluest and reddest
# wavelength values of each exposure within the model grid?
idx_lo = [n.argmin(n.abs(wave_vmodel - this_lim)) for this_lim in wavelim_lo]
idx_hi = [n.argmin(n.abs(wave_vmodel - this_lim)) for this_lim in wavelim_hi]

# Next we need to interpolate the instrumental sigmas onto the input space:
sigwave_input = [n.interp(wave_vmodel[idx_lo[k]:idx_hi[k]+1], wave_fib[k],
                          sigwave_fib[k]) for k in xrange(SpC.nspec_fib)]

#k = -1
#k+=1
#p.plot(wave_fib[k], sigwave_fib[k], hold=False)
#p.plot(wave_vmodel[idx_lo[k]:idx_hi[k]+1], sigwave_input[k], hold=True)
# Looks right...

# Now we can (finally!) build the instrumental projection matrices:

inst_proj_fib = [misc.gaussproj(wavebound_vmodel[idx_lo[k]:idx_hi[k]+2],
                                sigwave_input[k], wavebound_fib[k])
                 for k in xrange(SpC.nspec_fib)]

# (See if our packaged function returns the same thing:)
matrix_list, idx_list, nsamp_list = mf.multi_projector(wavebound_fib, sigwave_fib, coeff0, coeff1)
k = 1
p.plot(wave_fib[k], matrix_list[k] * big_grid[15,12,idx_list[k]:idx_list[k]+nsamp_list[k]], hold=False)
p.plot(wave_fib[k], inst_proj_fib[k] * big_grid[15,12,idx_lo[k]:idx_hi[k]+1], hold=True)
# Yes, seems to be correct!


# Eventually we want to loop over redshift-lags and
# velocity-dispersions, but for testing right now, we will
# just dial in the "known" values so that we can get some sort
# of fit up and running...

z_best = 0.63034
v_best = 172.
idx_v = n.argmin(n.abs(vdisp_base - v_best))

# Pixel lag within the models to give this redshift:
pixlag = int(round(n.log10(1. + z_best) / infodict['coeff1']))

# Shall we stuff the projections into a list?
# Yes, probably...
proj_grid = [n.zeros((n_age,len(this_flux)), dtype=float) for this_flux in SpC.flux_fib]

for i_exp in xrange(SpC.nspec_fib):
    #print i_exp
    for j_age in xrange(n_age):
        #print j_age
        proj_grid[i_exp][j_age] = inst_proj_fib[i_exp] * \
                                  big_grid[idx_v,j_age,idx_lo[i_exp]-pixlag:idx_hi[i_exp]+1-pixlag]


#hold_val = [True] * SpC.nspec_fib
#hold_val[0] = False
#j_age = 12
#for k in xrange(SpC.nspec_fib):
#    p.plot(wave_fib[k], proj_grid[k][j_age], hold=hold_val[k])
# Looks good!!

# For the polynomial terms, let's try quadratic for now:
# npoly_fib = [2] * SpC.nspec_fib
npoly = 3

# This will build the non-negative polynomial component grids for
# each of the exposures.  For now, I *think* we want the same polynomial
# amplitude for each of the exposures...
maxloglam = max([max(this_loglam) for this_loglam in SpC.loglam_fib])
minloglam = min([min(this_loglam) for this_loglam in SpC.loglam_fib])
normbase_fib = [(this_loglam - minloglam) / (maxloglam - minloglam)
                for this_loglam in SpC.loglam_fib]

npix_fib = [len(this_flux) for this_flux in SpC.flux_fib]
poly_grid = [n.zeros((2*npoly, npix_this), dtype=float) for npix_this in npix_fib]

for ipoly in xrange(npoly):
    for jfib in xrange(SpC.nspec_fib):
        poly_grid[jfib][2*ipoly] = normbase_fib[jfib]**ipoly
        poly_grid[jfib][2*ipoly+1] = - normbase_fib[jfib]**ipoly

# Now we prep everything for the amplitude fitting:
big_a = n.hstack(proj_grid)
big_data = n.hstack(SpC.flux_fib)
big_ivar = n.hstack(SpC.invvar_fib)
big_poly = n.hstack(poly_grid)
big_ap = n.vstack((big_a, big_poly))
# Following just for plotting reference:
big_wave = n.hstack(wave_fib)

# Scaling as necessary for scipy.optimize.nnls:
big_dscale = big_data * n.sqrt(big_ivar)
big_ascale = big_ap * n.sqrt(big_ivar).reshape((1,-1))

coeffs, rnorm = opt.nnls(big_ascale.T, big_dscale)
big_model = n.dot(big_ap.T, coeffs)

p.plot(big_wave, big_data, '.', hold=False)
p.plot(big_wave, big_model, '.', hold=True)

chisq = n.sum((big_data-big_model)**2 * big_ivar)
# So, "rnorm" from nnls is the square root of chi-squared...

# See if our velocity broadened grids match up
# to those from the expernal precomputation:
junk, bjunk, ijunk = io.read_ndArch('../templates/ndArch-ssp_hires_galaxy-v002.fits')
wjunk = 10.**(ijunk['coeff0'] + ijunk['coeff1'] * n.arange(ijunk['nwave']))
j_v = 25
i_a = 8
p.plot(wjunk, junk[j_v,i_a], hold=False)
p.plot(wave_vmodel, big_grid[j_v,i_a], hold=True)

# Yes, they are the same...








# We are close, but not there yet.
# We still need to:
#   1. Supplement with polynomials as appropriate
#   2. Flatten it all into the overall linear algebraic form
#   3. Get the input math right for scipy.optimize.nnls
#   4. Verify that scipy.optimize.nnls works right
#   5. Plot the resulting fit along with the data
# then assuming that all works,
#   6. Embed this in a loop over redshifts and vdisps
#   7. Encapsulate it all in a function or class
#   8. Test on other fibers
#   9. Apply also to coadds and see how it looks
# then assuming THAT works...
#   10. Add emission line components
#   11. Add star handling
#   12. Consider soaking up the fluxing nuisance vectors
#   13. Etc. etc.






####
#### Below here is obsolete code...
####




maxwave_fib = [max(wave_fib[iexp][SpC.disp_fib[i] > 0]) for iexp in xrange(SpC.nspec_fib)]

argminwave_sig = [min(n.where(this_disp > 0)[0]) for this_disp in SpC.disp_fib]
argmaxwave_sig = [max(n.where(this_disp > 0)[0]) for this_disp in SpC.disp_fib]
minwave_sig = [wave_fib[ispec][argminwave_sig[ispec]] for ispec in xrange(SpC.nspec_fib)]
maxwave_sig = [wave_fib[ispec][argmaxwave_sig[ispec]] for ispec in xrange(SpC.nspec_fib)]

argminwave_ivar = [min(n.where(this_ivar > 0)[0]) for this_ivar in SpC.invvar_fib]
argmaxwave_ivar = [max(n.where(this_ivar > 0)[0]) for this_ivar in SpC.invvar_fib]



argminwave_sig = [n.arange(len((SpC.disp_fib[i] > 0) for iexp in xrange(SpC.nspec_fib)]
argmaxwave_sig = [max(wave_fib[iexp][SpC.disp_fib[i] > 0]) for iexp in xrange(SpC.nspec_fib)]





# Code to get the individual exposures for this plate:
spf = os.getenv('BOSS_SPECTRO_REDUX') + '/' + os.getenv('RUN2D') + '/' + \
      str(plate).strip() + '/spPlate-' + str(plate).strip() + '-' + \
      str(mjd).strip() + '.fits'
hdr = fits.getheader(spf)

exp_keys = [this_key for this_key in hdr.keys() if this_key[:5] == 'EXPID']

exp_ids = [hdr[this_key][:11] for this_key in exp_keys]

# Build the individual file exposure names:
path_to_spectra = os.getenv('BOSS_SPECTRO_REDUX') + '/' + \
                  os.getenv('RUN2D') + '/' + str(plate).strip() + '/'
spCFrame_list = [path_to_spectra + 'spCFrame-' + this_id + '.fits' for this_id in exp_ids]

# Get the data:
data_list = [fits.getdata(this_file) for this_file in spCFrame_list]
invvar_list = [fits.getdata(this_file, 1) for this_file in spCFrame_list]
loglam_list = [fits.getdata(this_file, 3) for this_file in spCFrame_list]
sigma_list = [fits.getdata(this_file, 4) for this_file in spCFrame_list]
plug_list = [fits.getdata(this_file, 5) for this_file in spCFrame_list]

# Find the indices of the fiberid of interest:
i_exp = []
j_row = []
for i in xrange(len(plug_list)):
    wh_fib = n.where(plug_list[i].FIBERID == fiberid)[0]
    for j in xrange(len(wh_fib)):
        i_exp.append(i)
        j_row.append(wh_fib[j])

n_spec = len(i_exp)

# Pull out the bits we actually want:
data_list_one = [data_list[i_exp[k]][j_row[k]] for k in xrange(n_spec)]
invvar_list_one = [invvar_list[i_exp[k]][j_row[k]] for k in xrange(n_spec)]
loglam_list_one = [loglam_list[i_exp[k]][j_row[k]] for k in xrange(n_spec)]
sigma_list_one = [sigma_list[i_exp[k]][j_row[k]] for k in xrange(n_spec)]

# Derived vectors:
wave_list_one = [10.**loglam_this for loglam_this in loglam_list_one]
logbound_list_one = [misc.cen2bound(loglam_this) for loglam_this in loglam_list_one]
wavebound_list_one = [10.**logbound_this for logbound_this in logbound_list_one]
sigwave_list_one = [10.**(-4) * n.log(10.) * wave_list_one[i] * sigma_list_one[i]
                    for i in xrange(n_spec)]


# That basically does it for the data.
# Next thing is to get the models into a compatible form!

# Although perhaps we should encapsulate this spCFrame
# handling code into a slicker interface...

from redmonster.datamgr import sdss

SpC = sdss.SpCFrameAll(plate, mjd)

SpC.set_fiber(fiberid)

# Plot this and see if it makes sense:
hold_val = n_spec * [True]
hold_val[0] = False

for i in xrange(n_spec):
    p.plot(wave_list_one[i], data_list_one[i], hold=hold_val[i], color='k')

for i in xrange(SpC.nspec_fib):
    p.plot(10.**SpC.loglam_fib[i], SpC.flux_fib[i], hold=hold_val[i], color='k')



for i in xrange(n_spec):
    p.plot(wave_list_one[i], sigwave_list_one[i], hold=hold_val[i], color='k')





wh_fib = [n.where(this_plug.FIBERID == fiberid)[0] for this_plug in plug_list]

[
