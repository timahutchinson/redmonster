# Initial top-level script to run redmonster on a single plate, with options to use only a
# subset of fibers.  Each spPlate file is assumed to be in $BOSS_SPECTRO_REDUX/$RUN2D/pppp/
# where pppp is the 4-digit plate-id.
#
# Tim Hutchinson, July 2014
# t.hutchinson@utah.edu

import numpy as n
from astropy.io import fits
from redmonster.datamgr import spec
from redmonster.physics import zfinder, zfitter
from redmonster.sandbox import yanny as y
from time import gmtime, strftime
import matplotlib.pyplot as p
p.interactive(True)

''' Set plate, mjd, and fibers to be run.  If fiberid is not specified here and subsequently passed in during the next step, the default behavior is to run on all fibers. '''
plate = 3686
mjd = 55268
fiberid = [i for i in xrange(1000)] # fiberid must be a list, not a numpy array


''' Read spPlate file.  specs.flux, specs.ivar, specs.loglambda, are [nfibers, npix] arrays containing flux, inverse variances, and log-wavelength, respectively.  This step also flags sky fibers and masks pixels with unreasonable S/N. '''
specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)

''' Instantiate zfinder object that will do z-finding for the entire plate using a single template.  Here, fname is the template filename, assumed to be in $REDMONSTER_DIR/templates/ . npoly specifies number of polynomial terms to be used in finding, zmin and zmax are upper and lower bounds of the redshift range to be explored. '''
zssp = zfinder.Zfinder(fname='ndArch-ssp_em_galaxy-v000.fits', npoly=4, zmin=-0.01, zmax=1.2)

''' Run actual fitting routine on the object created above. Output is array of of minimum chi^2 values of shape [nfibers, ndim_1, ndim_2, ..., ndim_N, nredshifts], where ndim_i is the i'th dimension of the template file.'''
zchi2_ssp = zssp.zchi2(specs.flux, specs.loglambda, specs.ivar)

''' New object and fitting for a different template. '''
zstar = zfinder.Zfinder(fname='ndArch-spEigenStar-55734.fits', npoly=4, zmin=-.005, zmax=.005)
zchi2_star = zstar.zchi2(specs.flux, specs.loglambda, specs.ivar)

''' Instantiate zfitter to do subgrid fitting.  zchi2_ssp is chi^2 array from zfinder object above, and zssp.zbase is redshift-pixel baseline over the range explored by zfinder. '''
zfit_ssp = zfitter.Zfitter(zchi2_ssp, zssp.zbase)

''' Do actual subgrid refinement and fitting.  "Best fit" redshifts will be in [nfibers] shaped array in zfit_*****.best_z , and associated errors are in zfit_*****.z_err .  This routine also flags for small delta chi^2 and if the global minimum chi^2 is on an edge of the input chi^2 array. '''
zfit_ssp.z_refine()

''' Same as above for second template. '''
zfit_star = zfitter.Zfitter(zchi2_star, zstar.zbase)
zfit_star.z_refine()

''' Flagging throughout redmonster is done individually by the classes responsible for handling the relevant computations.  To have an 'overall' flag for each fiber, the individual flags need to be combined. '''
ssp_flags = n.zeros(len(fiberid))
star_flags = n.zeros(len(fiberid))
for ifiber in xrange(len(fiberid)):
    ssp_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zssp.zwarning[ifiber])) | int(zfit_ssp.zwarning[ifiber])
    star_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zstar.zwarning[ifiber])) | int(zfit_star.zwarning[ifiber])