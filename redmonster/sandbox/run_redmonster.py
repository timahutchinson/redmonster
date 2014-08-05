# Initial top-level script to run redmonster on a single plate, with options to use only a
# subset of fibers.  Each spPlate file is assumed to be in $BOSS_SPECTRO_REDUX/$RUN2D/pppp/
# where pppp is the 4-digit plate-id.
#
# Tim Hutchinson, University of Utah, July 2014
# t.hutchinson@utah.edu

import numpy as n
from astropy.io import fits
from redmonster.datamgr import spec
from redmonster.physics import zfinder, zfitter, zpicker
from redmonster.sandbox import yanny as y
from time import gmtime, strftime
import matplotlib.pyplot as p
p.interactive(True)

''' Set plate, mjd, and fibers to be run.  If fiberid is not specified here and subsequently passed in during the next step, the default behavior is to run on all fibers. '''
plate = 3686
mjd = 55268
fiberid = [i+100 for i in xrange(1)] # fiberid must be a list, not a numpy array


''' Read spPlate file.  specs.flux, specs.ivar, specs.loglambda, are [nfibers, npix] arrays containing flux, inverse variances, and log-wavelength, respectively.  This step also flags sky fibers and masks pixels with unreasonable S/N. '''
specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)

''' Instantiate zfinder object that will do z-finding for the entire plate using a single template.  Here, fname is the template filename, assumed to be in $REDMONSTER_DIR/templates/ . npoly specifies number of polynomial terms to be used in finding, zmin and zmax are upper and lower bounds of the redshift range to be explored. '''
zssp = zfinder.Zfinder(fname='ndArch-ssp_em_galaxy-v000.fits', type='GALAXY', npoly=4, zmin=-0.01, zmax=1.2)

''' Run actual fitting routine on the object created above. zssp.zchi2arr is array of of minimum chi^2 values of shape [nfibers, ndim_1, ndim_2, ..., ndim_N, nredshifts], where ndim_i is the i'th dimension of the template file.  Input flux need not be SDSS data - any spectra binned in constant log(lambda) will work.'''
zssp.zchi2(specs.flux, specs.loglambda, specs.ivar)

''' New object and fitting for a different template. '''
zstar = zfinder.Zfinder(fname='ndArch-spEigenStar-55734.fits', type='STAR', npoly=4, zmin=-.005, zmax=.005)
zstar.zchi2(specs.flux, specs.loglambda, specs.ivar)

''' Instantiate Zfitter to do subgrid fitting.  zchi2_ssp is chi^2 array from zfinder object above, and zssp.zbase is redshift-pixel baseline over the range explored by zfinder. '''
zfit_ssp = zfitter.Zfitter(zssp.zchi2arr, zssp.zbase)

''' Do actual subgrid refinement and fitting.  Best-fit and second-best-fit redshifts will be in [nfibers,2] shaped array in zfit_*****.z , and associated errors are in zfit_*****.z_err .  This routine also flags for small delta chi^2 and if the global minimum chi^2 is on an edge of the input chi^2 array. z_refine() includes an optional parameter, 'width', which specifies half the width of a window around the global chi2 minimum in which local minima are ignored during calculation of second-best redshift and flagging small delta-chi2.  If not specified, 'width' defaults to 15.'''
zfit_ssp.z_refine()

''' Same as above for second template. '''
zfit_star = zfitter.Zfitter(zstar.zchi2arr, zstar.zbase)
zfit_star.z_refine()

''' Compare chi2 surfaces from each template and classify each object accordingly. Arguments are number of pixels in DATA spectrum, followed by each object created by Zfinder.  This function can currently handle up to five objects from five separate templates.'''
zpick = zpicker.Zpicker(specs, zssp, zfit_ssp, zstar, zfit_star)

''' Flagging throughout redmonster is done individually by the classes responsible for handling the relevant computations.  To have an 'overall' flag for each fiber, the individual flags need to be combined. '''
ssp_flags = n.zeros(len(fiberid))
star_flags = n.zeros(len(fiberid))
for ifiber in xrange(len(fiberid)):
    ssp_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zssp.zwarning[ifiber])) | int(zfit_ssp.zwarning[ifiber])
    star_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zstar.zwarning[ifiber])) | int(zfit_star.zwarning[ifiber])

''' Write output file.  Arguments are zpick from above, and optionally dest and clobber, the path in which to write to file and whether or not to clobber old files with the same name, respectively.  See class documentation for more detail on Write_Redmonster behavior.'''
output = io.Write_Redmonster(zpick)