# Initial top-level script to run redmonster on a single plate, with options to use only a
# subset of fibers.  Each spPlate file is assumed to be in $BOSS_SPECTRO_REDUX/$RUN2D/pppp/
# where pppp is the 4-digit plate-id.  Templates are assumed to be in $REDMONSTER_DIR/templates/
#
# Tim Hutchinson, University of Utah, July 2014
# t.hutchinson@utah.edu

import numpy as n
from astropy.io import fits
from redmonster.datamgr import spec, io, io2
from redmonster.sandbox import donna_spec
from redmonster.physics import zfinder, zfitter, zpicker, zpicker2
from redmonster.sandbox import yanny as y
from redmonster.physics import misc
from time import gmtime, strftime
import matplotlib.pyplot as p
p.interactive(True)

''' Set plate, mjd, and fibers to be run.  If fiberid is not specified here and subsequently passed in during the next step, the default behavior is to run on all fibers. '''
plate = 7848
mjd = 56959
fiberid = [i for i in xrange(1000)] # fiberid must be a list, not a numpy array


''' Read spPlate file.  specs.flux, specs.ivar, specs.loglambda, are [nfibers, npix] arrays containing flux, inverse variances, and log-wavelength, respectively.  This step also flags sky fibers and masks pixels with unreasonable S/N. '''

specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)

#specs = donna_spec.Spec('/uufs/astro.utah.edu/common/home/u0814744/test.fits')
#skyfibers = n.where(specs.ebt1 == 0)[0][0:2]
#specs.flux = specs.flux[skyfibers]
#specs.ivar = specs.ivar[skyfibers]
#specs.dof = specs.dof[skyfibers]

''' Instantiate zfinder object that will do z-finding for the entire plate using a single template.  Here, fname is the template filename, assumed to be in $REDMONSTER_DIR/templates/ . npoly specifies number of polynomial terms to be used in finding, zmin and zmax are upper and lower bounds of the redshift range to be explored. Optionally, npixstep can specify the width of pixel steps in doing the cross-correlation.  If left blank, it defaults to 1.'''

zssp1 = zfinder.Zfinder(fname='ndArch-ssp_galaxy_cont-v002.fits', npoly=4, zmin=-0.01, zmax=1.2)

''' Run actual fitting routine on the object created above. zssp.zchi2arr is array of of minimum chi^2 values of shape [nfibers, ndim_1, ndim_2, ..., ndim_N, nredshifts], where ndim_i is the i'th dimension of the template file.  Input flux need not be SDSS data - any spectra binned in constant log(lambda) will work.'''

zssp1.zchi2(specs.flux, specs.loglambda, specs.ivar, npixstep=2)

''' New objects and fitting for different templates. '''

zssp2 = zfinder.Zfinder(fname='ndArch-ssp_galaxy_emit-v002.fits', npoly=4, zmin=-0.01, zmax=1.2)
zssp2.zchi2(specs.flux, specs.loglambda, specs.ivar, npixstep=2)
zstar = zfinder.Zfinder(fname='ndArch-all-CAP-grids.fits', npoly=4, zmin=-.005, zmax=.005)
zstar.zchi2(specs.flux, specs.loglambda, specs.ivar)
zqso = zfinder.Zfinder(fname='ndArch-QSO-V003.fits', npoly=4, zmin=.4, zmax=3.5)
zqso.zchi2(specs.flux, specs.loglambda, specs.ivar, npixstep=4)

''' Instantiate Zfitter to do subgrid fitting.  zchi2_ssp is chi^2 array from zfinder object above, and zssp.zbase is redshift-pixel baseline over the range explored by zfinder. '''

zfit_ssp = zfitter.Zfitter(zssp.zchi2arr, zssp.zbase)

''' Do actual subgrid refinement and fitting.  Best-fit and second-best-fit redshifts will be in [nfibers,2] shaped array in zfit_*****.z , and associated errors are in zfit_*****.z_err .  This routine also flags for small delta chi^2 and if the global minimum chi^2 is on an edge of the input chi^2 array. z_refine() includes an optional parameter, 'width', which specifies half the width of a window around the global chi2 minimum in which local minima are ignored during calculation of second-best redshift and flagging small delta-chi2.  If not specified, 'width' defaults to 15.'''

zfit_ssp.z_refine2()

''' Same as above for second template. '''

zfit_star = zfitter.Zfitter(zstar.zchi2arr, zstar.zbase)
zfit_star.z_refine2()
zfit_qso = zfitter.Zfitter(zqso.zchi2arr, zqso.zbase)
zfit_qso.z_refine2()

''' Flagging throughout redmonster is done individually by the classes responsible for handling the relevant computations.  To have an 'overall' flag for each fiber, the individual flags need to be combined. '''

ssp_flags = misc.comb_flags(specs, zssp, zfit_ssp)
star_flags = misc.comb_flags(specs, zstar, zfit_star)
qso_flags = misc.comb_flags(specs, zqso, zfit_qso)

''' Compare chi2 surfaces from each template and classify each object accordingly. Arguments are data object (in a format identical to that created by Spec), followed by each object created by Zfinder, Zfitter, and flags, in that order.  This function can currently handle up to five objects from five separate templates. If specs is a user created data object rather than one created by redmonster.datamgr.spec, it must contain specs.npix, the number of pixels in a single spectrum.'''

#zpick = zpicker.Zpicker(specs, zssp, zfit_ssp, ssp_flags, zstar, zfit_star, star_flags, zqso, zfit_qso, qso_flags)
zfindobjs=[]
zfindobjs.append(zssp)
zfindobjs.append(zstar)
zfindobjs.append(zqso)
zfitobjs=[]
zfitobjs.append(zfit_ssp)
zfitobjs.append(zfit_star)
zfitobjs.append(zfit_qso)
flags = []
flags.append(ssp_flags)
flags.append(star_flags)
flags.append(qso_flags)

zpick = zpicker2.Zpicker(specs, zfindobjs, zfitobjs, flags)

#zpick.plate = 0000
#zpick.mjd = 00000
#zpick.fiberid = [0]

''' Write output file.  Arguments are zpick object from above, and optionally dest and clobber, the path in which to write to file and whether or not to clobber old files with the same name, respectively.  See class documentation for more detail on Write_Redmonster behavior.'''

output = io2.Write_Redmonster(zpick)
output.write_fiber()

# Things left to do
#
# DONE 1. Incorporate flags into, probably, zpicker and subseqently Write_Redmonster
# 2. Incorporate Adam's spCFrame fittings somewhere
# 3. Function to turn variable resolution data into coadded log(lambda) data?
# DONE 4. Quasar templates?
# DONE 5. BOSS CMASS data testing
# DONE      Look into delta chi2 thresholds to trigger failure
# DONE      Look at completeness vs. purity
# 6. WISE selected LRG targets (SEQUELS/SDSS-IV)
