# Make a plot of the "best fit" spectra to any object, of arbitrary
# template type.  Only works for a single fiber.
#
# Tim Hutchinson, University of Utah, Oct 2014
# t.hutchinson@utah.edu

from time import gmtime, strftime

import numpy as n
from astropy.io import fits
import matplotlib.pyplot as p
p.interactive(True)
from astropy.convolution import convolve, Box1DKernel

from redmonster.datamgr import spec, io
from redmonster.physics import zfinder, zfitter, zpicker
from redmonster.sandbox import yanny as y
from redmonster.physics import misc
from redmonster.physics.misc import poly_array


plate = 3686 # Set plate, mjd, and fiberid here
mjd = 55268
fiberid = [89,100,102] #935, 937

specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)

# Use Charlie's SSPs
ztemp = zfinder.ZFinder(fname='ndArch-ssp_em_galaxy-v000.fits', npoly=4,
                        zmin=-0.01, zmax=1.2)
# Use Carlos' stellar templates
#ztemp = zfinder.ZFinder(fname='ndArch-all-CAP-grids.fits', npoly=4,
                         #zmin=-.005, zmax=.005)
# Use spEigenstars from IDL pipeline
#ztemp = zfinder.ZFinder(fname='ndArch-spEigenStar-55734.fits', npoly=4,
                         #zmin=-.005, zmax=.005)
# Use Nao's quasars
#ztemp = zfinder.ZFinder(fname='ndArch-QSO-V003.fits', npoly=4, zmin=.4,
                         #zmax=3.5)

ztemp.zchi2(specs.flux, specs.loglambda, specs.ivar, npixstep=1)
zfit_temp = zfitter.ZFitter(ztemp.zchi2arr, ztemp.zbase)
zfit_temp.z_refine()
#temp_flags = misc.comb_flags(specs, ztemp, zfit_temp)
#zpick = zpicker.ZPicker(specs, ztemp, zfit_temp)

# Solve for parameters, create model
import pdb; pdb.set_trace()
minloc = n.unravel_index( ztemp.zchi2arr.argmin(), ztemp.zchi2arr.shape )
pmat = n.zeros( (specs.flux.shape[-1],ztemp.npoly+1) )
this_temp = ztemp.templates[minloc[1:-1]]
pmat[:,0] = this_temp[(minloc[-1]*ztemp.npixstep)+ztemp.pixoffset:
                      (minloc[-1]*ztemp.npixstep)+ztemp.pixoffset+
                      specs.flux.shape[-1]]
polyarr = poly_array(ztemp.npoly, specs.flux.shape[-1])
pmat[:,1:] = n.transpose(polyarr)
ninv = n.diag(specs.ivar[0])
f = n.linalg.solve( n.dot(n.dot(n.transpose(pmat),ninv),pmat),
                   n.dot( n.dot(n.transpose(pmat),ninv),specs.flux[0]) )
model = n.dot(pmat, f)

# Make plot
p.plot(10**specs.loglambda, specs.flux[0],'r', label='Data')
p.plot(10**specs.loglambda, model, 'k', label='Model', hold=True)
p.title('Plate %s Fiber %s, z=%.4f' % (plate, fiberid[0], zfit_temp.z[0][0]),
        size=18)
p.xlabel(r'Wavelength ($\AA$)', size=16)
p.ylabel(r'Flux ($10^{-17}$ erg s$^{-1}$cm$^{-2}$$\AA^{-1}$)', size=16)
p.legend()



# Optionally, smooth data and then plot
'''
smoothed = convolve(specs.flux[0], Box1DKernel(5))
p.plot(10**specs.loglambda, smoothed,'r', label='Data')
p.plot(10**specs.loglambda, model, 'k', label='Model', hold=True)
p.title('Plate %s Fiber %s, z=%.4f' % (plate, fiberid[0], zfit_temp.z[0][0]), 
        size=18)
p.xlabel(r'Wavelength ($\AA$)', size=16)
p.ylabel(r'Flux ($10^{-17}$ erg s$^{-1}$cm$^{-2}$$\AA^{-1}$)', size=16)
p.legend()
'''












