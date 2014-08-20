import numpy as n
import matplotlib.pyplot as p
p.interactive(True)
from redmonster.sandbox import yanny as y
from astropy.io import fits
from redmonster.datamgr import spec, io
from redmonster.physics import zfinder, zfitter, zpicker
from time import gmtime, strftime


startfib = 100
nfib = 50

# Read yanny file
x = y.yanny(filename='/Users/boltonlab3/Downloads/spInspect_alltest_bolton.par.txt', np=True)

# Get fibers, zpipe, zperson for each plate
args = n.where(x['BOSSOBJECT']['plate'] == 3686)[0]
fibers3686 = []
zpipe3686 = []
zperson3686 = []
for i in args:
    fibers3686.append( x['BOSSOBJECT'][i][2])
    zpipe3686.append( x['BOSSOBJECT'][i][5])
    zperson3686.append( x['BOSSOBJECT'][i][6])

# Only keep 50 fibers
fiberid = fibers3686[startfib:startfib+nfib]
zpipe = zpipe3686[startfib:startfib+nfib]
zperson = zperson3686[startfib:startfib+nfib]

plate = 3686
mjd = 55268

threshold_vals = []
completeness = []
purity = []

# Loop over various dchi2 thresholds
for inc in xrange(10):
    this_thresh = 41.6 + inc
    threshold_vals.append(this_thresh)
    specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)
    zssp = zfinder.Zfinder(fname='ndArch-ssp_em_galaxy-v000.fits', type='GALAXY', npoly=4, zmin=-0.01, zmax=1.2)
    zssp.zchi2(specs.flux, specs.loglambda, specs.ivar)
    zstar = zfinder.Zfinder(fname='ndArch-spEigenStar-55734.fits', type='STAR', npoly=4, zmin=-.005, zmax=.005)
    zstar.zchi2(specs.flux, specs.loglambda, specs.ivar)
    zfit_ssp = zfitter.Zfitter(zssp.zchi2arr, zssp.zbase)
    zfit_ssp.z_refine(threshold=this_thresh)
    zfit_star = zfitter.Zfitter(zstar.zchi2arr, zstar.zbase)
    zfit_star.z_refine(threshold=this_thresh)
    zpick = zpicker.Zpicker(specs, zssp, zfit_ssp, zstar, zfit_star)
    ssp_flags = n.zeros(len(fiberid))
    star_flags = n.zeros(len(fiberid))
    for ifiber in xrange(len(fiberid)):
        ssp_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zssp.zwarning[ifiber])) | int(zfit_ssp.zwarning[ifiber])
        star_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zstar.zwarning[ifiber])) | int(zfit_star.zwarning[ifiber])
    purity.append( (len(n.where(ssp_flags == 0)))/float(nfib) )
    completeness.append( (len(n.where(abs(zpick.z[:,0]-zperson) <= .0001)))/float(nfib) )