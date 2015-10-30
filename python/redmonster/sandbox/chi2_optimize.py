import numpy as n
import matplotlib.pyplot as p
p.interactive(True)
from redmonster.sandbox import yanny as y
from astropy.io import fits
from redmonster.datamgr import spec, io
from redmonster.physics import zfinder, zfitter, zpicker
from redmonster.physics import misc
from time import gmtime, strftime
import multiprocessing as mp


startfib = 100
nfib = 1

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

completeness = []
purity = []

#-------------------------------------------------------------------------------------------------------


# Loop over various dchi2 thresholds
for inc in xrange(2):
    this_thresh = 41 + inc
    threshold_vals.append(this_thresh)
    specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)
    zssp = zfinder.ZFinder(fname='ndArch-ssp_em_galaxy-v000.fits', type='GALAXY', npoly=4, zmin=-0.01, zmax=1.2)
    zssp.zchi2(specs.flux, specs.loglambda, specs.ivar)
    zstar = zfinder.ZFinder(fname='ndArch-spEigenStar-55734.fits', type='STAR', npoly=4, zmin=-.005, zmax=.005)
    zstar.zchi2(specs.flux, specs.loglambda, specs.ivar)
    zfit_ssp = zfitter.ZFitter(zssp.zchi2arr, zssp.zbase)
    zfit_ssp.z_refine(threshold=this_thresh)
    zfit_star = zfitter.ZFitter(zstar.zchi2arr, zstar.zbase)
    zfit_star.z_refine(threshold=this_thresh)
    zpick = zpicker.ZPicker(specs, zssp, zfit_ssp, zstar, zfit_star)
    ssp_flags = n.zeros(len(fiberid))
    star_flags = n.zeros(len(fiberid))
    for ifiber in xrange(len(fiberid)):
        ssp_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zssp.zwarning[ifiber])) | int(zfit_ssp.zwarning[ifiber])
        star_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zstar.zwarning[ifiber])) | int(zfit_star.zwarning[ifiber])
    
    purity.append( (len(n.where(ssp_flags == 0)))/float(nfib) )
    completeness.append( (len(n.where(abs(zpick.z[:,0]-zperson) <= .0001)))/float(nfib) )

#-------------------------------------------------------------------------------------------------------


#threshold_vals = [41,42]
threshold_vals = [35.+i for i in xrange(10)]

def run_redmonster(this_thresh):
    #global completeness
    #global purity
    global nfib
    completeness = []
    purity = []
    specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)
    zssp = zfinder.ZFinder(fname='ndArch-ssp_em_galaxy-v000.fits', npoly=4, zmin=-0.01, zmax=1.2)
    zssp.zchi2(specs.flux, specs.loglambda, specs.ivar)
    zstar = zfinder.ZFinder(fname='ndArch-spEigenStar-55734.fits', npoly=4, zmin=-.005, zmax=.005)
    zstar.zchi2(specs.flux, specs.loglambda, specs.ivar)
    zfit_ssp = zfitter.ZFitter(zssp.zchi2arr, zssp.zbase)
    zfit_ssp.z_refine(threshold=this_thresh)
    zfit_star = zfitter.ZFitter(zstar.zchi2arr, zstar.zbase)
    zfit_star.z_refine(threshold=this_thresh)
    zpick = zpicker.ZPicker(specs, zssp, zfit_ssp, zstar, zfit_star)
    ssp_flags = n.zeros(len(fiberid))
    star_flags = n.zeros(len(fiberid))
    for ifiber in xrange(len(fiberid)):
        ssp_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zssp.zwarning[ifiber])) | int(zfit_ssp.zwarning[ifiber])
        star_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zstar.zwarning[ifiber])) | int(zfit_star.zwarning[ifiber])
    #purity.append( (len(n.where(ssp_flags == 0)))/float(nfib) )
    #completeness.append( (len(n.where(abs(zpick.z[:,0]-zperson) <= .0001)))/float(nfib) )
    purity = (len(n.where(ssp_flags == 0)))/float(nfib)
    completeness = (len(n.where(abs(zpick.z[:,0]-zperson) <= .0001)))/float(nfib)
    return completeness, purity


# Without multi
print strftime("%Y-%m-%d %H:%M:%S", gmtime())
result = map(run_redmonster,threshold_vals)
print strftime("%Y-%m-%d %H:%M:%S", gmtime())

# With multi
print strftime("%Y-%m-%d %H:%M:%S", gmtime())
num_proc = 12
pool = mp.Pool(num_proc)
completeness, purity = pool.map(run_redmonster,threshold_vals)
print strftime("%Y-%m-%d %H:%M:%S", gmtime())

completeness = n.array(completeness)
purity = n.array(purity)




































