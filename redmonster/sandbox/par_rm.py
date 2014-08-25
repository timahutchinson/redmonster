# Run redmonster on a cluster using multiprocessesing
#
# Tim Hutchinson, University of Utah, August 2014
# t.hutchinson@utah.edu


import numpy as n
import matplotlib.pyplot as p
p.interactive(True)
from redmonster.sandbox import yanny as y
from astropy.io import fits
from redmonster.datamgr import spec, io
from redmonster.physics import zfinder, zfitter, zpicker
from time import gmtime, strftime
import multiprocessing as mp

def parallel_rm( (plate,mjd,fiberid) ):
    specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)
    zssp = zfinder.Zfinder(fname='ndArch-ssp_em_galaxy-v000.fits', npoly=4, zmin=-0.01, zmax=1.2)
    zssp.zchi2(specs.flux, specs.loglambda, specs.ivar)
    # Write chi2 file with zbase
    prihdu = fits.PrimaryHDU(zssp.zchi2arr)
    col1 = fits.Column(name='ZBASE', format='E', array=zssp.zbase)
    cols = fits.ColDefs([col1])
    tbhdu = fits.new_table(cols)
    thdulist = fits.HDUList([prihdu,tbhdu])
    thdulist.writeto('chi2arr-%s-%s.fits' % (plate, zssp.type), clobber=True)
    # ----
    #zstar = zfinder.Zfinder(fname='ndArch-spEigenStar-55734.fits', npoly=4, zmin=-.005, zmax=.005)
    #zstar.zchi2(specs.flux, specs.loglambda, specs.ivar)
    #zfit_ssp = zfitter.Zfitter(zssp.zchi2arr, zssp.zbase)
    #zfit_ssp.z_refine(threshold=this_thresh)
    #zfit_star = zfitter.Zfitter(zstar.zchi2arr, zstar.zbase)
    #zfit_star.z_refine(threshold=this_thresh)
    #zpick = zpicker.Zpicker(specs, zssp, zfit_ssp, zstar, zfit_star)
    #ssp_flags = n.zeros(len(fiberid))
    #star_flags = n.zeros(len(fiberid))
    #for ifiber in xrange(len(fiberid)):
    #    ssp_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zssp.zwarning[ifiber])) | int(zfit_ssp.zwarning[ifiber])
    #    star_flags[ifiber] = (int(specs.zwarning[ifiber]) | int(zstar.zwarning[ifiber])) | int(zfit_star.zwarning[ifiber])
    ## Write flags file
    #prihdu = fits.PrimaryHDU(ssp_flags)
    #thdulist = fits.HDUList([prihdu])
    #thdulist.writeto('flags-%s-%s.fits' % (plate, zssp.type), clobber=True)
    #output = io.Write_Redmonster(zpick, dest='~/scratch', clobber=True)

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

args = n.where(x['BOSSOBJECT']['plate'] == 3687)[0]
fibers3687 = []
zpipe3687 = []
zperson3687 = []
for i in args:
    fibers3687.append( x['BOSSOBJECT'][i][2])
    zpipe3687.append( x['BOSSOBJECT'][i][5])
    zperson3687.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3804)[0]
fibers3804 = []
zpipe3804 = []
zperson3804 = []
for i in args:
    fibers3804.append( x['BOSSOBJECT'][i][2])
    zpipe3804.append( x['BOSSOBJECT'][i][5])
    zperson3804.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3805)[0]
fibers3805 = []
zpipe3805 = []
zperson3805 = []
for i in args:
    fibers3805.append( x['BOSSOBJECT'][i][2])
    zpipe3805.append( x['BOSSOBJECT'][i][5])
    zperson3805.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3853)[0]
fibers3853 = []
zpipe3853 = []
zperson3853 = []
for i in args:
    fibers3853.append( x['BOSSOBJECT'][i][2])
    zpipe3853.append( x['BOSSOBJECT'][i][5])
    zperson3853.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3855)[0]
fibers3855 = []
zpipe3855 = []
zperson3855 = []
for i in args:
    fibers3855.append( x['BOSSOBJECT'][i][2])
    zpipe3855.append( x['BOSSOBJECT'][i][5])
    zperson3855.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3856)[0]
fibers3856 = []
zpipe3856 = []
zperson3856 = []
for i in args:
    fibers3856.append( x['BOSSOBJECT'][i][2])
    zpipe3856.append( x['BOSSOBJECT'][i][5])
    zperson3856.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3860)[0]
fibers3860 = []
zpipe3860 = []
zperson3860 = []
for i in args:
    fibers3860.append( x['BOSSOBJECT'][i][2])
    zpipe3860.append( x['BOSSOBJECT'][i][5])
    zperson3860.append( x['BOSSOBJECT'][i][6])

args1= [(3686, 55268, fibers3686),(3687, 55269, fibers3687),(3804, 55267, fibers3804),(3805, 55269, fibers3805),(3853, 55268, fibers3853),(3855, 55268, fibers3855),(3856, 55269, fibers3856),(3860, 55269, fibers3860)]

'''
num_proc = 8
pool = mp.Pool(num_proc)
result = pool.map(parallel_rm,args1)
'''


















