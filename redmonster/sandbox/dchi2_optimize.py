# Optimize dchi2 threshold in zfitter for best purity/completeness

import numpy as n
import matplotlib.pyplot as p
p.interactive(True)
from redmonster.sandbox import yanny as y
from astropy.io import fits
from redmonster.datamgr import spec, io
from redmonster.physics import zfinder, zfitter, zpicker
from redmonster.math import misc
from time import gmtime, strftime

# Read yanny file
x = y.yanny(filename='spInspect_alltest_bolton.par.txt', np=True)

# Get fibers, zpipe, zperson for each plate
args = n.where(x['BOSSOBJECT']['plate'] == 3686)[0]
fibers3686 = []
zpipe3686 = []
zperson3686 = []
for i in args:
    fibers3686.append( x['BOSSOBJECT'][i][2]-1)
    zpipe3686.append( x['BOSSOBJECT'][i][5])
    zperson3686.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3687)[0]
fibers3687 = []
zpipe3687 = []
zperson3687 = []
for i in args:
    fibers3687.append( x['BOSSOBJECT'][i][2]-1 )
    zpipe3687.append( x['BOSSOBJECT'][i][5])
    zperson3687.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3804)[0]
fibers3804 = []
zpipe3804 = []
zperson3804 = []
for i in args:
    fibers3804.append( x['BOSSOBJECT'][i][2]-1)
    zpipe3804.append( x['BOSSOBJECT'][i][5])
    zperson3804.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3805)[0]
fibers3805 = []
zpipe3805 = []
zperson3805 = []
for i in args:
    fibers3805.append( x['BOSSOBJECT'][i][2]-1)
    zpipe3805.append( x['BOSSOBJECT'][i][5])
    zperson3805.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3853)[0]
fibers3853 = []
zpipe3853 = []
zperson3853 = []
for i in args:
    fibers3853.append( x['BOSSOBJECT'][i][2]-1)
    zpipe3853.append( x['BOSSOBJECT'][i][5])
    zperson3853.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3855)[0]
fibers3855 = []
zpipe3855 = []
zperson3855 = []
for i in args:
    fibers3855.append( x['BOSSOBJECT'][i][2]-1)
    zpipe3855.append( x['BOSSOBJECT'][i][5])
    zperson3855.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3856)[0]
fibers3856 = []
zpipe3856 = []
zperson3856 = []
for i in args:
    fibers3856.append( x['BOSSOBJECT'][i][2]-1)
    zpipe3856.append( x['BOSSOBJECT'][i][5])
    zperson3856.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3860)[0]
fibers3860 = []
zpipe3860 = []
zperson3860 = []
for i in args:
    fibers3860.append( x['BOSSOBJECT'][i][2]-1)
    zpipe3860.append( x['BOSSOBJECT'][i][5])
    zperson3860.append( x['BOSSOBJECT'][i][6])


args1= [(3686, 55268, fibers3686, zperson3686),(3687, 55269, fibers3687, zperson3687),(3804, 55267, fibers3804, zperson3804),(3805, 55269, fibers3805, zperson3805),(3853, 55268, fibers3853, zperson3853),(3855, 55268, fibers3855, zperson3855),(3856, 55269, fibers3856, zperson3856),(3860, 55269, fibers3860, zperson3860)]

#threshold_vals = [41,42]
threshold_vals = [35.+i for i in xrange(10)]

completeness = []
purity = []

this_thresh = 46.6

'''
print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
thiscomp, thispur = find_comp_purity(this_thresh, [args1[0]])
print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
'''

# ----------------------------------------------------------------------------------------------------------

def find_comp_purity(this_thresh, args):
    purity = []
    completeness = []
    for iarg in args:
        plate = iarg[0]
        mjd = iarg[1]
        fiberid = iarg[2]
        zperson = iarg[3]
        import pdb; pdb.set_trace()
        specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)
        
        hdu = fits.open('/uufs/astro.utah.edu/common/home/u0814744/scratch/screens/chi2arr-%s-ssp_em_galaxy.fits' % plate)
        #hdu = fits.open('/Users/boltonlab3/scratch/chi2arr-%s-ssp_em_galaxy.fits' % plate)
        sspchi2arr = hdu[0].data
        zbasessp = hdu[1].data.ZBASE
        
        hdu = fits.open('/uufs/astro.utah.edu/common/home/u0814744/scratch/screens/chi2arr-%s-spEigenStar.fits' % plate)
        #hdu = fits.open('/Users/boltonlab3/scratch/chi2arr-%s-spEigenStar.fits' % plate)
        starchi2arr = hdu[0].data
        zbasestar = hdu[1].data.ZBASE
        
        zfit_ssp = zfitter.Zfitter(sspchi2arr, zbasessp)
        zfit_ssp.z_refine(threshold=this_thresh)
        
        zfit_star = zfitter.Zfitter(starchi2arr, zbasestar)
        zfit_star.z_refine(threshold=this_thresh)
        
        ssp_flags = misc.comb_flags_2(specs, zfit_ssp.zwarning)
        star_flags = misc.comb_flags_2(specs, zfit_star.zwarning)
        
        zpick = Hacked_zpicker(specs, sspchi2arr, zfit_ssp, ssp_flags, starchi2arr, zfit_star, star_flags)
        
        purity.append( (len(n.where(zpick.zwarning == 0)))/float(len(fiberid)) )
        
        completeness.append( (len(n.where(abs(zpick.z[:,0]-zperson) <= .0001)))/float(len(fiberid)) )

    this_comp = n.mean(completeness)
    this_pur = n.mean(purity)

    return this_comp, this_pur




# ----------------------------------------------------------------------------------------------------------

class Hacked_zpicker:
    
    def __init__(self, specobj, zchi2arr1, zfit1, flags1, zchi2arr2=None, zfit2=None, flags2=None):
        self.npixflux = specobj.npix
        self.type = []
        self.minvector = []
        self.zwarning = []
        self.z = n.zeros( (zchi2arr1.shape[0],2) )
        if zchi2arr2 == None: self.nclass = 1
        else: self.nclass = 2
        self.minrchi2 = n.zeros( (zchi2arr1.shape[0], self.nclass) )
        self.classify_obj(zchi2arr1, zfit1, flags1, zchi2arr2, zfit2, flags2)
    
    def classify_obj(self, zchi2arr1, zfit1, flags1, zchi2arr2, zfit2, flags2):
        flag_val = int('0b100',2) # From BOSS zwarning flag definitions
        for ifiber in xrange(zchi2arr1.shape[0]):
            self.minrchi2[ifiber,0] = n.min(zchi2arr1[ifiber]) / (self.npixflux) # Calculate reduced chi**2 values to compare templates of differing lengths
            if zchi2arr2 != None: self.minrchi2[ifiber,1] = n.min(zchi2arr2[ifiber]) / (self.npixflux)
            minpos = self.minrchi2[ifiber].argmin() # Location of best chi2 of array of best (individual template) chi2s
            
            if minpos == 0: # Means overall chi2 minimum came from template 1
                self.type.append('GALAXY')
                self.minvector.append(zfit1.minvector[ifiber])
                minloc = n.unravel_index(zchi2arr1[ifiber].argmin(), zchi2arr1[ifiber].shape)[:-1]
                self.zwarning = n.append(self.zwarning, flags1[ifiber])
                argsort = self.minrchi2[ifiber].argsort()
                if len(argsort) > 1:
                    if argsort[1] == 1:
                        if ( n.min(zchi2arr2[ifiber]) - n.min(zchi2arr1[ifiber]) ) < zfit1.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val #THIS THRESHOLD PROBABLY ISN'T RIGHT AND NEEDS TO BE CHANGED
            
            elif minpos == 1: # Means overall chi2 minimum came from template 2
                self.type.append('STAR')
                self.minvector.append(zfit2.minvector[ifiber])
                minloc = n.unravel_index(zchi2arr2[ifiber].argmin(), zchi2arr2[ifiber].shape)[:-1]
                self.zwarning = n.append(self.zwarning, flags2[ifiber])
                argsort = self.minrchi2[ifiber].argsort()
                if argsort[1] == 0:
                    if ( n.min(zchi2arr1[ifiber]) - n.min(zchi2arr2[ifiber]) ) < zfit2.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val































