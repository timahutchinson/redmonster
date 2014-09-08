# Optimize dchi2 threshold in zfitter for best purity/completeness

import numpy as n
import matplotlib.pyplot as p
p.interactive(True)
from redmonster.sandbox import yanny as y
from astropy.io import fits
from redmonster.datamgr import spec, io
from redmonster.physics import zfinder, zfitter, zpicker
from redmonster.math import misc
from time import gmtime, strftime
from os.path import join
from os import environ
from scipy.optimize import curve_fit
from os.path import join
from os import environ

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



# Remove fiber numbers from each plate that have zperson = -9.
args = []
for i in xrange(8):
    badzfibs = n.where(n.asarray(args1[i][3]) == -9.)[0]
    keepfibs = n.delete(n.asarray(args1[i][2]),badzfibs).tolist()
    keepzperson = n.delete(n.asarray(args1[i][3]),badzfibs).tolist()
    args.append( (args1[i][0],args1[i][1],keepfibs,keepzperson) )

numgals = 0
for i in xrange(8):
    numgals += len(args[i][2])
print 'Total number of galaxies:' + str(numgals)

threshold_vals = [5.+(.2*i) for i in xrange(100)]
completeness = []
purity = []

# ----------------------------------------------------------------------------------------------------------

def find_comp_purity(this_thresh, args):
    purity = []
    completeness = []
    run = 1
    for iarg in args:
        print 'Running plate %s of 8' % (run)
        run += 1
        plate = iarg[0]
        mjd = iarg[1]
        fiberid = iarg[2]
        zperson = iarg[3]
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
        
        #completeness.append( (len(n.where(zpick.zwarning == 0)[0]))/float(len(fiberid)) )
        completeness.append( (len(n.where((zpick.zwarning.astype(int) & 4) == 0)[0]))/float(len(fiberid)) )
        
        purity_set = zpick.z[n.where((zpick.zwarning.astype(int) & 4) == 0)[0],0]
        purity_zperson = n.asarray(zperson)[n.where((zpick.zwarning.astype(int) & 4) == 0)[0]]
        #purity.append( (len(n.where(abs(zpick.z[:,0]-zperson) <= .0005)[0]))/float(len(fiberid)) )
        purity.append( (len(n.where(abs(purity_set-purity_zperson) <= .0005)[0]))/float(len(purity_set)) )

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
        self.z_err = n.zeros( (zchi2arr1.shape[0],2) )
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
                self.z[ifiber] = zfit1.z[ifiber]
                self.z_err[ifiber] = zfit1.z_err[ifiber]
                minloc = n.unravel_index(zchi2arr1[ifiber].argmin(), zchi2arr1[ifiber].shape)[:-1]
                self.zwarning = n.append(self.zwarning, flags1[ifiber])
                argsort = self.minrchi2[ifiber].argsort()
                if len(argsort) > 1:
                    if argsort[1] == 1:
                        if ( n.min(zchi2arr2[ifiber]) - n.min(zchi2arr1[ifiber]) ) < zfit1.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val #THIS THRESHOLD PROBABLY ISN'T RIGHT AND NEEDS TO BE CHANGED
            
            elif minpos == 1: # Means overall chi2 minimum came from template 2
                self.type.append('STAR')
                self.minvector.append(zfit2.minvector[ifiber])
                self.z[ifiber] = zfit2.z[ifiber]
                self.z_err[ifiber] = zfit2.z_err[ifiber]
                minloc = n.unravel_index(zchi2arr2[ifiber].argmin(), zchi2arr2[ifiber].shape)[:-1]
                self.zwarning = n.append(self.zwarning, flags2[ifiber])
                argsort = self.minrchi2[ifiber].argsort()
                if argsort[1] == 0:
                    if ( n.min(zchi2arr1[ifiber]) - n.min(zchi2arr2[ifiber]) ) < zfit2.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val

#--------------------------------------

print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
threshnum = 1
for this_thresh in threshold_vals:
    print 'Running threshold %s of %s'% (threshnum,len(threshold_vals))
    threshnum += 1
    thiscomp, thispur = find_comp_purity(this_thresh, args)
    completeness.append(thiscomp)
    purity.append(thispur)
print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing     


# Write output to fits file
prihdu = fits.PrimaryHDU()
col1 = fits.Column(name='COMPLETENESS', format='E', array=n.asarray(completeness))
col2 = fits.Column(name='PURITY', format='E', array=n.asarray(purity))
col3 = fits.Column(name='THRESHOLDS', format='E', array=n.asarray(threshold_vals))
cols = fits.ColDefs([col1,col2,col3])
tbhdu = fits.new_table(cols)
thdulist = fits.HDUList([prihdu,tbhdu])
thdulist.writeto('/uufs/astro.utah.edu/common/home/u0814744/scratch/comp_purity_5-25.fits', clobber=True)



'''
# Read and concatenate files
comp = n.array([])
pur = n.array)[])
thresh = n.array([])
comp = n.append(comp, fits.open('comp_purity_5-25.fits')[1].data.COMPLETENESS)
pur = n.append(pur, fits.open('comp_purity_5-25.fits')[1].data.PURITY)
thresh = n.append(thresh, fits.open('comp_purity_5-25.fits')[1].data.THRESHOLDS)
comp = n.append(comp, fits.open('comp_purity_25-45.fits')[1].data.COMPLETENESS)
pur = n.append(pur, fits.open('comp_purity_25-45.fits')[1].data.PURITY)
thresh = n.append(thresh, fits.open('comp_purity_25-45.fits')[1].data.THRESHOLDS)
comp = n.append(comp, fits.open('comp_purity_45-65.fits')[1].data.COMPLETENESS)
pur = n.append(pur, fits.open('comp_purity_45-65.fits')[1].data.PURITY)
thresh = n.append(thresh, fits.open('comp_purity_45-65.fits')[1].data.THRESHOLDS)

# Scatter plot of completeness vs purity
p.scatter(pur,comp,c=thresh)
p.plot(1,1,'rx',label='Ideal')
p.axis([.945,1.005,.93,1.005])
p.xlabel('Purity',size=14)
p.ylabel('Completeness',size=14)
p.title('Purity vs. Completeness for 4864 CMASS Galaxies',size=16)
cbar = p.colorbar()
cbar.set_label(r'$\delta \chi^2$ Threshold',size=14)
p.gca().figure.canvas.draw()
p.legend(loc=4)
print 'Completeness at optimal point is ' + str(comp[dist.argmin()])
print 'Purity at optimal point is ' + str(pur[dist.argmin()])
# Optional connecting line
def func(x,a,b):
    return a*x+b

xdata = n.array([pur[dist.argmin()],1.0])
ydata=n.array([comp[dist.argmin()],1.])
popt,pcov = curve_fit(func,xdata,ydata)
xfit = n.linspace(xdata[0],xdata[1],100)
yfit = func(xfit,popt[0],popt[1])
p.plot(xfit,yfit,color='black',label='Shortest Distance',hold=True)
p.legend(loc=4)

# Plot of threshold vs distance
dist = n.sqrt( (1-comp)**2 + (1-pur)**2 )
p.plot(thresh,dist,'b.')
p.xlabel(r'$\delta \chi^2$ Threshold',size=14)
p.ylabel('Distance',size=14)
p.title(r'$\delta \chi^2$ Threshold vs. Distance from (1.0,1.0)',size=16)

# Zoom in around global minimum and fit quadratic, then plot
ydata = dist[dist.argmin()-2:dist.argmin()+4]
xdata = thresh[dist.argmin()-2:dist.argmin()+4]
def func(x,a,b,c):
    return a*(x**2)+b*x+c

popt,pcov = curve_fit(func,xdata,ydata)
xfit = n.linspace(xdata[0],xdata[-1],100)
yfit = func(xfit,popt[0],popt[1],popt[2])
p.plot(xdata,ydata,'ko',label='Data')
p.plot(xfit,yfit,color='red',label='Fit')
p.xlabel(r'$\delta \chi^2$ Threshold',size=14)
p.ylabel('Distance',size=14)
p.title(r'$\delta \chi^2$ Threshold vs. Distance from (1.0,1.0) Near Global Minimum',size=16)
p.legend()

# Use minimum of quadratic as 'best' overall dchi2 threshold
print 'Best dchi2 threshold is ' + str(xfit[yfit.argmin()])




# Make same plot for IDL outputs
thresh = [5+(.2*j) for j in xrange(300)]
pur_idl = []
comp_idl = []
for this_thresh in thresh:
    purity = []
    completeness = []
    print 'Thresh' + str(this_thresh)
    for i in xrange(8):
        plate = args[i][0]
        mjd = args[i][1]
        fibers = args[i][2]
        zperson = args[i][3]
        hdu = fits.open( join(environ['BOSS_SPECTRO_REDUX'],environ['RUN2D'],str(plate), environ['RUN1D'],'spZbest-%s-%s.fits' % (plate,mjd)) )
        dof = hdu[1].data.DOF[fibers]
        rchi2diff = hdu[1].data.RCHI2DIFF_NOQSO[fibers]
        z = hdu[1].data.Z_NOQSO[fibers]
        flags = hdu[1].data.ZWARNING_NOQSO[fibers]
        chi2diff = rchi2diff*dof
        completeness.append( len( n.where(chi2diff > this_thresh)[0]) / float(len(fibers)) )
        purity_set = z[n.where(chi2diff > this_thresh)[0]]
        purity_zperson = n.asarray(zperson)[n.where(chi2diff > this_thresh)[0]]
        purity.append( (len(n.where(abs(purity_set-purity_zperson) <= .0005)[0]))/float(len(purity_set)) )
    this_comp = n.mean(completeness)
    this_pur = n.mean(purity)
    comp_idl.append(this_comp)
    pur_idl.append(this_pur)
'''


















