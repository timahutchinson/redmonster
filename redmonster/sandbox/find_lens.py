# Attempt at finding gravitational lenses using chi2 surfaces from zfinder to find superpositions
# of galaxy spectra at different redshifts.
#
# Arguments:
#   specobj = data object created by redmonster.datamgr.spec (or a similar user-created object)
#   chi2arr = chi**2 surface in shape of (nfiber,npar1,...,nparN,nredshift), designed to be from zfinder
#   dof = array of degrees of freedom for each fiber
#   fname = name of template file, assumed to be in $REDMONSTER_DIR/templates/
#   threshold = only fibers with chi2min > (threshold*dof) will be tested
#   width = number of pixels on either side of global minimum chi2 that will be ignored when searching
#           for a secondary minimum
#   npoly = number of polynomial terms (up to degree n-1) to fit alongside templates.  This should
#           match the number used when creating chi2arr
#
# Tim Hutchinson, University of Utah, September 2014
# t.hutchinson@utah.edu
#
# Thought: perhaps this would be better if, instead of some absolute threshold on min(chi2) based
# on the dof, I could find something like an average chi2 across the fiber, and then have the
# threshold be some amount below that.  That would probably better handle fibers with really high S/N?

import numpy as n
import matplotlib.pyplot as p
p.interactive(True)
from os import environ
from os.path import join, exists
from astropy.io import fits
from redmonster.math.misc import poly_array
from time import gmtime, strftime

class Find_Lens:

    def __init__(self, specobj, chi2arr, dof, fname, threshold=1.5, width=15, npoly=4):
        self.flux = specobj.flux
        self.ivar = specobj.ivar
        self.dof = dof
        self.fname = fname
        self.threshold = threshold
        self.width = width
        self.npoly = npoly
        self.bestchi2vecs = n.zeros( (chi2arr.shape[0],chi2arr.shape[-1]) )
        self.chi2mins = n.zeros( chi2arr.shape[0] )
        for i in xrange(chi2arr.shape[0]):
            print 'Checking fiber #' + str(i+1) + ' of ' + str(chi2arr.shape[0])
            for j in xrange(chi2arr.shape[-1]):
                self.bestchi2vecs[i,j] = n.min(chi2arr[i,...,j])
            self.chi2mins[i] = n.min(self.bestchi2vecs[i])
        self.checklocs = n.where( self.chi2mins >= (float(threshold)*dof) )[0]
        print '%s fibers have possible lenses!' % len(self.checklocs)
        print n.asarray(specobj.fiberid)[n.asarray(self.checklocs)]
        chi2s = n.zeros( (len(self.checklocs),300,300) )
        # For each fiber in checklocs, refit a linear combination of templates at best and second best z, separated from global min by width
        ind = 0
        for loc in self.checklocs:
            #print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
            print 'THIS INDEX NUMBER IS ' + str(ind) + '!!!!!!!!!!!!!!!!!!!!!!!!!'
            chi2s[ind] = self.fit_linear_combo(loc)
            ind += 1
            #print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing

    def fit_linear_combo(self, loc):
        # Find locations of best and second best z from chi2 surfaces
        chi2vec = self.bestchi2vecs[loc]
        minloc = chi2vec.argmin()
        if (minloc > self.width) & (minloc < (chi2vec.shape[0]-self.width)):
            minlessarg = chi2vec[:minloc-self.width].argmin()
            minmorearg = chi2vec[minloc+self.width:].argmin() + (minloc+self.width) # Argmin gives index within the slice; a value of 2 from chi2vec[minloc:].argmin() means minloc+2 was the minimum
            secondminloc = minlessarg if (chi2vec[minlessarg] < chi2vec[minmorearg]) else minmorearg
        elif (minloc < self.width):
            secondminloc = chi2vec[minloc+self.width:].argmin() + (minloc+width)
        else:
            secondminloc = chi2vec[:minloc-self.width].argmin()
        # Read in template
        temps = fits.open(join(environ['REDMONSTER_DIR'],'templates',eval(self.fname[loc])))[0].data
        flattemps = n.reshape(temps, (-1,temps.shape[-1]))
        fitchi2s = n.zeros( (flattemps.shape[0],flattemps.shape[0]) )
        data = self.flux[loc]
        ivar = self.ivar[loc]
        sn2_data = n.sum( (data**2)*ivar )
        polyarr = poly_array(self.npoly,data.shape[0])
        # Start setting up matrices for fitting
        pmat = n.zeros( (data.shape[0],self.npoly+2) )
        pmat[:,2], pmat[:,3], pmat[:,4], pmat[:,5] = polyarr[0],polyarr[1],polyarr[2],polyarr[3]
        ninv = n.diag(ivar)
        for i in xrange(flattemps.shape[0]): #changed from flattemps.shape[0]
            print i
            #import pdb; pdb.set_trace()
            pmat[:,0] = flattemps[i,minloc:minloc+data.shape[0]]
            for j in xrange(flattemps.shape[0]): #changed from flattemps.shape[0]
                pmat[:,1] = flattemps[j,secondminloc:secondminloc+data.shape[0]]
                pmattrans = n.transpose(pmat)
                a = n.dot(n.dot(pmattrans,ninv),pmat)
                b = n.dot(n.dot(pmattrans,ninv),data)
                f = n.linalg.solve(a,b) # Maybe this should use nnls instead of solve to preserve physicality?
                fitchi2s[i,j] = sn2_data - n.dot(n.dot(f,a),f)
        return fitchi2s
#import pdb; pdb.set_trace()


#2014-09-19 21:11:23
