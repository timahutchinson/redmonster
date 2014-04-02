import numpy as n
from os import environ
from os.path import join, exists
from redmonster.datamgr.ssp_prep import SSP_Prep
from redmonster.physics.misc import poly_array

class zfinder:

    def __init__(self, config=None, npoly=None):
        self.config = config
        try: self.specdir = environ['IDLSPEC2D_DIR']
        except: self.specdir = None
        if self.config.lower() == 'ssp': self.set_SSP(npoly=npoly)

    def set_SSP(self, npoly=None):
        self.npoly = npoly if npoly else 3
        ssp_stuff = SSP_Prep(velmin=100, velstep=100, nvel=3) # THIS MAY NOT BE THE BEST WAY TO DO THIS
        self.templates = ssp_stuff.specs

    def zchi2(self, specs, ivar):
        shape = (specs.shape[0],) + self.templates.shape[:-1]
        zchi2arr = n.zeros(shape) # Create chi2 array of shape (# of fibers, template_parameter1, ... , template_parameterN)
        polyarr = poly_array(self.npoly, specs.shape[1]) # Compute poly terms, noting that they will stay fixed with the data - assumes data is passed in as shape (nfibers, npix)
        num_z = self.templates.shape[-1] - specs.shape[-1] + 1 # Number of pixels to be fitted in redshift
        bvec = n.zeros(shape=(npoly+1,num_z))
        amat = n.zeros(shape=(npoly+1,npoly+1,num_z)) # Amat(z).f=bvec(z)
        for i in range(specs.shape[0]): # Loop over fibers
            for j in range(self.template.shape[0]): # Loop over template-ages
                for k in range(self.template.shape[1]): # Loop over template vdisps
                    bvec[0] = n.convolve( spec[i]/(ivar[i]**2), self.template[j,k], mode='valid' )
                    for l in range(npoly):
                        bvec[l+1] = n.sum( (spec[i]*polyarr[l])/(ivar[i]**2) )
                    for m in range(npoly+1): # Fill in entries of A matrix
                        for n in range(npoly+1):
                            if m == 0 and n == 0: amat[m,n] = n.convolve( (self.template[j,k]*self.template[j,k]), ivar[i]**2, mode='valid' )
                            if m == 0 and n != 0:
                                amat[m,n] = n.convolve( self.template[j,k], polyarr[n-1] / (ivar[i]**2), mode='valid' )
                                amat[n,m] = amat[m,n]
                            if m != 0 and n != 0:
                                amat[m,n] = n.sum( (polyarr[m-1]*polyarr[n-1])/(ivar[i]**2) )
                                amat[n,m] = amat[m,n]
                    for l in range(num_z): # Solve Amat.f = bvec at all pixel-redshifts l and computer chi2 at those redshifts
                        f = n.dot(n.linalg.inv(amat[:,:,l]),bvec[:,l])
# now use f to create model at redshift=l, compute chi^2 for that model, and plug it into zchi2arr array