import numpy as n
from os import environ
from os.path import join, exists
from redmonster.datamgr.ssp_prep import SSP_Prep
from redmonster.physics.misc import poly_array
from time import gmtime, strftime

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
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
        num_z = self.templates.shape[-1] - specs.shape[-1] + 1 # Number of pixels to be fitted in redshift
        shape = (specs.shape[0],) + self.templates.shape[:-1] + (num_z,)
        zchi2arr = n.zeros(shape) # Create chi2 array of shape (# of fibers, template_parameter_1, ... , template_parameter_N, # of redshifts)
        polyarr = poly_array(self.npoly, specs.shape[1]) # Compute poly terms, noting that they will stay fixed with the data - assumes data is passed in as shape (nfibers, npix)
        bvec = n.zeros((self.npoly+1,num_z))
        amat = n.zeros((self.npoly+1,self.npoly+1,num_z)) # Amat(z).f=bvec(z)
        model_mat = n.zeros((self.templates.shape[-1],self.npoly+1))
        #model_mat[:,1:] = polyarr
        for i in range(specs.shape[0]): # Loop over fibers
            bvec[1:] = n.sum( specs[i]*polyarr*(ivar[i]**2) )
            for j in range(self.templates.shape[0]): # Loop over template-ages
                for k in range(self.templates.shape[1]): # Loop over template vdisps
                    bvec[0] = n.convolve( specs[i]*(ivar[i]**2), self.templates[j,k], mode='valid' )
                    #model_mat[:,0] = self.templates[j,k]
                    for l in range(self.npoly+1): # Fill in entries of A matrix
                        for m in range(self.npoly+1):
                            if l == 0 and m == 0: amat[l,m] = n.convolve( (self.templates[j,k]*self.templates[j,k]), ivar[i]**2, mode='valid' )
                            elif l == 0 and m != 0:
                                amat[l,m] = n.convolve( self.templates[j,k], polyarr[m-1] * (ivar[i]**2), mode='valid' )
                                amat[m,l] = amat[l,m]
                            elif l != 0 and m != 0:
                                amat[l,m] = n.sum( (polyarr[l-1]*polyarr[m-1])*(ivar[i]**2) )
                    for l in range(num_z): # Solve Amat.f = bvec at all pixel-redshifts l and compute chi2 at those redshifts
                        f = n.dot(n.linalg.inv(amat[:,:,l]),bvec[:,l])
                        #model = n.dot(model_mat, f)
                        #zchi2arr[i,j,k,l] = n.sum( ((specs[i]-model[l:l+specs.shape[-1]])**2)*(ivar[i]**2) )
                        model = f[0]*self.templates[j,k]
                        polymodel = n.zeros(specs.shape[-1])
                        for m in range(self.npoly):
                            polymodel += f[m+1]*polyarr[m]
                        zchi2arr[i,j,k,l] = n.sum( ( ( specs[i] - model[l:l+specs.shape[-1]] - polymodel)**2)*(ivar[i]**2) )
                        #import pdb; pdb.set_trace()
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
        return zchi2arr, amat
# now use f to create model at redshift=l, compute chi^2 for that model, and plug it into zchi2arr array