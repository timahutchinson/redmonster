import numpy as n
from os import environ
from os.path import join, exists
from redmonster.datamgr.ssp_prep import SSP_Prep
from redmonster.physics.misc import poly_array
from time import gmtime, strftime

import matplotlib as m
from matplotlib import pyplot as p

class zfinder:

    def __init__(self, config=None, npoly=None):
        self.config = config
        try: self.specdir = environ['IDLSPEC2D_DIR']
        except: self.specdir = None
        if self.config.lower() == 'ssp': self.set_SSP(npoly=npoly)

    def set_SSP(self, npoly=None):
        self.npoly = npoly if npoly else 3
        ssp_stuff = SSP_Prep(velmin=100, velstep=100, nvel=3) # THIS MAY NOT BE THE BEST WAY TO DO THIS
        self.templates, self.wave = ssp_stuff.specs, ssp_stuff.wave

    def zchi2(self, specs, ivar, poffset=0, pspace=0, pmin=0, pmax=0):
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

    def zchi3(self, specs, ivar, poffset=0, pspace=0, pmin=0, pmax=0):
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
        num_z = self.templates.shape[-1] - specs.shape[-1] + 1 # Number of pixels to be fitted in redshift
        shape = (specs.shape[0],) + self.templates.shape[:-1] + (num_z,)
        
        # Create arrays for use in routine
        zchi2arr = n.zeros(shape) # Create chi2 array of shape (# of fibers, template_parameter_1, ... , template_parameter_N, # of redshifts)
        polyarr = poly_array(self.npoly, specs.shape[1]) # Compute poly terms, noting that they will stay fixed with the data - assumes data is passed in as shape (nfibers, npix)
        bvec = n.zeros((self.npoly+1,num_z))
        amat = n.zeros((self.npoly+1,self.npoly+1,num_z)) # Amat(z).f=bvec(z)

        # Pre-compute some matrix elements
        for i in range(self.npoly):
            for j in range(self.npoly):
                amat[i+1,j+1] = n.sum(polyarr[i]*polyarr[j])

    # Do z computation for all fibers
        for i in range(specs.shape[0]): # Loop over fibers
            bvec[1:] = n.sum( specs[i]*polyarr*(ivar[i]**2) )
            amat[1:,1:] = amat[1:,1:] * n.sum(ivar[i]**2)
            for j in range(self.templates.shape[0]): # Loop over template_parameter1
                for k in range(self.templates.shape[1]): # Loop over template_parameter2
                    bvec[0] = n.convolve(specs[i]*(ivar[i]**2), self.templates[j,k], mode='valid')
                    amat[0,0] = n.convolve( self.templates[j,k]**2, ivar[i]**2, mode='valid' )
                    amat[1,0] = n.convolve( self.templates[j,k], polyarr[0]*(ivar[i]**2), mode='valid' )
                    amat[2,0] = n.convolve( self.templates[j,k], polyarr[1]*(ivar[i]**2), mode='valid' )
                    amat[3,0] = n.convolve( self.templates[j,k], polyarr[2]*(ivar[i]**2), mode='valid' )
                    amat[0,1], amat[0,2], amat[0,3] = amat[1,0], amat[2,0], amat[3,0]
                    for l in range(num_z): # Solve Amat.f = bvec at all pixel-redshifts l and compute chi2 at those redshifts
                        f = n.dot(n.linalg.inv(amat[:,:,l]),bvec[:,l])
                        model = f[0]*self.templates[j,k]
                        polymodel = n.zeros(specs.shape[01])
                        for m in range(self.npoly):
                            polymodel += f[m+1]*polyarr[m]
                        zchi2arr[i,j,k,l] = n.sum((( specs[i] - (model[l:l+specs.shape[-1]] + polymodel) )**2)*(ivar[i]**2))
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
        
        
        i = 0; j = 32; k = 2; l = 1787
        bvec[1:] = n.sum( specs[i]*polyarr*(ivar[i]**2) )
        amat[1:,1:] = amat[1:,1:] * n.sum(ivar[i]**2)
        bvec[0] = n.convolve(specs[i]*(ivar[i]**2), self.templates[j,k], mode='valid')
        amat[0,0] = n.convolve( self.templates[j,k]**2, ivar[i]**2, mode='valid' )
        amat[1,0] = n.convolve( self.templates[j,k], polyarr[0]*(ivar[i]**2), mode='valid' )
        amat[2,0] = n.convolve( self.templates[j,k], polyarr[1]*(ivar[i]**2), mode='valid' )
        amat[3,0] = n.convolve( self.templates[j,k], polyarr[2]*(ivar[i]**2), mode='valid' )
        amat[0,1], amat[0,2], amat[0,3] = amat[1,0], amat[2,0], amat[3,0]
        f = n.dot(n.linalg.inv(amat[:,:,l]),bvec[:,l])
        model = f[0]*self.templates[j,k]
        polymodel = n.zeros(specs.shape[01])
        for m in range(self.npoly):
            polymodel += f[m+1]*polyarr[m]
        p.plot(specs[i])
        p.plot(model[l:l+specs.shape[-1]]+polymodel)
        p.show()



        return zchi2arr, amat