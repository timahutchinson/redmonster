import numpy as n
from os import environ
from os.path import join, exists
from redmonster.datamgr.ssp_prep import SSP_Prep
from redmonster.physics.misc import poly_array
from time import gmtime, strftime

import matplotlib as m
from matplotlib import pyplot as p

class Zfinder:
    
    def __init__(self, config=None, npoly=None):
        self.config = config
        try: self.specdir = environ['IDLSPEC2D_DIR']
        except: self.specdir = None
        if self.config.lower() == 'ssp': self.set_SSP(npoly=npoly)
    
    def set_SSP(self, npoly=None):
        self.npoly = npoly if npoly else 3
        ssp_stuff = SSP_Prep(velmin=100, velstep=100, nvel=3) # THIS MAY NOT BE THE BEST WAY TO DO THIS
        self.templates, self.tempwave, self.coeff1 = ssp_stuff.specs, ssp_stuff.wave, ssp_stuff.coeff1
    
    def zchi3(self, specs, ivar, poffset=0, pspace=0, pmin=0, pmax=0):
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
        num_z = self.templates.shape[-1] - specs.shape[-1] + 1 # Number of pixels to be fitted in redshift
        shape = (specs.shape[0],) + self.templates.shape[:-1] + (num_z,)
        fftnaxis1 = 2**13 if self.templates.shape[-1] <= 2**13 else 2**14

        # Create arrays for use in routine
        zchi2arr = n.zeros(shape) # Create chi2 array of shape (# of fibers, template_parameter_1, ... , template_parameter_N, # of redshifts)
        polyarr = poly_array(self.npoly, specs.shape[1]) # Compute poly terms, noting that they will stay fixed with the data - assumes data is passed in as shape (nfibers, npix)
        pmat = n.zeros( (self.npoly+1, self.npoly+1, fftnaxis1), dtype=float)
        bvec = n.zeros( (self.npoly+1, fftnaxis1), dtype=float )
        
        # Pad data and SSPs to a power of 2 for faster FFTs
        temp_pad = n.zeros(self.templates.shape[:-1] + (fftnaxis1,), dtype=float)
        temp_pad[...,:self.templates.shape[-1]] = self.templates
        data_pad = n.zeros(specs.shape[:-1] + (fftnaxis1,), dtype=float)
        data_pad[...,:specs.shape[-1]] = specs
        ivar_pad = n.zeros(ivar.shape[:-1] + (fftnaxis1,), dtype=float)
        ivar_pad[...,:specs.shape[-1]] = ivar
        poly_pad = n.zeros((self.npoly,fftnaxis1), dtype=float)
        poly_pad[...,:polyarr.shape[-1]] = polyarr
        
        # Pre-compute FFTs for use in convolutions
        t_fft = n.fft.fft(temp_pad)
        t2_fft = n.fft.fft(temp_pad**2)
        data_fft = n.fft.fft(data_pad * ivar_pad)
        ivar_fft = n.fft.fft(ivar_pad)
        poly_fft = n.zeros((ivar_pad.shape[0], self.npoly, fftnaxis1), dtype=complex)
        for i in xrange(self.npoly):
            poly_fft[:,i,:] = n.fft.fft(poly_pad[i] * ivar_pad)
        
        # Compute z for all fibers
        for i in xrange(specs.shape[0]): # Loop over fibers
            bvec[1:] = n.sum( poly_pad * data_pad[i] * ivar_pad[i])
            sn_data = n.sum( (specs[i]**2)*ivar[i] )
            for ipos in xrange(self.npoly):
                for jpos in xrange(self.npoly): pmat[ipos+1,jpos+1] = n.sum( poly_pad[ipos] * poly_pad[jpos] * ivar_pad[i])
            for j in xrange(self.templates.shape[0]): # Loop over template_param1
                for k in xrange(self.templates.shape[1]): # Loop over template_param2
                    pmat[0,0] = n.fft.ifft(t2_fft[j,k] * ivar_fft[i].conj()).real
                    bvec[0] = n.fft.ifft(t_fft[j,k] * data_fft[i].conj()).real
                    for ipos in xrange(self.npoly): pmat[ipos+1,0] = pmat[0,ipos+1] = n.fft.ifft(t_fft[j,k] * poly_fft[i,ipos].conj()).real
                    for l in range(num_z): # Solve pmat.f = bvec at all pixel-redshifts l and compute chi2 at those redshifts
                        f = n.linalg.solve(pmat[:,:,l],bvec[:,l])
                        zchi2arr[i,j,k,l] = sn_data - n.dot(n.dot(f,pmat[:,:,l]),f)
                        print sn_data - n.dot(n.dot(f,pmat[:,:,l]),f)
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
        bestl = n.where(zchi2arr == n.min(zchi2arr))[3][0]
        thisz = ((10**(specloglam[0]))/self.tempwave[bestl])-1
        print thisz
        return zchi2arr
        
    
    def zchi4(self, specs, specloglam, ivar):
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
        num_z = self.templates.shape[-1] - specs.shape[-1] + 1 # Number of pixels to be fitted in redshift
        shape = (specs.shape[0],) + self.templates.shape[:-1] + (num_z,)
        self.poff = n.around( (specloglam[0] - n.log10(self.tempwave[0]))/self.coeff1 )
        
        # Create arrays for use in routine
        zchi2arr = n.zeros(shape) # Create chi2 array of shape (# of fibers, template_parameter_1, ... , template_parameter_N, # of redshifts)
        polyarr = poly_array(self.npoly, specs.shape[1]) # Compute poly terms, noting that they will stay fixed with the data - assumes data is passed in as shape (nfibers, npix)
        bvec = n.zeros((self.npoly+1,num_z))
        pmat = n.zeros((self.npoly+1,self.npoly+1,num_z))
        
        # Pre-compute some matrix elements
        #for i in range(self.npoly):
        #    for j in range(self.npoly):
        #        pmat[i+1, j+1] = n.sum(polyarr[i]*polyarr[j])
        
        # Compute z for all fibers
        for i in range(specs.shape[0]): # Loop over fibers
            bvec[1:] = n.sum(specs[i]*polyarr*ivar[i])
            sn_data = n.sum( (specs[i]**2)*ivar[i] )
            for irow in range(self.npoly):
                for icol in range(self.npoly):
                    pmat[irow+1,icol+1] = n.sum(polyarr[irow]*polyarr[icol]*ivar[i])
            for j in range(self.templates.shape[0]): # Loop over template_param1
                for k in range(self.templates.shape[1]): # Loop over template_param2
                    bvec[0] = n.convolve(specs[i]*ivar[i], self.templates[j,k], mode='valid')
                    pmat[0,0] = n.convolve( self.templates[j,k]**2, ivar[i], mode='valid')
                    pmat[1,0] = n.convolve( self.templates[j,k], polyarr[0]*ivar[i], mode='valid')
                    pmat[2,0] = n.convolve( self.templates[j,k], polyarr[1]*ivar[i], mode='valid')
                    pmat[3,0] = n.convolve( self.templates[j,k], polyarr[2]*ivar[i], mode='valid')
                    pmat[0,1], pmat[0,2], pmat[0,3] = pmat[1,0], pmat[2,0], pmat[3,0]
                    for l in range(num_z): # Solve Pmat.f = bvec at all pixel-redshifts l and compute chi2 at those redshifts
                        #f = n.dot(n.linalg.inv(pmat[:,:,l]),bvec[:,l])
                        f = n.linalg.solve(pmat[:,:,l],bvec[:,l])
                        model = f[0]*self.templates[j,k]
                        polymodel = n.zeros(specs.shape[-1])
                        for m in range(self.npoly): polymodel += f[m+1]*polyarr[m]
                        #coordtrans = n.zeros((4,))
                        #zchi2arr[i,j,k,l] = n.sum((( specs[i] - (model[l:l+specs.shape[-1]] + polymodel) )**2)*(ivar[i]))
                        #zchi2arr[i,j,k,l] = n.sum( (specs[i]**2)*ivar[i]) - n.dot(n.dot(n.transpose(f),pmat[:,:,l]),f)
                        print zchi2arr[i,j,k,l]
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
        bestl = n.where(zchi2arr == n.min(zchi2arr))[3][0]
        thisz = ((10**(specloglam[0]))/self.tempwave[bestl])-1
        print thisz
        
        
        
        return zchi2arr






