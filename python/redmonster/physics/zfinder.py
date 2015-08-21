# First pass redshift finder of redmonster.  Able to search entire parameter-space and redshift-space
# of models in increments of d(loglam)/d(pixel)
#
# Tim Hutchinson, University of Utah @ IAC, May 2014
# t.hutchinson@utah.edu

import numpy as n
from os import environ, makedirs
from os.path import join, exists
from redmonster.datamgr.ssp_prep import SSP_Prep
from redmonster.physics.misc import poly_array, two_pad
from time import gmtime, strftime
from astropy.io import fits
from redmonster.datamgr.io2 import read_ndArch, write_chi2arr
from scipy.optimize import nnls

import matplotlib as m
from matplotlib import pyplot as p

# Assumes all templates live in $REDMONSTER_DIR/templates/

class Zfinder:
    
    def __init__(self, fname=None, npoly=None, zmin=None, zmax=None):
        self.fname = fname
        self.npoly = npoly if npoly else 4
        self.zmin = float(zmin)
        self.zmax = float(zmax)
        self.pixoffset = None
        self.zchi2arr = None
        try: self.templatesdir = environ['REDMONSTER_TEMPLATES_DIR']
        except: self.templatesdir = None
        self.read_template()
        self.npars = len(self.templates.shape) - 1
        self.templates_flat = n.reshape(self.templates, (-1,self.fftnaxis1))
        self.t_fft = n.fft.fft(self.templates_flat)
        self.t2_fft = n.fft.fft(self.templates_flat**2)
    

    def read_template(self):
        self.templates, self.baselines, self.infodict = read_ndArch(join(self.templatesdir,self.fname))
        self.type = self.infodict['class']
        self.origshape = self.templates.shape
        self.ntemps = self.templates[...,0].size
        self.fftnaxis1 = two_pad(self.origshape[-1])
        self.tempwave = 10**(self.infodict['coeff0'] + n.arange(self.infodict['nwave'])*self.infodict['coeff1'])
        templates_pad = n.zeros( self.origshape[:-1]+(self.fftnaxis1,) )
        templates_pad[...,:self.origshape[-1]] = self.templates
        self.templates = templates_pad
    
    
    def create_z_baseline(self, loglam0):
        self.zbase = ((10**loglam0)/self.tempwave) - 1
    
    
    def conv_zbounds(self):
        zmaxpix = n.where( abs((self.zbase-self.zmin)) == n.min(abs(self.zbase-self.zmin)) )[0][0]
        zminpix = n.where( abs((self.zbase-self.zmax)) == n.min(abs(self.zbase-self.zmax)) )[0][0]
        return zminpix, zmaxpix
    
    
    def zchi2(self, specs, specloglam, ivar, npixstep=1, chi2file=False, plate=None, mjd=None, fiberid=None):
        self.chi2file = chi2file
        self.npixstep = npixstep
        self.zwarning = n.zeros(specs.shape[0])
        flag_val_unplugged = int('0b10000000',2)
        flag_val_neg_model = int('0b1000',2)
        #print strftime("%Y-%m-%d %H:%M:%S", gmtime()) # For timing while testing
        self.create_z_baseline(specloglam[0])
        if (self.zmin != None) and (self.zmax != None) and (self.zmax > self.zmin):
            bounds_set = True
            zminpix, zmaxpix = self.conv_zbounds()
            self.pixoffset = zminpix
            num_z = int(n.floor( (zmaxpix - zminpix) / npixstep )) #num_z = zmaxpix - zminpix + 1 # Number of pixels to be fitted in redshift
            zinds = zminpix + n.arange(num_z)*npixstep
            self.zbase = self.zbase[zinds] #self.zbase = self.zbase[zminpix:zminpix+num_z]
        else:
            bounds_set = False
            num_z = int(n.floor( (zself.origshape[-1] - specs.shape[-1]) / npixstep )) #num_z = self.origshape[-1] - specs.shape[-1] + 1 # Number of pixels to be fitted in redshift
        
        # Create arrays for use in routine
        zchi2arr = n.zeros((specs.shape[0], self.templates_flat.shape[0], num_z)) # Create chi2 array of shape (# of fibers, template_parameter_1, ... , template_parameter_N, # of redshifts)
        temp_zwarning = n.zeros(zchi2arr.shape)
        polyarr = poly_array(self.npoly, specs.shape[1]) # Compute poly terms, noting that they will stay fixed with the data - assumes data is passed in as shape (nfibers, npix)
        pmat = n.zeros( (self.npoly+1, self.npoly+1, self.fftnaxis1), dtype=float)
        bvec = n.zeros( (self.npoly+1, self.fftnaxis1), dtype=float)
        
        # Pad data and SSPs to a power of 2 for faster FFTs
        data_pad = n.zeros(specs.shape[:-1] + (self.fftnaxis1,), dtype=float)
        data_pad[...,:specs.shape[-1]] = specs
        ivar_pad = n.zeros(ivar.shape[:-1] + (self.fftnaxis1,), dtype=float)
        ivar_pad[...,:specs.shape[-1]] = ivar
        poly_pad = n.zeros((self.npoly, self.fftnaxis1), dtype=float)
        poly_pad[...,:polyarr.shape[-1]] = polyarr
        
        # Pre-compute FFTs for use in convolutions
        data_fft = n.fft.fft(data_pad * ivar_pad)
        ivar_fft = n.fft.fft(ivar_pad)
        poly_fft = n.zeros((ivar_pad.shape[0], self.npoly, self.fftnaxis1), dtype=complex)
        for i in xrange(self.npoly):
            poly_fft[:,i,:] = n.fft.fft(poly_pad[i] * ivar_pad)
        
        # Compute z for all fibers
        
        for i in xrange(specs.shape[0]): # Loop over fibers
            print 'Fitting fiber %s of %s for template %s' % (i+1, specs.shape[0], self.fname) #print i, self.fname
            #write_to_log(
            if len(n.where(specs[i] != 0.)[0]) == 0: # If flux is all zeros, flag as unplugged according to BOSS zwarning flags and don't bother with doing fit
                self.zwarning[i] = int(self.zwarning[i]) | flag_val_unplugged
            else: # Otherwise, go ahead and do fit
                for ipos in xrange(self.npoly): bvec[ipos+1] = n.sum( poly_pad[ipos] * data_pad[i] * ivar_pad[i])
                sn2_data = n.sum( (specs[i]**2)*ivar[i] )
                for ipos in xrange(self.npoly):
                    for jpos in xrange(self.npoly): pmat[ipos+1,jpos+1] = n.sum( poly_pad[ipos] * poly_pad[jpos] * ivar_pad[i])
                for j in xrange(self.templates_flat.shape[0]): # Loop over templates
                    pmat[0,0] = n.fft.ifft(self.t2_fft[j] * ivar_fft[i].conj()).real
                    bvec[0] = n.fft.ifft(self.t_fft[j] * data_fft[i].conj()).real
                    for ipos in xrange(self.npoly): pmat[ipos+1,0] = pmat[0,ipos+1] = n.fft.ifft(self.t_fft[j] * poly_fft[i,ipos].conj()).real
                    if bounds_set:
                        for l in n.arange(num_z)*self.npixstep: #for l in range(num_z):
                            try:
                                f = n.linalg.solve(pmat[:,:,l+zminpix],bvec[:,l+zminpix])
                                zchi2arr[i,j,(l/self.npixstep)] = sn2_data - n.dot(n.dot(f,pmat[:,:,l+zminpix]),f)
                                if (f[0] < 0):
                                    temp_zwarning[i,j,(l/self.npixstep)] = int(temp_zwarning[i,j,(l/self.npixstep)]) | flag_val_neg_model #if (f[0] < 0): temp_zwarning[i,j,l] = int(temp_zwarning[i,j,l]) | flag_val_neg_model
                                    zchi2arr[i,j,(l/self.npixstep)] = 1e9
                            except Exception as e:
                                print "Exception: %r" % e
                                zchi2arr[i,j,(l/self.npixstep)] = 1e9
                    else:
                        for l in n.arange(num_z)*self.npixstep: #for l in range(num_z):
                            try:
                                f = n.linalg.solve(pmat[:,:,l],bvec[:,l])
                                zchi2arr[i,j,(l/self.npixstep)] = sn2_data - n.dot(n.dot(f,pmat[:,:,l]),f)
                                if (f[0] < 0.):
                                    temp_zwarning[i,j,(l/self.npixstep)] = int(temp_zwarning[i,j,(l/self.npixstep)]) | flag_val_neg_model #if (f[0] < 0.): temp_zwarning[i,j,l] = int(temp_zwarning[i,j,l]) | flag_val_neg_model
                                    zchi2arr[i,j,(l/self.npixstep)] = 1e9
                            except Exception as e:
                                print "Exception: %r" % e
                                zchi2arr[i,j,(l/self.npixstep)] = 1e9
        for i in xrange(self.zwarning.shape[0]): # Use only neg_model flag from best fit model/redshift and add it to self.zwarning
            minpos = ( n.where(zchi2arr[i] == n.min(zchi2arr[i]))[0][0] , n.where(zchi2arr[i] == n.min(zchi2arr[i]))[1][0] )
            self.zwarning[i] = int(self.zwarning[i]) | int(temp_zwarning[i,minpos[0],minpos[1]])
        zchi2arr = n.reshape(zchi2arr, (specs.shape[0],) + self.origshape[:-1] + (num_z,) )
        bestl = n.where(zchi2arr == n.min(zchi2arr))[-1][0]
        if bounds_set:
            thisz = ((10**(specloglam[0]))/self.tempwave[bestl+zminpix])-1
        else:
            thisz = ((10**(specloglam[0]))/self.tempwave[bestl])-1
        #return zchi2arr
        self.zchi2arr = zchi2arr
        #self.store_models(specs, ivar)
        if self.chi2file is True:
            if (plate is not None) & (mjd is not None) & (fiberid is not None):
                write_chi2arr(plate, mjd, fiberid, self.zchi2arr)
            else:
                print 'Plate/mjd/fiberid not given - unable to write chi2 file!'
        else:
            print 'Not writing chi2'


    def store_models(self, specs, ivar):
        self.models = n.zeros( (specs.shape) )
        for i in xrange(self.models.shape[0]):
            minloc = n.unravel_index( self.zchi2arr[i].argmin(), self.zchi2arr[i].shape )
            pmat = n.zeros( (specs.shape[-1],self.npoly+1) )
            this_temp = self.templates[minloc[:-1]]
            pmat[:,0] = this_temp[(minloc[-1]*self.npixstep)+self.pixoffset:(minloc[-1]*self.npixstep)+self.pixoffset+specs.shape[-1]]
            polyarr = poly_array(self.npoly, specs.shape[-1])
            pmat[:,1:] = n.transpose(polyarr)
            ninv = n.diag(ivar[i])
            try: # Some eBOSS spectra have ivar[i] = 0 for all i
                '''
                f = n.linalg.solve( n.dot(n.dot(n.transpose(pmat),ninv),pmat), n.dot( n.dot(n.transpose(pmat),ninv),specs[i]) )
                '''
                f = nnls( n.dot(n.dot(n.transpose(pmat),ninv),pmat), n.dot( n.dot(n.transpose(pmat),ninv),specs[i]) ); f = n.array(f)[0]
                self.models[i] = n.dot(pmat, f)
            except:
                self.models[i] = n.zeros(specs.shape[-1])


    def write_chi2arr(self, plate, mjd, fiberid):
        prihdu = fits.PrimaryHDU(self.zchi2arr)
        thdulist = fits.HDUList([prihdu])
        try:
            rsr = environ['REDMONSTER_SPECTRO_REDUX']
            run2d = environ['RUN2D']
            run1d = environ['RUN1D']
            if (rsr is not None) & (run2d is not None) & (run1d is not None):
                testpath = join(rsr, run2d, '%s' % plate, run1d)
                if exists(testpath):
                    dest = testpath
                else:
                    try:
                        makedirs(testpath)
                        dest = testpath
                    except:
                        dest = None
        except:
            dest = None
        if dest is not None:
            try:
                thdulist.writeto(join(dest, '%s' % 'chi2arr-%s-%s-%s-%03d.fits' % (self.type, plate, mjd, fiberid)), clobber=True)
                print 'Writing chi2 file to %s' % join(dest, '%s' % 'chi2arr-%s-%s-%s-%03d.fits' % (self.type, plate, mjd, fiberid))
            except:
                print 'Environment variables not set or path does not exist - not writing chi2 file!'
        else:
            print 'Environment variables not set or path does not exist - not writing chi2 file!'





























