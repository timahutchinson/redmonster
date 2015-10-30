# Prepare Charlie Conroy's SSP templates for use in redmonster
#
# Tim Hutchinson, April 2014
# t.hutchinson@utah.edu
#
# May 2014 - now probably defunct, replaced by redmonster.templates.ssp_ztemplate.py


from os import environ
from os.path import join,exists
from time import gmtime, strftime

import numpy as n
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter1d

from redmonster.physics.misc import cen2bound, bound2cen


# Assumes SSPs are a fits file located in $IDLSPEC2D_DIR/templates/SSPs/

class SSPPrep:
    
    def __init__(self, coeff1 = .0001, velmin=None, velstep=None, nvel=None,
                 ssp_file = 'SSP_Padova_RRLIB_Kroupa_Z0.0190.fits'):
        self.ssp_path = None
        self.hdr = None
        self.specs = None
        self.wave = None
        self.wavebounds = None
        self.npix = None
        self.coeff0 = None
        self.coeff1 = coeff1
        self.velmin = velmin
        self.velstep = float(velstep) if velstep else None
        nvel = nvel
        ssp_file = ssp_file
        try: self.specdir = environ['IDLSPEC2D_DIR']
        except: self.specdir = None
        if self.specdir:
            self.ssp_path = join(self.specdir,"%s" % "templates","%s" % "SSPs",
                                 "%s" % ssp_file)
            if exists(self.ssp_path):
                self.read_ssp()
                self.reduce_num() # THIS IS TEMPORARY JUST DURING TESTING
                self.unit_conv()
                self.rebin_ssp()
                if self.velmin and self.velstep and nvel:
                    self.add_velo_disp(self.velmin, self.velstep, nvel)
            else: print "%s IS NOT A VALID PATH" % self.ssp_path
    
    def read_ssp(self):
        if exists(self.ssp_path):
            hdu = fits.open(self.ssp_path)
            self.hdr = hdu[0].header
            self.specs = n.reshape(hdu[1].data.SPEC,(188,-1))
            self.wave = hdu[1].data.LAMBDA[0]
            self.wavebounds = cen2bound(self.wave)
            logwavebounds = n.log10(self.wavebounds)
            self.coeff0 = logwavebounds[0]
            self.npix = int( n.floor( (logwavebounds[-1]-
                                       logwavebounds[0]) / self.coeff1 ) )
    
    # SSPs come in units of L_sun / Hz -- convert this to erg/s/Angstrom
    # to match BOSS spectra
    def unit_conv(self):
        self.specs = ( (3*10**18)/(self.wave**2) ) * (3.839*10**33) * self.specs
    
    def rebin_ssp(self):
        self.specs = gaussian_filter1d(self.specs,
                                       sigma=(self.specs.shape[1]/ \
                                              float(self.npix)), order=0)
           # Binned wavebounds
        bwb = 10**(self.coeff0 + n.arange(self.npix)*self.coeff1)
        binspecs = n.zeros(shape=(self.specs.shape[0],self.npix))
        pixlowbound = 0
        parts = n.zeros(shape=(4,self.npix))
        for i in range(self.npix-1):
            parts[0,i] = max(n.where( (bwb[i+1] - self.wavebounds) > 0)[0])
            parts[1,i] = (bwb[i+1]-self.wavebounds[parts[0,i]]) / \
                    (self.wavebounds[parts[0,i]+1]-self.wavebounds[parts[0,i]])
            parts[2,i] = (self.wavebounds[pixlowbound]-bwb[i]) / \
                    (self.wavebounds[pixlowbound]-self.wavebounds[pixlowbound-1])
            parts[3,i] = pixlowbound
            pixlowbound = parts[0,i]+1
        for j in range(self.specs.shape[0]):
            for i in range(self.npix-1):
                binspecs[j,i] = ( sum(self.specs[j,parts[3,i]:parts[0,i]]) + (parts[1,i]*self.specs[j,parts[0,i]]) + (parts[2,i]*self.specs[j,parts[3,i]-1]) ) / ( (parts[0,i]-parts[3,i]) + parts[1,i] + parts[2,i] )
        self.specs = binspecs * (10**-31) # Rough normalization so SSP fluxes aren't ~50 orders of magnitude larger than BOSS fluxes
        self.wave = bound2cen(bwb)
    
    def add_velo_disp(self, velmin, velstep, nvel):
        velspecs = n.zeros(shape=(self.specs.shape[0], nvel, self.npix))
        for i in range(nvel): velspecs[:,i,:] = gaussian_filter1d(self.specs, sigma=(velmin + velstep*i)/69., order=0) # Divide by 69 due to d(loglam)/dpix = .0001 corresponding to dv/dpix ~ 69 km/s
        self.specs = velspecs
    
    def reduce_num(self):
        newtemps = n.zeros(shape=(self.specs.shape[0]/4., self.specs.shape[1]))
        for i in range(newtemps.shape[0]):
            newtemps[i] = self.specs[4*i]
        self.specs = newtemps