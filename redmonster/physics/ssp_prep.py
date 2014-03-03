import numpy as n
from os import environ
from os.path import join,exists
from astropy.io import fits
from redmonster.physics.misc import cen2bound
#from time import gmtime, strftime

class ssp_prep:

    def __init__(self, coeff1 = .0001):
        self.ssp_path = None
        self.specs = None
        self.wave = None
        self.wavebounds = None
        self.npix = None
        self.coeff0 = None
        self.coeff1 = coeff1
        try: self.specdir = environ['IDLSPEC2D_DIR']
        except: self.specdir = None
        if self.specdir: self.ssp_path = join(self.specdir,"%s" % "templates","%s" % "SSPs","%s" % "SSP_Padova_RRLIB_Kroupa_Z0.0190.fits")
        if self.ssp_path and exists(self.ssp_path):
            bwb = self.read_ssp()
            self.unit_conv()
            #print strftime("%Y-%m-%d %H:%M:%S", gmtime())
            self.rebin_ssp(bwb)
            #print strftime("%Y-%m-%d %H:%M:%S", gmtime())

    def read_ssp(self):
        if exists(self.ssp_path):
            hdu = fits.open(self.ssp_path)
            self.specs = n.reshape(hdu[1].data.SPEC,(188,-1))
            self.wave = hdu[1].data.LAMBDA[0]
            self.wavebounds = cen2bound(self.wave)
            logwavebounds = n.log10(self.wavebounds)
            self.coeff0 = logwavebounds[0]
            self.npix = int( n.floor( (logwavebounds[-1]-logwavebounds[0]) / self.coeff1 ) )
            bwb = 10**(self.coeff0 + n.arange(self.npix)*self.coeff1)
            return bwb

    # SSPs come in units of L_sun / Hz -- convert this to erg/s/Angstrom to match BOSS spectra
    def unit_conv(self):
        self.specs = ( (3*10**18)/(self.wave**2) ) * (3.839*10**33) * self.specs

    def rebin_ssp(self, bwb = None):
        if bwb != None:
            binspecs = n.zeros(shape=(self.specs.shape[0],self.npix))
            pixlowbound = 0
            parts = n.zeros(shape=(4,self.npix))
            for i in range(self.npix-1):
                parts[0,i] = max(n.where( (bwb[i+1] - self.wavebounds) > 0)[0])
                parts[1,i] = (bwb[i+1]-self.wavebounds[parts[0,i]]) / (self.wavebounds[parts[0,i]+1]-self.wavebounds[parts[0,i]])
                parts[2,i] = (self.wavebounds[pixlowbound]-bwb[i]) / (self.wavebounds[pixlowbound]-self.wavebounds[pixlowbound-1])
                parts[3,i] = pixlowbound
                pixlowbound = parts[0,i]+1
            for j in range(self.specs.shape[0]):
                for i in range(self.npix-1):
                    binspecs[j,i] = ( sum(self.specs[j,parts[3,i]:parts[0,i]]) + (parts[1,i]*self.specs[j,parts[0,i]]) + (parts[2,i]*self.specs[j,parts[3,i]-1]) ) / ( (parts[0,i]-parts[3,i]) + parts[1,i] + parts[2,i] )
            self.specs = binspecs
            self.wavebounds = bwb