from os import environ
from os.path import exists, join
from astropy.io import fits
import numpy as n

class Spec:

    def __init__(self, plate=None, mjd=None):
        self.flux = None
        self.ivar = None
        self.loglambda = None
        self.ormask = None
        self.plugmap = None
        self.skyflux = None
        try: self.topdir = environ['BOSS_SPECTRO_REDUX']
        except: self.topdir = None
        try: self.run2d = environ['RUN2D']
        except: self.run2d = None
        self.set_plate_mjd(plate=plate, mjd=mjd)
    
    def set_plate_mjd(self, plate=None, mjd=None):
        self.plate = plate
        self.mjd = mjd
        if self.topdir and self.run2d and self.plate and self.mjd:
            self.platepath = join(self.topdir,self.run2d,"%s" % self.plate,"spPlate-%s-%s.fits" % (self.plate,self.mjd))
            self.set_data()
    
    def set_data(self):
        if self.platepath and exists(self.platepath): hdu = fits.open(self.platepath)
        else: print "Missing path to %r" % self.platepath
        try:
            self.flux = hdu[0].data
            self.ivar = hdu[1].data
            self.loglambda = hdu[0].header['COEFF0'] + n.arange(hdu[0].header['NAXIS1']) * hdu[0].header['COEFF1']
            self.andmask = hdu[2].data
            self.ormask = hdu[3].data
            self.plugmap = hdu[5].data
            self.skyflux = hdu[6].data
        except Exception as e: print "Exception: %r" % e
