import os
from astropy.io import fits
import numpy as n

class Spec:

    def __init__(self, plate=None, mjd=None):
        self.set_plate_mjd(plate=plate, mjd=mjd)
        self.set_platepath()
        global hdu
        hdu = fits.open(self.platepath)
        self.set_flux()
        self.set_ivar()
        self.set_loglambda()
    
    def set_plate_mjd(self, plate=None, mjd=None):
        self.plate = plate
        self.mjd = mjd
    
    def set_platepath(self):
        self.platepath = os.environ['BOSS_SPECTRO_REDUX'] + os.environ['RUN2D'] + str(self.plate) + '/' + 'spPlate-' + str(self.plate) + '-' + str(self.mjd) + '.fits'
    
    def set_flux(self):
        self.flux = hdu[0].data
    
    def set_ivar(self):
        self.ivar = hdu[1].data
        
    def set_loglambda(self):
        self.loglambda = hdu[0].header['COEFF0'] + n.arange(hdu[0].header['NAXIS1']) * hdu[0].header['COEFF1']

    def set_andmask(self):
        self.andmask = hdu[2].data

    def set_ormask(self):
        self.ormask = hdu[3].data

    def set_plugmap(self):
        self.plugmap = hdu[5].data

    def set_skyflux(self):
        self.skyflux = hdu[6].data