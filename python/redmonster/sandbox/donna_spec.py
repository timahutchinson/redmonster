from os import environ
from os.path import exists, join
from astropy.io import fits
import numpy as n
from math import ceil, floor
from redmonster.physics.misc import flux_check
from redmonster.datamgr.io import remove_log, write_to_log

class Spec:

    def __init__(self, filepath):
        hdu = fits.open(filepath)
        self.flux = hdu[0].data
        self.ivar = hdu[1].data
        self.loglambda = hdu[0].header['COEFF0'] + n.arange(hdu[0].header['NAXIS1']) * hdu[0].header['COEFF1']
        self.dof = n.array([self.flux.shape[-1]]*self.flux.shape[0])
        self.ebt0 = hdu[5].data.EBOSS_TARGET0
        self.ebt1 = hdu[5].data.EBOSS_TARGET1
        self.npix = hdu[0].header['NAXIS1']
