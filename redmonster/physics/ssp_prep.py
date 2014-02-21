import numpy as numpy
from os import environ
from os.path import join,exists
from astropy.io import fits

class ssp_prep:

    def __init__(self):
        try: self.specdir = environ['IDLSPEC2D_DIR']
        except: self.specdir = None
        if self.specdir: self.ssp_path = join(self.specdir,"%s" % "templates","%s" % "SSPs","%s" % "SSP_Padova_RRLIB_Kroupa_Z0.0190.fits")
        if self.ssp_path and exists(self.ssp_path): hdu = fits.open(ssp_path)

    def rebin(self):
        