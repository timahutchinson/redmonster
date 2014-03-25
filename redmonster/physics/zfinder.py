import numpy as n
from os import environ
from os.path import join, exists
from redmonster.datamgr.ssp_prep import SSP_Prep

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
        #self.eigendir = join(self.specdir,"%s" % "templates", "%s" % "SSPs") if self.specdir else None

    def zchi2(self, specs, ivar):
        #ndim = len(self.templates.shape) - 1
        zchi2 = n.append(n.zeros(specs.shape[0]),n.zeros(shape=self.templates.shape[:-1]))
        return zchi2