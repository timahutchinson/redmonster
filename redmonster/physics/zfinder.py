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

    def zchi2(self, specs, ivar):
        zchi2 = n.zeros(shape=(n.append(specs.shape[0],self.templates.shape[:-1]))) # Create chi2 array of shape (# of fibers, template_parameter1, ... , template_parameterN)
        # NEED TO ADD LINE TO CREATE POLYNOMIAL TERMS HERE
        for i in range(specs.shape[0]):
            for j in range(self.template.shape[0]):
                for k in range(self.template.shape[1]):
                    conv_vector = n.convolve( spec[i]/ivar[i], self.template[j,k], mode='valid' )
                    conv_matrix = n.convolve( )