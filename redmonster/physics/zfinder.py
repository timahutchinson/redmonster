import numpy as n
from os import environ
from os.path import join

class zfinder:

    def __init__(self, specobj=None, npoly=None):
        self.npoly = npoly
        try: self.specdir = environ['IDLSPEC2D_DIR']
        except: self.specdir = None
        if self.specdir: self.eigendir = join(self.specdir,"%s" % "templates")