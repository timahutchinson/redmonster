# Combine individual fiber files into plate files for all eBOSS plates
#
# Tim Hutchinson, University of Utah, December 2014
# t.hutchinson@utah.edu

#import numpy as n
#from astropy.io import fits
from os import environ, makedirs, getcwd
from os.path import exists, join, basename
# time import gmtime, strftime
from glob import iglob
from redmonster.datamgr import io

try: topdir = environ['REDMONSTER_SPECTRO_REDUX']
except: topdir = None
try: run2d = environ['RUN2D']
except: run2d = None
try: run1d = environ['RUN1D']
except: run1d = None
platedir = join( topdir, run2d, '*') if topdir and run2d else None

plates = []
    
if platedir:
    for path in iglob(platedir):
        plates.append( basename(path) )
        plates.sort()
    for plate in plates:
        mjds = []
        try:
            for x in iglob( join( topdir, run2d, str(plate), run1d, 'redmonster-%s-*-000.fits' % plate) ):
                if mjds is not basename(x)[16:21]:
                    mjds.append( basename(x)[16:21] )
                else:
                    mjds.append( basename(x)[16:21] )
        except: mjds = None
        for mjd in mjds:
            x = io.Merge_Redmonster(plate, mjd)
            x.merge_fibers()