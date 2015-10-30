# Combine individual fiber files into plate files for all eBOSS plates
#
# Tim Hutchinson, University of Utah, December 2014
# t.hutchinson@utah.edu


from os import environ, makedirs, getcwd
from os.path import exists, join, basename
from time import gmtime, strftime

import numpy as n
from astropy.io import fits
from glob import iglob

from redmonster.datamgr import io2

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
        #import pdb; pdb.set_trace()
        mjds = []
        files = []
        try:
            for x in iglob( join( topdir, run2d, str(plate), run1d,
                                 'redmonster-%s-*-*.fits' % plate) ):
                files.append(x)
            for file in files:
                '''
                if mjds is not basename(file)[16:21]:
                    mjds.append( basename(file)[16:21] )
                else:
                    mjds.append( basename(file)[16:21] )
                    '''
                if basename(file)[16:21] not in mjds:
                    mjds.append( basename(file)[16:21] )
        except: mjds = None
        if mjds is not []:
            for mjd in mjds:
                print 'Merging fibers for plate %s, mjd %s' % (plate, mjd)
                x = io2.MergeRedmonster(plate, mjd)
                x.merge_fibers2()