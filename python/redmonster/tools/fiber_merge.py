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

try:
    topdir = environ['REDMONSTER_SPECTRO_REDUX']
except KeyError as e:
    topdir = None
    print("Environmental variable 'REDMONSTER_SPECTRO_REDUX' is not set: %r" % e)
try:
    rmver = environ['REDMONSTER_VER']
except KeyError as e:
    rmver = None
    print "Environmental variable 'REDMONSTER_VER' is not set: %r" % e
try:
    run2d = environ['RUN2D']
except KeyError as e:
    run2d = None
    print("Environmental variable 'RUN2D' is not set: %r" % e)
try:
    run1d = environ['RUN1D']
except KeyError as e:
    run1d = None
    print("Environmental variable 'RUN1D' is not set: %r" % e)

platedir = join( topdir, run2d, rmver, '*') if topdir and run2d and rmver else None

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
            for x in iglob( join( topdir, run2d, rmver, str(plate),
                                 'redmonster-%s-*-*.fits' % plate) ):
                files.append(x)
            for file in files:
                '''
                if mjds is not basename(file)[16:21]:
                    mjds.append( basename(file)[16:21] )
                else:
                    mjds.append( basename(file)[16:21] )
                    '''
                if len(plate) == 4:
                    if basename(file)[16:21] not in mjds:
                        mjds.append( basename(file)[16:21] )
                elif len(plate) == 5:
                    if basename(file)[17:22] not in mjds:
                        mjds.append( basename(file)[17:22] )
        except Exception as e:
            mjds = None
            print("Exception: %r" % e)
        if mjds is not []:
            for mjd in mjds:
                print('Merging fibers for plate %s, mjd %s' % (plate, mjd))
                x = io2.MergeRedmonster(plate, mjd)
                x.merge_fibers2()
