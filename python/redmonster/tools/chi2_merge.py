# Combine individual fiber chi2 files into plate chi2 files for all plates under a given tag.
# fiber_merge.py must be run for each plate BEFORE this can be run.
#
# Tim Hutchinson, University of Utah, February 2015
# t.hutchinson@utah.edu

#import numpy as n
#from astropy.io import fits
from os import environ, makedirs, getcwd
from os.path import exists, join, basename
import re
from time import gmtime, strftime

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
        mjds = []
        try:
            for x in iglob( join( topdir, run2d, str(plate), run1d,
                                 'redmonster-%s-*.fits' % plate) ):
                if mjds is not basename(x)[16:21]:
                    mjds.append( basename(x)[16:21] )
                else:
                    mjds.append( basename(x)[16:21] )
        except: mjds = None
        for mjd in mjds:
            temps = []
            for x in iglob( join( topdir, run2d, str(plate), run1d,
                                 'chi2arr-*-%s-%s-*.fits' % (plate,mjd)) ):
                m = re.search(r'chi2arr-(\D+)-%s-%s-\d+.fits' % (plate,mjd),
                              basename(x))
                if m.group(1) not in temps: temps.append(m.group(1))
            for temp in temps:
                print 'Merging chi2 files for plate %s, mjd %s, template %s' % \
                        (plate, mjd, temp)
                x = io.Merge_Redmonster(plate, mjd, temp)
                x.merge_chi2()