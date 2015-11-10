# Read SDSS spPlate files for use in redmonster
#
# Tim Hutchinson, University of Utah @ IAC, May 2014
# t.hutchinson@utah.edu

from os import environ
from os.path import exists, join
from math import ceil, floor

import numpy as n
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel

from redmonster.physics.misc import flux_check
from redmonster.datamgr.io import remove_log, write_to_log



class Spec:

    def __init__(self, plate=None, mjd=None, fiberid=None, data_range=None):
        print 'plate %s mjd %s fiberid %s' % (plate, mjd, fiberid)
        self.hdr = None
        self.flux = None
        self.ivar = None
        self.loglambda = None
        self.andmask = None
        self.ormask = None
        self.plugmap = None
        self.skyflux = None
        self.nobj = None
        self.npix = None
        self.coeff0 = None
        self.coeff1 = None
        self.dof = None
        try:
            self.topdir = environ['BOSS_SPECTRO_REDUX']
        except KeyError as e:
            self.topdir = None
            print "Enviromental variable 'BOSS_SPECTRO_REDUX' not set: %r" % e
        try:
            self.run2d = environ['RUN2D']
        except KeyError:
            self.run2d = None
            print "Enviromental variable 'RUN2D' not set: %r" % e
        self.set_plate_mjd(plate=plate, mjd=mjd, fiberid=fiberid,
                           data_range=data_range)
        #for i in xrange(self.flux.shape[0]):
            #self.flux[i] = convolve(self.flux[i], Box1DKernel(5))
        if exists(self.platepath):
            self.ivar, self.dof = flux_check(self.flux, self.ivar, plate, mjd)
        else:
            # This should probably be logged eventually as well
            print '%s does not exist.' % self.platepath
            # Added by TH 21 July 2015
            #write_to_log(plate, mjd, '%s does not exist.' % self.platepath)

    def set_plate_mjd(self, plate=None, mjd=None, fiberid=None,
                      data_range=None):
        self.plate = plate
        self.mjd = mjd
        if self.topdir and self.run2d and self.plate and self.mjd:
            self.platepath = join(self.topdir, self.run2d, "%s" % self.plate,
                                  "spPlate-%s-%s.fits" % (self.plate,self.mjd))
            if exists(self.platepath):
                self.set_data()
                if data_range: self.chop_data(data_range)
                if fiberid != None: self.set_fibers(fiberid)
                else: self.fiberid = [i for i in xrange(1000)]
                self.flag_sky_fibers()
    
    def set_data(self):
        if self.platepath and exists(self.platepath):
            hdu = fits.open(self.platepath)
        else:
            print "Missing path to %r" % self.platepath
            # Added by TH 21 July 2015
            #write_to_log(self.plate, self.mjd,
                          #'Missing path to %r' % self.platepath)
        try:
            self.hdr = hdu[0].header
            self.flux = hdu[0].data
            self.ivar = hdu[1].data
            self.loglambda = (hdu[0].header['COEFF0'] +
                              n.arange(hdu[0].header['NAXIS1']) *
                              hdu[0].header['COEFF1'])
            self.andmask = hdu[2].data
            self.ormask = hdu[3].data
            self.plugmap = hdu[5].data
            try:
                self.boss_target1 = hdu[5].data.BOSS_TARGET1
            except AttributeError:
                pass
            try:
                self.eboss_target0 = hdu[5].data.EBOSS_TARGET0
            except AttributeError:
                pass
            try:
                self.eboss_target1 = hdu[5].data.EBOSS_TARGET1
            except AttributeError:
                pass
            self.skyflux = hdu[6].data
            # For plate files before Spectro-2D v5, there are no sky vectors and
            # hdu[6] is something else
            if self.skyflux.shape != self.flux.shape: self.skyflux = 0
            self.nobj = hdu[0].header['NAXIS2']
            self.npix = hdu[0].header['NAXIS1']
            self.coeff0 = hdu[0].header['COEFF0']
            self.coeff1 = hdu[0].header['COEFF1']
        except Exception as e:
            print "Exception: %r" % e
            # Added by TH 21 July 2015
            #write_to_log(self.plate, self.mjd, 'Exception: %r' % e')

    def chop_data(self, data_range):
        self.data_range = data_range
        i1 = ceil( (n.log10(data_range[0]) - self.coeff0) / self.coeff1 )
        i2 = floor( (n.log10(data_range[1]) - self.coeff0) / self.coeff1 )
        if i1 >= 0: self.ivar[:,:i1] = 0
        if i2 <= self.npix: self.ivar[:,i2:] = 0
        # CHANGE PRINT STATEMENT TO LOG - done by TH 21 July 2015
        print 'Trim wavelength range to %s' % data_range
        #write_to_log(self.plate, self.mjd,
                     #'Trim wavelength range to %s' % data_range)

    def set_fibers(self, fiberid):
        if min(fiberid) < 0 or max(fiberid) > self.nobj:
            # CHANGE THIS TO LOG INSTEAD OF PRINT - done TH
            print ('Invalid value for FIBERID: must be between 0 and %s' %
                   hdu[0].header['NAXIS1'])
            #write_to_log(self.plate, self.mjd,
                         #'Invalid value for FIBERID: must be between 0 and %s'%
                         #hdu[0].header['NAXIS1'])
        else:
            self.fiberid = fiberid
            self.flux = self.flux[fiberid]
            self.ivar = self.ivar[fiberid]
            self.andmask = self.andmask[fiberid]
            self.ormask = self.ormask[fiberid]
            self.plugmap = self.plugmap[fiberid]
            try:
                self.boss_target1 = self.boss_target1[fiberid]
            except NameError:
                pass
            try:
                self.eboss_target0 = self.eboss_target0[fiberid]
            except NameError:
                pass
            try:
                self.eboss_target1 = self.eboss_target1[fiberid]
            except NameError:
                pass
            self.nobj = len(fiberid)
            if self.skyflux.shape[0] != 1: self.skyflux = self.skyflux[fiberid]

    def flag_sky_fibers(self):
        self.zwarning = n.zeros( len(self.fiberid) )
        flag_val = int('0b1',2) # From BOSS zwarning flag definitions
        for i in xrange(self.plugmap.shape[0]):
            if ( self.plugmap[i]['OBJTYPE'].lower() == 'sky'):
                self.zwarning[i] = int(self.zwarning[i]) | flag_val












