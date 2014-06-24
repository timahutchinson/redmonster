#
# redmonster/datamgr/sdss.py
#
# SDSS-specific data interface tools
#

import numpy as n
import matplotlib as m
m.interactive(True)
from matplotlib import pyplot as p
from astropy.io import fits
import os
from redmonster.math import misc


class SpCFrameAll:
    """
    Simultaneous interface to all spCFrame files for a given plate/mjd
    """
    def __init__(self, plate, mjd, topdir=None, run2d=None):
        if topdir is None:
            topdir = os.getenv('BOSS_SPECTRO_REDUX')
        if run2d is None:
            run2d = os.getenv('RUN2D')
        pstring = str(int(plate)).strip()
        mstring = str(int(mjd)).strip()
        self.path_to_spectra = topdir + '/' + run2d + '/' + pstring + '/'
        self.spPlate_file = self.path_to_spectra + 'spPlate-' + pstring \
                            + '-' + mstring + '.fits'
        hdr = fits.getheader(self.spPlate_file)
        exp_keys = [this_key for this_key in hdr.keys() if this_key[:5] == 'EXPID']
        exp_ids = [hdr[this_key][:11] for this_key in exp_keys]
        self.spCFrame_files = [self.path_to_spectra + 'spCFrame-' + this_id +
                               '.fits' for this_id in exp_ids]
        self.flux_list = [fits.getdata(this_file) for this_file in self.spCFrame_files]
        self.invvar_list = [fits.getdata(this_file, 1) for this_file in self.spCFrame_files]
        self.mask_list = [fits.getdata(this_file, 2) for this_file in self.spCFrame_files]
        self.loglam_list = [fits.getdata(this_file, 3) for this_file in self.spCFrame_files]
        self.disp_list = [fits.getdata(this_file, 4) for this_file in self.spCFrame_files]
        self.plug_list = [fits.getdata(this_file, 5) for this_file in self.spCFrame_files]
        self.sky_list = [fits.getdata(this_file, 6) for this_file in self.spCFrame_files]
        self.xpos_list = [fits.getdata(this_file, 7) for this_file in self.spCFrame_files]
        self.superflat_list = [fits.getdata(this_file, 8) for this_file in self.spCFrame_files]
        self.n_exposures = len(self.spCFrame_files)
    def set_fiber(self, fiberid):
        self.fiberid_now = fiberid
        self.i_exp = []
        self.j_row = []
        for i in xrange(self.n_exposures):
            wh_fib = n.where(self.plug_list[i].FIBERID == fiberid)[0]
            for j in xrange(len(wh_fib)):
                self.i_exp.append(i)
                self.j_row.append(wh_fib[j])
        self.nspec_fib = len(self.i_exp)
        self.flux_fib = [self.flux_list[self.i_exp[k]][self.j_row[k]]
                         for k in xrange(self.nspec_fib)]
        self.invvar_fib = [self.invvar_list[self.i_exp[k]][self.j_row[k]]
                           for k in xrange(self.nspec_fib)]
        self.mask_fib = [self.mask_list[self.i_exp[k]][self.j_row[k]]
                         for k in xrange(self.nspec_fib)]
        self.loglam_fib = [self.loglam_list[self.i_exp[k]][self.j_row[k]]
                           for k in xrange(self.nspec_fib)]
        self.disp_fib = [self.disp_list[self.i_exp[k]][self.j_row[k]]
                         for k in xrange(self.nspec_fib)]

