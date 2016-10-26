from __future__ import print_function

from os.path import join, basename, exists
from os import environ
from math import isnan
import time
from sys import stderr

import numpy as n
from scipy.integrate import trapz
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
from matplotlib.colors import LogNorm
from glob import iglob
from astropy.convolution import convolve, Box1DKernel
from scipy.optimize import curve_fit
from scipy import ndimage
import seaborn as sns

from redmonster.sandbox import yanny as y
from redmonster.datamgr import spec
from redmonster.physics import zfinder
from redmonster.datamgr.io import read_ndArch
from redmonster.physics.misc import poly_array

class VerifyRM:

    def __init__(self,version='v5_10_0',
                 plates=[3686,3687,3804,3805,3853,3855,3856,3860],
                 mjds={
                        3686:55268,3687:55269,3804:55267,3805:55269,3853:55268,
                        3855:55268,3856:55269,3860:55269
                 },
                 sns_pal='muted'):
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')
        self.version = version
        self.plates = plates
        self.mjds = mjds
        self.redmonster_spectro_redux = \
                join( environ['REDMONSTER_SPECTRO_REDUX'], '%s' % self.version)
        self.vifibers = None
        self.zperson = None
        self.zpipe = None
        self.vitype = None
        self.comments = None
        #self.yanny_to_arrays()
        #self.rm_z = []
        #self.rm_class = []
        #self.rm_zwarning = []
        #self.vis_z = []

    def yanny_to_arrays(self):
        # Convert yanny file to arrays
        # Read yanny file
        x = y.yanny(filename='/uufs/astro.utah.edu/common/home/u0814744/boss/\
                    spInspect_alltest_bolton.par.txt', np=True)
        # Get fibers, zpipe, zperson for each plate
        args = n.where(x['BOSSOBJECT']['plate'] == 3686)[0]
        self.fibers3686 = []
        self.zpipe3686 = []
        self.zperson3686 = []
        self.type3686 = []
        self.comments3686 = []
        for i in args:
            self.fibers3686.append( x['BOSSOBJECT'][i][2])
            self.zpipe3686.append( x['BOSSOBJECT'][i][5])
            self.zperson3686.append( x['BOSSOBJECT'][i][6])
            self.type3686.append( x['BOSSOBJECT'][i][7])
            self.comments3686.append( x['BOSSOBJECT'][i][8])
        args = n.where(x['BOSSOBJECT']['plate'] == 3687)[0]
        self.fibers3687 = []
        self.zpipe3687 = []
        self.zperson3687 = []
        self.type3687 = []
        self.comments3687 = []
        for i in args:
            self.fibers3687.append( x['BOSSOBJECT'][i][2])
            self.zpipe3687.append( x['BOSSOBJECT'][i][5])
            self.zperson3687.append( x['BOSSOBJECT'][i][6])
            self.type3687.append( x['BOSSOBJECT'][i][7])
            self.comments3687.append( x['BOSSOBJECT'][i][8])
        args = n.where(x['BOSSOBJECT']['plate'] == 3804)[0]
        self.fibers3804 = []
        self.zpipe3804 = []
        self.zperson3804 = []
        self.type3804 = []
        self.comments3804 = []
        for i in args:
            self.fibers3804.append( x['BOSSOBJECT'][i][2])
            self.zpipe3804.append( x['BOSSOBJECT'][i][5])
            self.zperson3804.append( x['BOSSOBJECT'][i][6])
            self.type3804.append( x['BOSSOBJECT'][i][7])
            self.comments3804.append( x['BOSSOBJECT'][i][8])
        args = n.where(x['BOSSOBJECT']['plate'] == 3805)[0]
        self.fibers3805 = []
        self.zpipe3805 = []
        self.zperson3805 = []
        self.type3805 = []
        self.comments3805 = []
        for i in args:
            self.fibers3805.append( x['BOSSOBJECT'][i][2])
            self.zpipe3805.append( x['BOSSOBJECT'][i][5])
            self.zperson3805.append( x['BOSSOBJECT'][i][6])
            self.type3805.append( x['BOSSOBJECT'][i][7])
            self.comments3805.append( x['BOSSOBJECT'][i][8])
        args = n.where(x['BOSSOBJECT']['plate'] == 3853)[0]
        self.fibers3853 = []
        self.zpipe3853 = []
        self.zperson3853 = []
        self.type3853 = []
        self.comments3853 = []
        for i in args:
            self.fibers3853.append( x['BOSSOBJECT'][i][2])
            self.zpipe3853.append( x['BOSSOBJECT'][i][5])
            self.zperson3853.append( x['BOSSOBJECT'][i][6])
            self.type3853.append( x['BOSSOBJECT'][i][7])
            self.comments3853.append( x['BOSSOBJECT'][i][8])
        args = n.where(x['BOSSOBJECT']['plate'] == 3855)[0]
        self.fibers3855 = []
        self.zpipe3855 = []
        self.zperson3855 = []
        self.type3855 = []
        self.comments3855 = []
        for i in args:
            self.fibers3855.append( x['BOSSOBJECT'][i][2])
            self.zpipe3855.append( x['BOSSOBJECT'][i][5])
            self.zperson3855.append( x['BOSSOBJECT'][i][6])
            self.type3855.append( x['BOSSOBJECT'][i][7])
            self.comments3855.append( x['BOSSOBJECT'][i][8])
        args = n.where(x['BOSSOBJECT']['plate'] == 3856)[0]
        self.fibers3856 = []
        self.zpipe3856 = []
        self.zperson3856 = []
        self.type3856 = []
        self.comments3856 = []
        for i in args:
            self.fibers3856.append( x['BOSSOBJECT'][i][2])
            self.zpipe3856.append( x['BOSSOBJECT'][i][5])
            self.zperson3856.append( x['BOSSOBJECT'][i][6])
            self.type3856.append( x['BOSSOBJECT'][i][7])
            self.comments3856.append( x['BOSSOBJECT'][i][8])
        args = n.where(x['BOSSOBJECT']['plate'] == 3860)[0]
        self.fibers3860 = []
        self.zpipe3860 = []
        self.zperson3860 = []
        self.type3860 = []
        self.comments3860 = []
        for i in args:
            self.fibers3860.append( x['BOSSOBJECT'][i][2])
            self.zpipe3860.append( x['BOSSOBJECT'][i][5])
            self.zperson3860.append( x['BOSSOBJECT'][i][6])
            self.type3860.append( x['BOSSOBJECT'][i][7])
            self.comments3860.append( x['BOSSOBJECT'][i][8])


    def get_vifibers(self,plate):
        # Set self.fibers to yanny info for a given plate
        if plate == 3686: self.vifibers = self.fibers3686
        elif plate == 3687: self.vifibers = self.fibers3687
        elif plate == 3804: self.vifibers = self.fibers3804
        elif plate == 3805: self.vifibers = self.fibers3805
        elif plate == 3853: self.vifibers = self.fibers3853
        elif plate == 3855: self.vifibers = self.fibers3855
        elif plate == 3856: self.vifibers = self.fibers3856
        elif plate == 3860: self.vifibers = self.fibers3860


    def get_zperson(self,plate):
        # Set self.zperson to yanny info for a given plate
        if plate == 3686: self.zperson = self.zperson3686
        elif plate == 3687: self.zperson = self.zperson3687
        elif plate == 3804: self.zperson = self.zperson3804
        elif plate == 3805: self.zperson = self.zperson3805
        elif plate == 3853: self.zperson = self.zperson3853
        elif plate == 3855: self.zperson = self.zperson3855
        elif plate == 3856: self.zperson = self.zperson3856
        elif plate == 3860: self.zperson = self.zperson3860


    def get_zpipe(self,plate):
        # Set self.zpipe to yanny info for a given plate
        if plate == 3686: self.zpipe = self.zpipe3686
        elif plate == 3687: self.zpipe = self.zpipe3687
        elif plate == 3804: self.zpipe = self.zpipe3804
        elif plate == 3805: self.zpipe = self.zpipe3805
        elif plate == 3853: self.zpipe = self.zpipe3853
        elif plate == 3855: self.zpipe = self.zpipe3855
        elif plate == 3856: self.zpipe = self.zpipe3856
        elif plate == 3860: self.zpipe = self.zpipe3860


    def get_vitype(self,plate):
        # Set self.vitype to yanny info for a given plate
        if plate == 3686: self.vitype = self.type3686
        elif plate == 3687: self.vitype = self.type3687
        elif plate == 3804: self.vitype = self.type3804
        elif plate == 3805: self.vitype = self.type3805
        elif plate == 3853: self.vitype = self.type3853
        elif plate == 3855: self.vitype = self.type3855
        elif plate == 3856: self.vitype = self.type3856
        elif plate == 3860: self.vitype = self.type3860


    def get_comments(self,plate):
        # Set self.comments to yanny info for a given plate
        if plate == 3686: self.comments = self.comments3686
        elif plate == 3687: self.comments = self.comments3687
        elif plate == 3804: self.comments = self.comments3804
        elif plate == 3805: self.comments = self.comments3805
        elif plate == 3853: self.comments = self.comments3853
        elif plate == 3855: self.comments = self.comments3855
        elif plate == 3856: self.comments = self.comments3856
        elif plate == 3860: self.comments = self.comments3860

    def get_all_yanny(self,plate):
        # Call all of the above self.get_XXX() methods in one fell swoop
        self.get_vifibers(plate)
        self.get_zperson(plate)
        self.get_zpipe(plate)
        self.get_vitype(plate)
        self.get_comments(plate)


    def read_redmonster(self,plate):
        # Read in the redmonster output file for a given plate
        redmonsterpath = join( self.redmonster_spectro_redux, '%s' % plate,
                              '%s' % self.version, 'redmonster-%s-%s.fits' %
                              (plate,self.mjds[plate]) )
        hdu = fits.open(redmonsterpath)
        self.rm_z1 = hdu[1].data.Z1
        self.rm_zerr1 = hdu[1].data.Z_ERR1
        # +1 here because rm fibers are 0-based and idlspec2d are 1-based
        self.rm_fibers = hdu[1].data.FIBERID + 1
        self.rm_type = hdu[1].data.CLASS
        self.rm_zwarning = hdu[1].data.ZWARNING


    def read_spPlate(self,plate):
        # Read in the spPlate file for a given plate
        spPlatepath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version,
                           '%s' % plate, 'spPlate-%s-%s.fits' %
                           (plate, self.mjds[plate]) )
        hdu = fits.open(spPlatepath)
        self.boss_target1 = hdu[5].data.BOSS_TARGET1


    def read_spZbest(self,plate):
        # Read in the spZbest file for a given plate
        spZbestpath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version,
                           '%s' % plate, '%s' % self.version,
                           'spZbest-%s-%s.fits' % (plate, self.mjds[plate]) )
        hdu = fits.open(spZbestpath)
        self.sn_median = hdu[1].data.SN_MEDIAN[:,2:]
        # In i-band, note conversion from nanomaggies to magnitudes
        self.spectroflux = 22.5 - 2.5*n.log10(hdu[1].data.SPECTROFLUX)


    def get_cmass(self):
        # Return (0-based) indices of CMASS targets
        return n.where( self.boss_target1 & 2 == 2 )[0].tolist()


    def get_lowz(self):
        # Return (0-based indices) of LOWZ targets
        return n.where( self.boss_target1 & 1 == 1 )[0].tolist()


    def get_okay_cmass(self):
        # Return (0-based) indices of CMASS targets that have the
        # yanny comment 'v5_4_9 ok' and imag <= 21.5
        # self.get_fibers() and self.get_comments() need to have
        # already been called on this plate for this method to work properly
        #
        # -1 due to fibers being 1-based and python using 0-based
        okay_fibers = (n.asarray(self.vifibers)[n.where(n.asarray(self.comments)
                == 'v5_4_9 ok')[0].tolist()]-1).tolist()
        return n.asarray(okay_fibers)[n.where( (self.boss_target1[okay_fibers]
                                                & 2 == 2) &
                                              (self.spectroflux[okay_fibers]\
                                               [:,3] <= 21.5) )[0].tolist()]\
                                                .tolist()


    def get_okay_lowz(self):
        # Return (0-based) indices of LOWZ targets that have the yanny comment 'v5_4_9 ok' and imag <= 21.5
        # self.get_fibers() and self.get_comments() (or, equivalently, self.get_all_yanny() ) need to have already been called on this plate
        okay_fibers = (n.asarray(self.vifibers)[n.where(n.asarray(self.comments) == 'v5_4_9 ok')[0].tolist()]-1).tolist() # -1 due to fibers being 1-based and python using 0-based
        return n.asarray(okay_fibers)[n.where( (self.boss_target1[okay_fibers] & 1 == 1) & (self.spectroflux[okay_fibers][:,3] <= 21.5) )[0].tolist()].tolist()


    def count_total_targets(self):
        # Prints the total number of visually inspected targets
        count = 0
        for plate in self.plates:
            self.get_all_yanny(plate)
            count += len(self.vifibers)
        print(count)


    def cmass_completeness(self):
        # Prints percent of all CMASS targets with rm_zwarning == 0
        vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_cmass()
            vals.append( float(len(n.where( self.rm_zwarning[fibers] == 0
                                           )[0].tolist())) / float(len(fibers)))
        avg = n.sum(vals) / float(len(vals))
        print(avg)


    def lowz_completeness(self):
        # Prints percent of all LOWZ targets with rm_zwarning == 0
        vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_lowz()
            vals.append( float(len(n.where( self.rm_zwarning[fibers] == 0
                                           )[0].tolist())) / float(len(fibers)))
        avg = n.sum(vals) / float(len(vals))
        print(avg)


    def cmass_galaxy_completeness(self):
        # Prints percent of all CMASS targets that have rm_warning == 0
        # and were classified as 'ssp_galaxy_glob'
        #vals = []
        count = 0
        total = 0
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_cmass()
            for fiber in fibers:
                if self.rm_type[fiber] == 'ssp_galaxy_glob':
                    total += 1
                    if self.rm_zwarning[fiber] == 0:
                        count += 1
        #avg = n.sum(vals) / float(len(vals))
        avg = float(count) / float(total)
        print(count)
        print(total)
        print(avg)


    def lowz_galaxy_completeness(self):
        # Prints percent of all LOWZ targets that have rm_zwarning == 0
        # and were classified as 'ssp_galaxy_glob'
        #vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_lowz()
            for fiber in fibers:
                if self.rm_type[fiber] == 'ssp_galaxy_glob':
                    total += 1
                    if self.rm_zwarning[fiber] == 0:
                        count += 1
        #avg = n.sum(vals) / float(len(vals))
        avg = float(count) / float(total)
        print(count)
        print(total)
        print(avg)

    def count_okay_cmass_fibers(self):
        # Prints number of CMASS targets with yanny comment
        # 'v5_4_9 ok' and imag <= 21.5
        count = 0
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            self.get_all_yanny(plate)
            count += len(self.get_okay_cmass())
        print(count)

    def cmass_okay_completeness(self):
        # Prints fraction of CMASS targets having yanny comment
        # 'v5_4_9 ok' and imag <= 21.5 that have rm_zwarning == 0
        count = 0
        total = 0
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            self.get_all_yanny(plate)
            fibers = self.get_okay_cmass()
            total += len(fibers)
            count += len(n.where(self.rm_zwarning[fibers] == 0)[0].tolist())
        print('%s out of %s' % (count,total))
        print(float(count) / float(total))

    def cmass_okay_galaxy_completeness(self):
        # Prints fraction of targets classified by RM as 'ssp_em_galaxies'
        # in the subset of CMASS targets having yanny comment 'v5_4_9 ok'
        # and imag <= 21.5
        count = 0
        total = 0
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.get_all_yanny(plate)
            self.read_spZbest(plate)
            fibers = self.get_okay_cmass()
            total += len(fibers)
            count += len( n.where( self.rm_type[fibers] ==
                                  'ssp_galaxy_glob')[0].tolist() )
        print('%s out of %s' % (count,total))
        print(float(count) / float(total))


    def dz_to_dv(self, z, dz):
        # Convert redshift error dz to velocity error dv
        c_kms = 299792.458 # speed of light in km s^-1
        return (dz * c_kms) / (1 + z)


    def redshift_bin_fibers(self, fibers, zmin, zmax):
        # Return subset of fibers in redshift range [zmin,zmax]
        bin_fibers = []
        for fiber in fibers:
            if (self.rm_z1[fiber] >= zmin) & (self.rm_z1[fiber] <= zmax):
                bin_fibers.append(fiber)
        return bin_fibers


    def logdv_histos(self, nbins=12):
        # Make histograms of log10(dv) in redshift bins for
        # LOWZ and CMASS galaxies
        colors = ['purple', 'cyan', 'blue', 'red', 'gold', 'lime']
        labels = ['0.1<z<0.2','0.2<z<0.3','0.3<z<0.4','0.4<z<0.5']
        f = p.figure()
        ax1 = f.add_subplot(1,2,1)
        for j,zmin in enumerate(n.linspace(.1,.4,4)):
            zmax = zmin + .1
            errors = n.array([])
            count = 0
            for plate in self.plates:
                self.read_redmonster(plate)
                self.read_spPlate(plate)
                self.read_spZbest(plate)
                self.get_all_yanny(plate)
                fibers = self.get_okay_lowz()
                fibers = self.redshift_bin_fibers(fibers, zmin, zmax)
                count += len(fibers)
                errors = n.append(errors, self.rm_zerr1[fibers])
            errors = self.dz_to_dv(errors)
            errors = n.log10(errors)
            hist,binedges = n.histogram(errors, bins=nbins)
            bins = n.zeros(nbins)
            for i in range(nbins):
                bins[i] = (binedges[i+1]+binedges[i])/2.
            normhist = hist / float(count)
            p.plot(bins,normhist,drawstyle='steps-mid', color=colors[j],
                   label=labels[j])
        p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        p.title('LOWZ Sample', size=18)
        p.legend()
        p.axis([.55,2,0,.4])
        ax2 = f.add_subplot(1,2,2)
        labels = ['0.4<z<0.5','0.5<z<0.6','0.6<z<0.7','0.7<z<0.8']
        nbins = 25
        for j,zmin in enumerate(n.linspace(.4,.7,4)):
            #import pdb; pdb.set_trace()
            zmax = zmin + .1
            errors = n.array([])
            count = 0
            for plate in self.plates:
                self.read_redmonster(plate)
                self.read_spPlate(plate)
                self.read_spZbest(plate)
                self.get_all_yanny(plate)
                fibers = self.get_okay_cmass()
                fibers = self.redshift_bin_fibers(fibers, zmin, zmax)
                count += len(fibers)
                errors = n.append(errors,self.rm_zerr1[fibers])
                #errors.append(self.rm_zerr1[fibers].tolist())
            errors = self.dz_to_dv(errors)
            errors = n.log10(errors)
            hist,binedges = n.histogram(errors, bins=nbins)
            bins = n.zeros(nbins)
            for i in range(nbins):
                bins[i] = (binedges[i+1]+binedges[i])/2.
            normhist = hist / float(count)
            p.plot(bins,normhist,drawstyle='steps-mid', color=colors[j],
                   label=labels[j])
        p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        p.title('CMASS Sample', size=18)
        p.axis([.9,2.4,0,.3])
        p.legend()
        p.subplots_adjust(wspace = .35)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/\
                  dv_histo_cmass.pdf')


    def identify_catastrophic_failures(self):
        # Identify fibers in 'okay CMASS' sample with zwarning == 0
        # and abs(z_rm - z_person) > .005
        self.bad_fibers = []
        self.bad_plates = []
        self.bad_rm_z = []
        self.bad_zperson = []
        self.bad_type = []
        count_bad = 0
        total = 0
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            self.get_all_yanny(plate)
            fibers = self.get_okay_cmass()
            for fiber in fibers:
                # +1 to go from python indexing to boss fiber conventions
                if (fiber+1) in self.vifibers:
                    total += 1
                    vi_index = n.where( n.asarray(self.vifibers) ==
                                       (fiber+1) )[0][0]
                    if self.rm_zwarning[fiber] == 0:
                        if self.rm_type[fiber] == 'ssp_galaxy_glob':
                            if n.abs(self.rm_z1[fiber] -
                                     self.zperson[vi_index]) >= 0.005:
                                self.bad_plates.append(plate)
                                self.bad_fibers.append(fiber)
                                self.bad_rm_z.append(self.rm_z1[fiber])
                                self.bad_zperson.append(self.zperson[vi_index])
                                self.bad_type.append(self.rm_type[fiber])
                                count_bad += 1
        print('%s catastrophic failures out of %s fibers, or %s PERCENT \
                (not fraction!) of the total' % (count_bad,total,
                                                 (count_bad/float(total))*100))
        for i,fiber in enumerate(self.bad_fibers):
            print('Plate %s, fiber %s, redmonster z = %s, zperson = %s' % \
                    (self.bad_plates[i],fiber,self.bad_rm_z[i],
                     self.bad_zperson[i]))


    def identify_unclear_impurities(self):
        # Identify fibers that have zwarning == 0
        # but no confident visual redshift
        pass


    def identify_recoverable_incompleteness(self):
        # Identify fibers with confident visual redshift and 'galaxy'
        # classification but have zwarning != 0 or rm_type == 'star'
        # Also makes plot of z_visual vs z_rm for all identified fibers
        self.recoverable_fibers = []
        self.recoverable_plates = []
        self.recoverable_rm_z = []
        self.recoverable_rm_type = []
        self.recoverable_zperson = []
        count_recoverable = 0
        total = 0
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            self.get_all_yanny(plate)
            fibers = self.get_okay_cmass()
            for fiber in fibers:
                # +1 to go from python indexing to boss fiber conventions
                if (fiber+1) in self.vifibers:
                    total += 1
                    vi_index = n.where( n.asarray(self.vifibers) ==
                                       (fiber+1) )[0][0]
                    if (self.rm_zwarning[fiber] != 0) | \
                            (self.rm_type[fiber] != 'ssp_galaxy_glob'):
                        if (self.zperson[vi_index] != -9) & \
                                (self.vitype[vi_index] == 4):
                            self.recoverable_fibers.append(fiber)
                            self.recoverable_plates.append(plate)
                            self.recoverable_rm_z.append(self.rm_z1[fiber])
                            self.recoverable_rm_type.append(self.rm_type[fiber])
                            self.recoverable_zperson.append(
                                    self.zperson[vi_index])
                            count_recoverable += 1
        print('%s recoverable failures out of %s fibers, or %s PERCENT \
                (not fraction!) of the total' % (count_recoverable,total,
                                                 (count_recoverable /
                                                  float(total))*100))
        for i,fiber in enumerate(self.recoverable_fibers):
            print('Plate %s, fiber %s, redmonster z = %s, \
                    redmonster class = %s, zperson = %s' % \
                    (self.recoverable_plates[i],fiber,self.recoverable_rm_z[i],
                     self.recoverable_rm_type[i], self.recoverable_zperson[i]))
        big_diff_num = len( n.where( n.abs(n.asarray(self.recoverable_rm_z) -
                                           n.asarray(self.recoverable_zperson))
                                    >= .005 )[0] )
        f = p.figure()
        ax1 = f.add_subplot(1,1,1)
        p.plot(self.recoverable_rm_z,self.recoverable_zperson, 'k.')
        p.plot(n.linspace(0,1,1000),n.linspace(0,1,1000),'red')
        p.axis([-.5,3,0,1])
        p.xlabel(r'$z_{redmonster}$',size=16)
        p.ylabel(r'$z_{visual}$',size=16)
        p.title('Objects with "recoverable" redshifts', size=18)
        p.text(1.25, .2, '%s out of %s fibers with confident visual' %
               (count_recoverable,total), fontsize=10)
        p.text(1.25,.15, 'redshift and called "galaxy" but have', size=10)
        p.text(1.25, .1, 'zwarning > 0 or class != "galaxy". Of', size=10)
        p.text(1.25,.05, 'these, %s have $\delta z > 0.005$.' % (big_diff_num),
               size=10)
        #p.savefig('recov.pdf')


    def cmass_failure_vs_sn(self,sn_max=7,nbins=29):
        # Makes plot of CMASS failure rate (zwarning > 0) vs
        # median S/N in r-, i-, and z-bands
        f = p.figure()
        ax = f.add_subplot(1,1,1)
        total = 0
        bad_fibers = []
        bad_r_sn = []
        bad_i_sn = []
        bad_z_sn = []
        r_sn = []
        i_sn = []
        z_sn = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            self.get_all_yanny(plate)
            fibers = self.get_cmass()
            for fiber in fibers:
                if (self.sn_median[fiber,0] <= sn_max):
                    total += 1.
                    r_sn.append(self.sn_median[fiber,0])
                    i_sn.append(self.sn_median[fiber,1])
                    z_sn.append(self.sn_median[fiber,2])
                    if (self.rm_zwarning[fiber] > 0):
                        bad_fibers.append(fiber)
                        bad_r_sn.append(self.sn_median[fiber,0])
                        bad_i_sn.append(self.sn_median[fiber,1])
                        bad_z_sn.append(self.sn_median[fiber,2])
        nbinsarr = n.linspace(0,sn_max,nbins+1)
        rtotal,rbinedges = n.histogram(r_sn,bins=nbinsarr)
        itotal,ibinedges = n.histogram(i_sn,bins=nbinsarr)
        ztotal,zbinedges = n.histogram(z_sn,bins=nbinsarr)
        rhist,rbinedges = n.histogram(bad_r_sn,bins=nbinsarr)
        ihist,ibinedges = n.histogram(bad_i_sn,bins=nbinsarr)
        zhist,zbinedges = n.histogram(bad_z_sn,bins=nbinsarr)
        rbins = n.zeros(nbins)
        ibins = n.zeros(nbins)
        zbins = n.zeros(nbins)
        for i in range(nbins):
            rbins[i] = (rbinedges[i+1]+rbinedges[i])/2.
            ibins[i] = (ibinedges[i+1]+ibinedges[i])/2.
            zbins[i] = (zbinedges[i+1]+zbinedges[i])/2.
        rhist = rhist / rtotal.astype(float)
        ihist = ihist / itotal.astype(float)
        zhist = zhist / ztotal.astype(float)
        for i in range(nbins):
            if i != 0 and i != (nbins-1):
                if isnan(rhist[i]):
                    try:
                        rhist[i] = (rhist[i-1] + rhist[i+1]) / 2.
                    except:
                        rhist[i] = 0
                if isnan(ihist[i]):
                    try:
                        ihist[i] = (ihist[i-1] + ihist[i+1]) / 2.
                    except:
                        ihist[i] = 0
                if isnan(zhist[i]):
                    try:
                        zhist[i] = (zhist[i-1] + zhist[i+1]) / 2.
                    except:
                        zhist[i] = 0
        rhist = convolve(rhist,Box1DKernel(2))
        ihist = convolve(ihist,Box1DKernel(2))
        zhist = convolve(zhist,Box1DKernel(2))
        p.plot(rbins,rhist,color='purple',label='r-band')
        p.plot(ibins,ihist,color='blue',label='i-band')
        p.plot(zbins,zhist,color='cyan',label='z-band')
        ax.set_yscale('log')
        p.xlabel(r'Median S/N per 69 km s$^{-1}$ coadded pixel',size=14)
        p.ylabel(r'CMASS failure rate', size=14)
        #print rbins
        #print rhist
        #print rtotal
        p.legend()
        p.savefig('failure_vs_sn.pdf')


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---------------FULL SEQUELS METHODS ONLY BELOW THIS LINE ---------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


    def read_redmonster_all(self,plate, mjd):
        # Read a redmonster file in the context of looking
        #  at SEQUELS LRG dataset
        redmonsterpath = join( self.redmonster_spectro_redux, '%s' % plate,
                              '%s' % self.version, 'redmonster-%s-%s.fits' %
                              (plate,mjd) )
        #paths = []
        #for path in iglob(redmonsterpath):
        #    paths.append(path)
        #paths.sort()
        #hdu = fits.open(paths[0])
        hdu = fits.open(redmonsterpath)
        self.rm_z1 = hdu[1].data.Z1
        self.rm_zerr1 = hdu[1].data.Z_ERR1
        self.rm_fibers = hdu[1].data.FIBERID
        self.rm_type = hdu[1].data.CLASS1
        self.rm_type2 = hdu[1].data.CLASS2
        self.rm_zwarning = hdu[1].data.ZWARNING



    def read_redmonster_summary_file(self):
        # Read the redmonster summary file
        redmonster_spectro_redux = join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly1' % self.version)
        summary_path = join( redmonster_spectro_redux,
                            'redmonsterAll-%s.fits' % self.version )
        hdu = fits.open(summary_path)
        self.rm_z1 = hdu[1].data.Z
        self.rm_zerr1 = hdu[1].data.Z_ERR
        self.rm_type = hdu[1].data.CLASS
        self.rm_zwarning = hdu[1].data.ZWARNING
        self.rm_fibers_summary = hdu[1].data.FIBERID
        self.rm_plates_summary = hdu[1].data.PLATE
        self.rm_mjds_summary = hdu[1].data.MJD
        self.rm_rchi2s = hdu[1].data.MINRCHI2
        self.rm_dof = hdu[1].data.DOF
        self.rm_rchi2diff = hdu[1].data.RCHI2DIFF
        self.rm_chi2_null = hdu[1].data.CHI2NULL
        self.rm_sn2_data = hdu[1].data.SN2DATA


    def read_spPlate_all(self,plate, mjd=None):
        # Read in the spPlate file for a given plate in the context of
        # the entire DR10 dataset
        if mjd is not None:
            hdu = fits.open( join( environ['BOSS_SPECTRO_REDUX'], '%s' %
                                  self.version, '%s' % plate,
                                  'spPlate-%s-%s.fits' % (plate,mjd) ) )
            try: self.eboss_target0 = hdu[5].data.EBOSS_TARGET0
            except: pass
            try: self.eboss_target1 = hdu[5].data.EBOSS_TARGET1
            except: pass
        else:
            globpath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version,
                            '%s' % plate, 'spPlate-%s-*.fits' % plate )
            spPlatepaths = []
            for spPlatepath in iglob(globpath):
                spPlatepaths.append(spPlatepath)
            spPlatepaths.sort()
            hdu = fits.open(spPlatepaths[0])
            try: self.eboss_target0 = hdu[5].data.EBOSS_TARGET0
            except: pass
            try: self.eboss_target1 = hdu[5].data.EBOSS_TARGET1
            except: pass


    def read_spZbest_all(self,plate,mjd=None):
        # Read in the spZbest file for a given plate
        if mjd is not None:
            hdu = fits.open(join(environ['BOSS_SPECTRO_REDUX'],
                                 '%s' % self.version, '%s' % plate,
                                 '%s' % self.version,
                                 'spZbest-%s-%s.fits' % (plate,mjd) ))
            self.sn_median = hdu[1].data.SN_MEDIAN[:,2:]
            # In i-band, note conversion from nanomaggies to magnitudes
            self.spectroflux = 22.5 - 2.5*n.log10(hdu[1].data.SPECTROFLUX)
            self.idl_rchi2s = hdu[1].data.RCHI2
            self.idl_dof = hdu[1].data.DOF
            self.idl_rchi2diff = hdu[1].data.RCHI2DIFF_NOQSO
            #self.modelmag = hdu[1].data.MODELMAG[:,2:]
            #self.extinction = hdu[1].data.EXTINCTION[:,2:]
        else:
            globpath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version,
                            '%s' % plate, '%s' % self.version,
                            'spZbest-%s-*.fits' % plate )
            spZbestpaths = []
            for spZbestpath in iglob(globpath):
                spZbestpaths.append(spZbestpath)
            spZbestpaths.sort()
            hdu = fits.open(spZbestpaths[0])
            self.sn_median = hdu[1].data.SN_MEDIAN[:,2:]
            # In i-band, note conversion from nanomaggies to magnitudes
            self.spectroflux = 22.5 - 2.5*n.log10(hdu[1].data.SPECTROFLUX)
            self.idl_dof = hdu[1].data.DOF
            self.idl_rchi2diff = hdu[1].data.RCHI2DIFF
            #self.modelmag = hdu[1].data.MODELMAG[:,2:]
            #eself.extinction = hdu[1].data.EXTINCTION[:,2:]


    def sequels_completeness_all(self):
        # Prints percent of all SEQUELS LRG targets with rm_zwarning == 0
        count = 0
        total = 0
        self.read_redmonster_summary_file()
        for zwarn in self.rm_zwarning:
            total += 1
            if zwarn & 4 == 0:
                count += 1
        avg = float(count) / float(total)
        print(count)
        print(total)
        print(avg)


    def sequels_galaxy_completeness_all(self):
        # Prints percent of all DR10 CMASS targets that have
        # rm_warning == 0 and were classified as 'ssp_galaxy_glob'
        count = 0
        total = 0
        #globpath = join( self.redmonster_spectro_redux, '*')
        #for path in iglob(globpath):
        #    plate = basename(path)
        #    self.read_spPlate_all(plate)
        #    self.read_redmonster_all(plate)
        #    fibers = self.get_cmass()
        #    for fiber in fibers:
        #        total += 1
        #        if (self.rm_zwarning[fiber] == 0) & \
        #                (self.rm_type[fiber] == 'ssp_galaxy_glob'):
        #            count += 1
        self.read_redmonster_summary_file()
        for i,zwarn in enumerate(self.rm_zwarning):
            #if zwarn == 0:
            if (zwarn & 4) == 0:
                total += 1
                if self.rm_type[i] == 'ssp_galaxy_glob':
                    count += 1
        avg = float(count) / float(total)
        print(count)
        print(total)
        print(avg)


    def logdv_vs_z_histos(self, nbins=12):
        # Make histograms of log10(dv) in redshift bins for
        # LOWZ and CMASS galaxies
        colors = [
                  'tomato','sage','cornflowerblue','sandybrown',
                  'mediumpurple','grey'
                  ]
        labels = ['0.1<z<0.2','0.2<z<0.3','0.3<z<0.4','0.4<z<0.5']
        f = p.figure()
        '''
        ax1 = f.add_subplot(1,2,1)
        for j,zmin in enumerate(n.linspace(.1,.4,4)):
            zmax = zmin + .1
            errors = n.array([])
            count = 0
            for plate in self.plates:
                self.read_redmonster(plate)
                self.read_spPlate(plate)
                self.read_spZbest(plate)
                self.get_all_yanny(plate)
                fibers = self.get_okay_lowz()
                fibers = self.redshift_bin_fibers(fibers, zmin, zmax)
                count += len(fibers)
                errors = n.append(errors, self.rm_zerr1[fibers])
            errors = self.dz_to_dv(errors)
            errors = n.log10(errors)
            hist,binedges = n.histogram(errors, bins=nbins)
            bins = n.zeros(nbins)
            for i in xrange(nbins):
                bins[i] = (binedges[i+1]+binedges[i])/2.
            normhist = hist / float(count)
            p.plot(bins,normhist,drawstyle='steps-mid', color=colors[j],
                   label=labels[j])
        p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        p.title('LOWZ Sample', size=18)
        p.legend()
        p.axis([.55,2,0,.4])
        '''
        ax2 = f.add_subplot(1,1,1)
        labels = ['0.6<z<0.7','0.7<z<0.8','0.8<z<0.9','0.9<z<1.0']
        nbins = 25
        for j,zmin in enumerate(n.linspace(.6,.9,4)):
            #import pdb; pdb.set_trace()
            zmax = zmin + .1
            errors = n.array([])
            zs = n.array([])
            count = 0
            '''
            for plate in self.plates:
                self.read_redmonster(plate)
                #self.read_spPlate(plate)
                #self.read_spZbest(plate)
                #self.get_all_yanny(plate)
                fibers = self.get_okay_cmass()
                fibers = self.redshift_bin_fibers(fibers, zmin, zmax)
                count += len(fibers)
                errors = n.append(errors,self.rm_zerr1[fibers])
            '''
            self.read_redmonster_summary_file()
            for i,z in enumerate(self.rm_z1):
                if (z >= zmin) & (z <= zmax):
                    if (self.rm_type[i] == 'ssp_galaxy_glob') & \
                            (self.rm_zwarning[i] == 0) & (self.rm_zerr1[i] > 0):
                        count += 1
                        errors = n.append(errors,self.rm_zerr1[i])
                        zs = n.append(zs,z)
            #errors.append(self.rm_zerr1[fibers].tolist())
            errors = self.dz_to_dv(zs, errors)
            print(zmin, zmax, n.mean(errors), n.std(errors))
            errors = n.log10(errors)
            hist,binedges = n.histogram(errors, bins=nbins)
            bins = n.zeros(nbins)
            for i in range(nbins):
                bins[i] = (binedges[i+1]+binedges[i])/2.
            normhist = hist / float(count)
            p.plot(bins,normhist,drawstyle='steps-mid', color=colors[j],
                   label=labels[j])
        p.minorticks_on()
        p.xlabel(r'$\log_{10} \delta v$ (km s$^{-1}$)', size=16)
        p.ylabel(r'Fraction per bin in $\log_{10} \delta v$', size=16)
        p.title('SEQUELS LRGs', size=18)
        p.axis([.5,3.0,0,.3])
        p.legend()
        p.subplots_adjust(wspace = .35)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/dv_vs_z_histos.pdf')
        p.clf()


    def sequels_failure_vs_sn_all(self,sn_max=4.5,nbins=18):
        # Makes plot of SEQUELS LRG target failure rate (zwarning > 0)
        # vs median S/N in r-, i-, and z-bands
        f = p.figure()
        ax = f.add_subplot(1,1,1)
        total = 0
        bad_fibers = []
        bad_r_sn = []
        bad_i_sn = []
        bad_z_sn = []
        r_sn = []
        i_sn = []
        z_sn = []
        rmax = 0
        imax = 0
        zmax = 0
        globpath = join( self.redmonster_spectro_redux,'*')
        openplate = 0
        openmjd = 0
        self.read_redmonster_summary_file()
        for i,fiber in enumerate(self.rm_fibers_summary):
            plate = self.rm_plates_summary[i]
            mjd = self.rm_mjds_summary[i]
            print('%s-%s-%s' % (plate,fiber,mjd))
            if (openplate != plate) and (openmjd != mjd):
                self.read_spZbest_all(plate,mjd)
                openplate = plate
                openmjd = mjd
            if (self.sn_median[fiber,0] <= sn_max):
                total += 1
                r_sn.append(self.sn_median[fiber,0])
                if self.sn_median[fiber,0] > rmax:
                    rmax = self.sn_median[fiber,0]
                i_sn.append(self.sn_median[fiber,1])
                if self.sn_median[fiber,1] > imax:
                    imax = self.sn_median[fiber,1]
                z_sn.append(self.sn_median[fiber,2])
                if self.sn_median[fiber,2] > zmax:
                    zmax = self.sn_median[fiber,2]
                if (self.rm_zwarning[i] > 0):
                    bad_fibers.append(fiber)
                    bad_r_sn.append(self.sn_median[fiber,0])
                    bad_i_sn.append(self.sn_median[fiber,1])
                    bad_z_sn.append(self.sn_median[fiber,2])
        nbinsarr = n.linspace(0,sn_max,nbins+1)
        rtotal,rbinedges = n.histogram(r_sn,bins=nbinsarr)
        itotal,ibinedges = n.histogram(i_sn,bins=nbinsarr)
        ztotal,zbinedges = n.histogram(z_sn,bins=nbinsarr)
        rhist,rbinedges = n.histogram(bad_r_sn,bins=nbinsarr)
        ihist,ibinedges = n.histogram(bad_i_sn,bins=nbinsarr)
        zhist,zbinedges = n.histogram(bad_z_sn,bins=nbinsarr)
        rbins = n.zeros(nbins)
        ibins = n.zeros(nbins)
        zbins = n.zeros(nbins)
        for i in range(nbins):
            rbins[i] = (rbinedges[i+1]+rbinedges[i])/2.
            ibins[i] = (ibinedges[i+1]+ibinedges[i])/2.
            zbins[i] = (zbinedges[i+1]+zbinedges[i])/2.
        rhist = rhist / rtotal.astype(float)
        ihist = ihist / itotal.astype(float)
        zhist = zhist / ztotal.astype(float)
        for i in range(nbins):
            if i != 0 and i != (nbins-1):
                if isnan(rhist[i]):
                    try:
                        rhist[i] = (rhist[i-1] + rhist[i+1]) / 2.
                    except:
                        rhist[i] = 0
                if isnan(ihist[i]):
                    try:
                        ihist[i] = (ihist[i-1] + ihist[i+1]) / 2.
                    except:
                        ihist[i] = 0
                if isnan(zhist[i]):
                    try:
                        zhist[i] = (zhist[i-1] + zhist[i+1]) / 2.
                    except:
                        zhist[i] = 0

        p.plot(rbins,rhist,color='purple',label='r-band', drawstyle='steps-mid')
        p.plot(ibins,ihist,color='blue',label='i-band', drawstyle='steps-mid')
        p.plot(zbins,zhist,color='cyan',label='z-band', drawstyle='steps-mid')
        ax.set_yscale('log')
        p.xlabel(r'Median S/N per 69 km s$^{-1}$ coadded pixel',size=14)
        p.ylabel(r'SEQUELS LRG target failure rate', size=14)
        print(rbins)
        print(rhist)
        print(rtotal)
        print(total)
        print(rmax)
        print(imax)
        print(zmax)
        p.legend()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/failure_vs_sn.pdf')
        p.clf()


    def logdv_vs_sn_histos(self, nbins=25):
        # Make histograms of log10(dv) in S/N bins in bands
        # r,i,z for SEQUELS LRG targets
        colors = [
                  'tomato','sage','cornflowerblue','sandybrown','mediumpurple',
                  'grey'
                  ] #['purple', 'cyan', 'blue', 'lime', 'red', 'black']
        labels = [
                  '1.0<S/N<1.5','1.5<S/N<2.0','2.0<S/N<2.5','2.5<S/N<3.0',
                  '3.0<S/N<3.5','3.5<S/N<4.0','4.0<S/N<4.5'
                  ]
        f = p.figure()

        ax1 = f.add_subplot(3,1,1)
        errors1 = n.array([])
        errors2 = n.array([])
        errors3 = n.array([])
        errors4 = n.array([])
        errors5 = n.array([])
        errors6 = n.array([])
        z1 = n.array([])
        z2 = n.array([])
        z3 = n.array([])
        z4 = n.array([])
        z5 = n.array([])
        z6 = n.array([])
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        openplate = 0
        openmjd = 0
        self.read_redmonster_summary_file()
        for j,sn_min in enumerate(n.linspace(1,3.5,6)):
            sn_max = sn_min + .5
            for i,fiber in enumerate(self.rm_fibers_summary):
                plate = self.rm_plates_summary[i]
                mjd = self.rm_mjds_summary[i]
                #print '%s-%s-%s' % (plate,fiber,mjd)
                if (openplate != plate) or (openmjd != mjd):
                    #self.read_spZbest_all(plate,mjd)
                    hduzbest = fits.open(join('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14', '%s' % plate, 'test_dr14', 'spZbest-%s-%s.fits' % (plate, mjd)))
                    self.sn_median = hduzbest[1].data.SN_MEDIAN[:,2:]
                    #self.read_spPlate_all(plate,mjd)
                    hduspplate = fits.open(join('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14', '%s' % plate, 'spPlate-%s-%s.fits' % (plate, mjd)))
                    openplate = plate
                    openmjd = mjd
                if (self.rm_zwarning[i] == 0) & (self.rm_zerr1[i] > 0):
                    if (self.sn_median[fiber][0] >= sn_min) & \
                            (self.sn_median[fiber][0] <= sn_max):
                        if j == 0:
                            errors1 = n.append(errors1,self.rm_zerr1[i])
                            z1 = n.append(z1,self.rm_z1[i])
                            count1 += 1
                        elif j == 1:
                            errors2 = n.append(errors2,self.rm_zerr1[i])
                            z2 = n.append(z2,self.rm_z1[i])
                            count2 += 1
                        elif j == 2:
                            errors3 = n.append(errors3,self.rm_zerr1[i])
                            z3 = n.append(z3,self.rm_z1[i])
                            count3 += 1
                        elif j == 3:
                            errors4 = n.append(errors4,self.rm_zerr1[i])
                            z4 = n.append(z4,self.rm_z1[i])
                            count4 += 1
                        elif j == 4:
                            errors5 = n.append(errors5,self.rm_zerr1[i])
                            z5 = n.append(z5,self.rm_z1[i])
                            count5 += 1
                        elif j == 5:
                            errors6 = n.append(errors6,self.rm_zerr1[i])
                            z6 = n.append(z6,self.rm_z1[i])
                            count6 += 1

        errors1 = self.dz_to_dv(z1,errors1)
        print('r', labels[0], n.mean(errors1), n.std(errors1))
        errors1 = n.log10(errors1)
        hist1,binedges1 = n.histogram(errors1, bins=nbins)
        bins1 = n.zeros(nbins)
        for i in range(nbins):
            bins1[i] = (binedges1[i+1]+binedges1[i])/2.
        normhist1 = hist1 / float(count1)
        p.plot(bins1,normhist1,drawstyle='steps-mid', color=colors[0],
               label=labels[0])
        errors2 = self.dz_to_dv(z2,errors2)
        print('r', labels[1], n.mean(errors2), n.std(errors2))
        errors2 = n.log10(errors2)
        hist2,binedges2 = n.histogram(errors2, bins=nbins)
        bins2 = n.zeros(nbins)
        for i in range(nbins):
            bins2[i] = (binedges2[i+1]+binedges2[i])/2.
        normhist2 = hist2 / float(count2)
        p.plot(bins2,normhist2,drawstyle='steps-mid', color=colors[1],
               label=labels[1])
        errors3 = self.dz_to_dv(z3,errors3)
        print('r', labels[2], n.mean(errors3), n.std(errors3))
        errors3 = n.log10(errors3)
        hist3,binedges3 = n.histogram(errors3, bins=nbins)
        bins3 = n.zeros(nbins)
        for i in range(nbins):
            bins3[i] = (binedges3[i+1]+binedges3[i])/2.
        normhist3 = hist3 / float(count3)
        p.plot(bins3,normhist3,drawstyle='steps-mid', color=colors[2],
               label=labels[2])
        errors4 = self.dz_to_dv(z4,errors4)
        print('r', labels[3], n.mean(errors4), n.std(errors4))
        errors4 = n.log10(errors4)
        hist4,binedges4 = n.histogram(errors4, bins=nbins)
        bins4 = n.zeros(nbins)
        for i in range(nbins):
            bins4[i] = (binedges4[i+1]+binedges4[i])/2.
        normhist4 = hist4 / float(count4)
        p.plot(bins4,normhist4,drawstyle='steps-mid', color=colors[3],
               label=labels[3])
        errors5 = self.dz_to_dv(z5,errors5)
        print('r', labels[4], n.mean(errors5), n.std(errors5))
        errors5 = n.log10(errors5)
        hist5,binedges5 = n.histogram(errors5, bins=nbins)
        bins5 = n.zeros(nbins)
        for i in range(nbins):
            bins5[i] = (binedges5[i+1]+binedges5[i])/2.
        normhist5 = hist5 / float(count5)
        #p.plot(bins5,normhist5,drawstyle='steps-mid', color=colors[4],
                #label=labels[4])
        errors6 = self.dz_to_dv(z6,errors6)
        print('r', labels[5], n.mean(errors6), n.std(errors6))
        errors6 = n.log10(errors6)
        hist6,binedges6 = n.histogram(errors6, bins=nbins)
        bins6 = n.zeros(nbins)
        for i in range(nbins):
            bins6[i] = (binedges6[i+1]+binedges6[i])/2.
        normhist6 = hist6 / float(count6)
        #p.plot(bins6,normhist6,drawstyle='steps-mid', color=colors[5],
                #label=labels[5])
        p.text(0.8, 0.2, '$r$-band', fontsize=12)
        #p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        #p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        #p.title('r-band', size=18)
        p.axis([.5,2.5,0,.25])
        p.legend(prop={'size':6})
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print((count1+count2+count3+count4+count5+count6) / \
                float(self.rm_fibers_summary.shape[0]))

        ax2 = f.add_subplot(3,1,2)
        errors1 = n.array([])
        errors2 = n.array([])
        errors3 = n.array([])
        errors4 = n.array([])
        errors5 = n.array([])
        errors6 = n.array([])
        z1 = n.array([])
        z2 = n.array([])
        z3 = n.array([])
        z4 = n.array([])
        z5 = n.array([])
        z6 = n.array([])
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        openplate = 0
        openmjd = 0
        self.read_redmonster_summary_file()
        for j,sn_min in enumerate(n.linspace(1,3.5,6)):
            sn_max = sn_min + .5
            for i,fiber in enumerate(self.rm_fibers_summary):
                plate = self.rm_plates_summary[i]
                mjd = self.rm_mjds_summary[i]
                #print '%s-%s-%s' % (plate,fiber,mjd)
                if (openplate != plate) and (openmjd != mjd):
                    #self.read_spZbest_all(plate,mjd)
                    #self.read_spPlate_all(plate,mjd)
                    hduzbest = fits.open(join('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14', '%s' % plate, 'test_dr14', 'spZbest-%s-%s.fits' % (plate, mjd)))
                    self.sn_median = hduzbest[1].data.SN_MEDIAN[:,2:]
                    hduspplate = fits.open(join('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14', '%s' % plate, 'spPlate-%s-%s.fits' % (plate, mjd)))
                    openplate = plate
                    openmjd = mjd
                if (self.rm_zwarning[i] == 0) & (self.rm_zerr1[i] > 0):
                    if (self.sn_median[fiber][1] >= sn_min) & \
                            (self.sn_median[fiber][1] <= sn_max):
                        if j == 0:
                            errors1 = n.append(errors1,self.rm_zerr1[i])
                            z1 = n.append(z1,self.rm_z1[i])
                            count1 += 1
                        elif j == 1:
                            errors2 = n.append(errors2,self.rm_zerr1[i])
                            z2 = n.append(z2,self.rm_z1[i])
                            count2 += 1
                        elif j == 2:
                            errors3 = n.append(errors3,self.rm_zerr1[i])
                            z3 = n.append(z3,self.rm_z1[i])
                            count3 += 1
                        elif j == 3:
                            errors4 = n.append(errors4,self.rm_zerr1[i])
                            z4 = n.append(z4,self.rm_z1[i])
                            count4 += 1
                        elif j == 4:
                            errors5 = n.append(errors5,self.rm_zerr1[i])
                            z5 = n.append(z5,self.rm_z1[i])
                            count5 += 1
                        elif j == 5:
                            errors6 = n.append(errors6,self.rm_zerr1[i])
                            z6 = n.append(z6,self.rm_z1[i])
                            count6 += 1

        errors1 = self.dz_to_dv(z1,errors1)
        print('i', labels[0], n.mean(errors1), n.std(errors1))
        errors1 = n.log10(errors1)
        hist1,binedges1 = n.histogram(errors1, bins=nbins)
        bins1 = n.zeros(nbins)
        for i in range(nbins):
            bins1[i] = (binedges1[i+1]+binedges1[i])/2.
        normhist1 = hist1 / float(count1)
        p.plot(bins1,normhist1,drawstyle='steps-mid', color=colors[0],
               label=labels[0])
        errors2 = self.dz_to_dv(z2,errors2)
        print('i', labels[1], n.mean(errors2), n.std(errors2))
        errors2 = n.log10(errors2)
        hist2,binedges2 = n.histogram(errors2, bins=nbins)
        bins2 = n.zeros(nbins)
        for i in range(nbins):
            bins2[i] = (binedges2[i+1]+binedges2[i])/2.
        normhist2 = hist2 / float(count2)
        p.plot(bins2,normhist2,drawstyle='steps-mid', color=colors[1],
               label=labels[1])
        errors3 = self.dz_to_dv(z3,errors3)
        print('i', labels[2], n.mean(errors3), n.std(errors3))
        errors3 = n.log10(errors3)
        hist3,binedges3 = n.histogram(errors3, bins=nbins)
        bins3 = n.zeros(nbins)
        for i in range(nbins):
            bins3[i] = (binedges3[i+1]+binedges3[i])/2.
        normhist3 = hist3 / float(count3)
        p.plot(bins3,normhist3,drawstyle='steps-mid', color=colors[2],
               label=labels[2])
        errors4 = self.dz_to_dv(z4,errors4)
        print('i', labels[3], n.mean(errors4), n.std(errors4))
        errors4 = n.log10(errors4)
        hist4,binedges4 = n.histogram(errors4, bins=nbins)
        bins4 = n.zeros(nbins)
        for i in range(nbins):
            bins4[i] = (binedges4[i+1]+binedges4[i])/2.
        normhist4 = hist4 / float(count4)
        p.plot(bins4,normhist4,drawstyle='steps-mid', color=colors[3],
               label=labels[3])
        errors5 = self.dz_to_dv(z5,errors5)
        print('i', labels[4], n.mean(errors5), n.std(errors5))
        errors5 = n.log10(errors5)
        hist5,binedges5 = n.histogram(errors5, bins=nbins)
        bins5 = n.zeros(nbins)
        for i in range(nbins):
            bins5[i] = (binedges5[i+1]+binedges5[i])/2.
        normhist5 = hist5 / float(count5)
        p.plot(bins5,normhist5,drawstyle='steps-mid', color=colors[4],
               label=labels[4])
        errors6 = self.dz_to_dv(z6,errors6)
        print('i', labels[5], n.mean(errors6), n.std(errors6))
        errors6 = n.log10(errors6)
        hist6,binedges6 = n.histogram(errors6, bins=nbins)
        bins6 = n.zeros(nbins)
        for i in range(nbins):
            bins6[i] = (binedges6[i+1]+binedges6[i])/2.
        normhist6 = hist6 / float(count6)
        p.plot(bins6,normhist6,drawstyle='steps-mid', color=colors[5],
               label=labels[5])
        p.text(0.8, 0.2, '$i$-band', fontsize=12)
        #p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        #p.title('r-band', size=18)
        p.axis([.5,2.5,0,.25])
        p.legend(prop={'size':6})
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print((count1+count2+count3+count4+count5+count6) / \
                float(self.rm_fibers_summary.shape[0]))

        ax3 = f.add_subplot(3,1,3)
        errors1 = n.array([])
        errors2 = n.array([])
        errors3 = n.array([])
        errors4 = n.array([])
        errors5 = n.array([])
        errors6 = n.array([])
        z1 = n.array([])
        z2 = n.array([])
        z3 = n.array([])
        z4 = n.array([])
        z5 = n.array([])
        z6 = n.array([])
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        openplate = 0
        openmjd = 0
        self.read_redmonster_summary_file()
        for j,sn_min in enumerate(n.linspace(1.5,4.0,6)):
            sn_max = sn_min + .5
            for i,fiber in enumerate(self.rm_fibers_summary):
                plate = self.rm_plates_summary[i]
                mjd = self.rm_mjds_summary[i]
                #print '%s-%s-%s' % (plate,fiber,mjd)
                if (openplate != plate) and (openmjd != mjd):
                    #self.read_spZbest_all(plate,mjd)
                    #self.read_spPlate_all(plate,mjd)
                    hduzbest = fits.open(join('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14', '%s' % plate, 'test_dr14', 'spZbest-%s-%s.fits' % (plate, mjd)))
                    self.sn_median = hduzbest[1].data.SN_MEDIAN[:,2:]
                    hduspplate = fits.open(join('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14', '%s' % plate, 'spPlate-%s-%s.fits' % (plate, mjd)))
                    openplate = plate
                    openmjd = mjd
                if (self.rm_zwarning[i] == 0) & (self.rm_zerr1[i] > 0):
                    if (self.sn_median[fiber][2] >= sn_min) & \
                            (self.sn_median[fiber][2] <= sn_max):
                        if j == 0:
                            errors1 = n.append(errors1,self.rm_zerr1[i])
                            z1 = n.append(z1,self.rm_z1[i])
                            count1 += 1
                        elif j == 1:
                            errors2 = n.append(errors2,self.rm_zerr1[i])
                            z2 = n.append(z2,self.rm_z1[i])
                            count2 += 1
                        elif j == 2:
                            errors3 = n.append(errors3,self.rm_zerr1[i])
                            z3 = n.append(z3,self.rm_z1[i])
                            count3 += 1
                        elif j == 3:
                            errors4 = n.append(errors4,self.rm_zerr1[i])
                            z4 = n.append(z4,self.rm_z1[i])
                            count4 += 1
                        elif j == 4:
                            errors5 = n.append(errors5,self.rm_zerr1[i])
                            z5 = n.append(z5,self.rm_z1[i])
                            count5 += 1
                        elif j == 5:
                            errors6 = n.append(errors6,self.rm_zerr1[i])
                            z6 = n.append(z6,self.rm_z1[i])
                            count6 += 1

        errors1 = self.dz_to_dv(z1,errors1)
        print('z', labels[0], n.mean(errors1), n.std(errors1))
        errors1 = n.log10(errors1)
        hist1,binedges1 = n.histogram(errors1, bins=nbins)
        bins1 = n.zeros(nbins)
        for i in range(nbins):
            bins1[i] = (binedges1[i+1]+binedges1[i])/2.
        normhist1 = hist1 / float(count1)
        p.plot(bins1,normhist1,drawstyle='steps-mid', color=colors[0],
               label=labels[1])
        errors2 = self.dz_to_dv(z2,errors2)
        print('z', labels[1], n.mean(errors2), n.std(errors2))
        errors2 = n.log10(errors2)
        hist2,binedges2 = n.histogram(errors2, bins=nbins)
        bins2 = n.zeros(nbins)
        for i in range(nbins):
            bins2[i] = (binedges2[i+1]+binedges2[i])/2.
        normhist2 = hist2 / float(count2)
        p.plot(bins2,normhist2,drawstyle='steps-mid', color=colors[1],
               label=labels[2])
        errors3 = self.dz_to_dv(z3,errors3)
        print('z', labels[2], n.mean(errors3), n.std(errors3))
        errors3 = n.log10(errors3)
        hist3,binedges3 = n.histogram(errors3, bins=nbins)
        bins3 = n.zeros(nbins)
        for i in range(nbins):
            bins3[i] = (binedges3[i+1]+binedges3[i])/2.
        normhist3 = hist3 / float(count3)
        p.plot(bins3,normhist3,drawstyle='steps-mid', color=colors[2],
               label=labels[3])
        errors4 = self.dz_to_dv(z4,errors4)
        print('z', labels[3], n.mean(errors4), n.std(errors4))
        errors4 = n.log10(errors4)
        hist4,binedges4 = n.histogram(errors4, bins=nbins)
        bins4 = n.zeros(nbins)
        for i in range(nbins):
            bins4[i] = (binedges4[i+1]+binedges4[i])/2.
        normhist4 = hist4 / float(count4)
        p.plot(bins4,normhist4,drawstyle='steps-mid', color=colors[3],
               label=labels[4])
        errors5 = self.dz_to_dv(z5,errors5)
        print('z', labels[4], n.mean(errors5), n.std(errors5))
        errors5 = n.log10(errors5)
        hist5,binedges5 = n.histogram(errors5, bins=nbins)
        bins5 = n.zeros(nbins)
        for i in range(nbins):
            bins5[i] = (binedges5[i+1]+binedges5[i])/2.
        normhist5 = hist5 / float(count5)
        p.plot(bins5,normhist5,drawstyle='steps-mid', color=colors[4],
               label=labels[5])
        errors6 = self.dz_to_dv(z6,errors6)
        print('z', labels[5], n.mean(errors6), n.std(errors6))
        errors6 = n.log10(errors6)
        hist6,binedges6 = n.histogram(errors6, bins=nbins)
        bins6 = n.zeros(nbins)
        for i in range(nbins):
            bins6[i] = (binedges6[i+1]+binedges6[i])/2.
        normhist6 = hist6 / float(count6)
        p.plot(bins6,normhist6,drawstyle='steps-mid', color=colors[5],
               label=labels[6])
        p.text(0.8, 0.28, '$z$-band', fontsize=12)
        p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        #p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        p.axis([.5,2.5,0,.35])
        p.legend(prop={'size':6})
        p.subplots_adjust(hspace = .5)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/dv_vs_sn_histos.pdf')
        p.clf()
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print((count1+count2+count3+count4+count5+count6) / \
                float(self.rm_fibers_summary.shape[0]))


    def sequels_failure_vs_imag_all(self,imin=18,imax=24,nbins=21):
        # Makes plot of SEQUELS LRG failure rate (zwarning > 0)
        # vs i-band magnitude
        f = p.figure()
        ax = f.add_subplot(1,1,1)
        total = 0
        bad_i_mag = []
        i_mag = []
        openplate = 0
        openmjd = 0
        self.read_redmonster_summary_file()
        for i,fiber in enumerate(self.rm_fibers_summary):
            plate = self.rm_plates_summary[i]
            mjd = self.rm_mjds_summary[i]
            print('%s-%s-%s' % (plate,mjd,fiber))
            if (openplate != plate) and (openmjd != mjd):
                self.read_spZbest_all(plate,mjd)
                self.read_spPlate_all(plate,mjd)
                openplate = plate
                openmjd = mjd
            if (self.spectroflux[fiber,3] <= imax):
                total += 1.
                i_mag.append(self.spectroflux[fiber,3])
                if (self.rm_zwarning[i] & 4 > 0):
                    bad_i_mag.append(self.spectroflux[fiber,3])
        nbinsarr = n.linspace(imin,imax,nbins+1)
        itotal,ibinedges = n.histogram(i_mag,bins=nbinsarr)
        ihist,ibinedges = n.histogram(bad_i_mag,bins=nbinsarr)
        ibins = n.zeros(nbins)
        for i in range(nbins):
            ibins[i] = (ibinedges[i+1]+ibinedges[i])/2.
        ihist = ihist / itotal.astype(float)
        for i in range(nbins):
            if i != 0 and i != (nbins-1):
                if isnan(ihist[i]):
                    try:
                        ihist[i] = (ihist[i-1] + ihist[i+1]) / 2.
                    except:
                        ihist[i] = 0
        p.plot(ibins,ihist,color='blue',drawstyle='steps-mid',label='i-band')
        p.axis([imin,imax,.01,1])
        ax.set_yscale('log')
        p.axvline(21.8,linestyle='--',color='k')
        p.xlabel(r'$i$-band magnitude',size=14)
        p.ylabel(r'Failure rate', size=14)
        #print rbins
        #print rhist
        #print rtotal
        #p.legend()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/failure_vs_imag.pdf')
        p.clf()


    def sequels_logdv_vs_z_scatter_all(self,nobjs=100000):
    # Makes a scatterplot nobjs CMASS targets of redshift vs log(dv)
        self.read_redmonster_summary_file()
        errors = []
        zs = []
        for i in range(nobjs):
            if (self.rm_zwarning[i] == 0) & \
                    (self.rm_type[i] == 'ssp_galaxy_glob') & \
                    (self.rm_zerr1[i] != -1):
                errors.append(self.rm_zerr1[i])
                zs.append(self.rm_z1[i])
        errors = self.dz_to_dv(n.asarray(errors))
        logerrs = n.log10(errors)
        p.scatter(zs,logerrs, marker='.')
        p.axhline(2.48,linestyle='--',color='k')
        p.xlabel('Redshift',size=16)
        p.ylabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/\
                  dv_vs_z_scatter.pdf')
        p.clf()


    def sequels_chi2_histos(self,nbins=50, rchi2=True):
        # Makes histogram of SEQUELS chi2 values for redmonster
        # and idlspec1d (in chi2 or rchi2)
        rm_rchi2s = []
        idl_rchi2s = []
        openplate = 0
        openmjd = 0
        total = 0
        self.read_redmonster_summary_file()
        for i,fiber in enumerate(self.rm_fibers_summary):
            plate = self.rm_plates_summary[i]
            mjd = self.rm_mjds_summary[i]
            print('%s-%s-%s' % (plate,mjd,fiber))
            if (openplate != plate) and (openmjd != mjd):
                self.read_spZbest_all(plate,mjd)
                self.read_spPlate_all(plate,mjd)
                openplate = plate
                openmjd = mjd
            if (self.rm_rchi2s[i] < 2) and (self.idl_rchi2s[fiber] < 2):
                total += 1
                if rchi2:
                    rm_rchi2s.append(self.rm_rchi2s[i])
                    idl_rchi2s.append(self.idl_rchi2s[fiber])
                else:
                    rm_rchi2s.append(self.rm_rchi2s[i] * self.rm_dof[i])
                    idl_rchi2s.append(self.idl_rchi2s[fiber] *
                                      self.idl_dof[fiber])
                #rm_rchi2s.append(self.rm_rchi2diff[i])
                #idl_rchi2s.append(self.idl_rchi2diff[fiber])
        rmhist,rmbinedges = n.histogram(rm_rchi2s,nbins)
        rmbins = n.zeros(nbins)
        for i in range(nbins):
            rmbins[i] = (rmbinedges[i+1]+rmbinedges[i])/2.
        rmhist = rmhist / float(total)
        idlhist, idlbinedges = n.histogram(idl_rchi2s,nbins)
        idlbins = n.zeros(nbins)
        for i in range(nbins):
            idlbins[i] = (idlbinedges[i+1]+idlbinedges[i])/2.
        idlhist = idlhist / float(total)
        p.plot(rmbins, rmhist, color='red', drawstyle='steps-mid',
               label='redmonster')
        p.plot(idlbins, idlhist, color='blue', drawstyle='steps-mid',
               label='idlspec1d')
        p.xlabel(r'$\chi_r^2$', size=16)
        p.ylabel(r'Fraction per bin', size=16)
        p.legend()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/rchi2_histos.pdf')
        p.clf()


    def sequels_drchi2_histos(self, drchi2max=.02, nbins=50, rchi2=True):
        # Makes histogram of SEQUELS delta-rchi2 (or delta-chi2) values
        # for redmonster and idlspec1d
        rm_drchi2s = []
        idl_drchi2s = []
        openplate = 0
        openmjd = 0
        total = 0
        self.read_redmonster_summary_file()
        for i,fiber in enumerate(self.rm_fibers_summary):
            plate = self.rm_plates_summary[i]
            mjd = self.rm_mjds_summary[i]
            print('%s-%s-%s' % (plate,mjd,fiber))
            if (openplate != plate) and (openmjd != mjd):
                self.read_spZbest_all(plate,mjd)
                self.read_spPlate_all(plate,mjd)
                openplate = plate
                openmjd = mjd
            if (self.rm_rchi2diff[i] < drchi2max) and \
                    (self.idl_rchi2diff[fiber] < drchi2max):
                total += 1
                if rchi2:
                    rm_drchi2s.append( self.rm_rchi2diff[i] )
                    idl_drchi2s.append( self.idl_rchi2diff[fiber] )
                else:
                    rm_drchi2s.append( self.rm_rchi2diff[i] * self.rm_dof[i] )
                    idl_drchi2s.append( self.idl_rchi2diff[fiber] *
                                       self.idl_dof[fiber] )
        rmhist, rmbinedges = n.histogram(rm_drchi2s,nbins)
        rmbins = n.zeros(nbins)
        for i in range(nbins):
            rmbins[i] = (rmbinedges[i+1]+rmbinedges[i])/2.
        rmhist = rmhist / float(total)
        idlhist, idlbinedges = n.histogram(idl_drchi2s,nbins)
        idlbins = n.zeros(nbins)
        for i in range(nbins):
            idlbins[i] = (idlbinedges[i+1]+idlbinedges[i])/2.
        idlhist = idlhist / float(total)
        p.plot(rmbins, rmhist, color='red', drawstyle='steps-mid',
               label='redmonster')
        p.plot(idlbins, idlhist, color='blue', drawstyle='steps-mid',
               label='idlspec1d')
        p.axvline(.005,linestyle='--',color='red')
        p.axvline(.01,linestyle='--',color='blue')
        p.xlabel(r'$\Delta\chi_r^2$', size=16)
        p.ylabel(r'Fraction per bin', size=16)
        p.legend()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/drchi2_histos.pdf')
        p.clf()


    def dchi2_failure_diff_function(self, diff, drchi2max=.05):
        # Helper function for sequels_failure_vs_dchi2()
        rm_failures = 0
        idl_failures = 0
        openplate = 0
        openmjd = 0
        total = 0
        self.read_redmonster_summary_file()
        for i,fiber in enumerate(self.rm_fibers_summary):
            plate = self.rm_plates_summary[i]
            mjd = self.rm_mjds_summary[i]
            #print '%s-%s-%s' % (plate,mjd,fiber)
            if (openplate != plate) and (openmjd != mjd):
                #self.read_spZbest_all(plate,mjd)
                #self.read_spPlate_all(plate,mjd)
                hduidl = fits.open(join(environ['BOSS_SPECTRO_REDUX'], 'test/bautista/test_dr14', '%s' % plate, 'test_dr14', 'spZbest-%s-%s.fits' % (plate, mjd)))
                self.idl_rchi2diff = hduidl[1].data.RCHI2DIFF_NOQSO
                self.idl_dof = hduidl[1].data.DOF
                openplate = plate
                openmjd = mjd
            #if (self.rm_rchi2diff[i] < drchi2max) and \
                    #(self.idl_rchi2diff[fiber] < drchi2max):
            total += 1
            if self.rm_rchi2diff[i] < diff: rm_failures += 1.
            if self.idl_rchi2diff[fiber] < diff: idl_failures += 1.
        return (rm_failures/total), (idl_failures/total)


    def sequels_failure_vs_dchi2(self, drchi2max=.02, npoints=150):
        # Makes a plot of SEQUELS LRG failure rate as a function of
        # dchi2 threshold for redmonster and idlspec1d
        rm_data = []
        idl_data = []
        diffs = n.linspace(0,drchi2max,npoints)
        for i,diff in enumerate(diffs):
            print('%s of %s' % (i+1,npoints))
            rm_point, idl_point = self.dchi2_failure_diff_function(diff)
            rm_data.append(rm_point)
            idl_data.append(idl_point)
        p.plot(diffs, rm_data, 'red', label='redmonster')
        p.plot(diffs, idl_data, 'blue', label='idlspec1d')
        p.xlabel(r'$\Delta\chi_{r}^2$ threshold', size=16)
        p.ylabel(r'Cumulative fraction below threshold', size=16)
        p.grid(b=True, which='major', color='black', linestyle='--')
        p.legend(loc=2)
        p.axis([0,.02,0,.7])
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/drchi2_vs_failure.pdf')
        p.clf()


    def sequels_reobs_errors(self, nbins=25):
        # Makes a histogram of (z2-z1)/sqrt(dz1**2 + dz2**2)
        # with best fit Gaussian overplotted for all SEQUELS LRG
        # targets with repeat observations
        globpath = join( self.redmonster_spectro_redux,'*')
        z1 = []
        z2 = []
        zerr1 = []
        zerr2 = []
        for path in iglob(globpath):
            plate = basename(path)
            if plate != 'redmonsterAll-%s.fits' % self.version:
                print(plate)
                mjds = []
                mjdglobpath = join( self.redmonster_spectro_redux, '%s' % plate,
                                   '%s' % self.version,
                                   'redmonster-%s-*.fits' % plate)
                for mjdpath in iglob(mjdglobpath):
                    mjd = basename(mjdpath)[16:21]
                    if mjd not in mjds:
                        mjds.append(mjd)
                if len(mjds) > 1:
                    print('Plate %s has multiple MJDs' % plate)
                    hdu1 = fits.open( join( self.redmonster_spectro_redux,
                                           plate, self.version,
                                           'redmonster-%s-%s.fits' %
                                           (plate,mjds[0]) ) )
                    hdu2 = fits.open( join( self.redmonster_spectro_redux,
                                           plate, self.version,
                                           'redmonster-%s-%s.fits' %
                                           (plate,mjds[1]) ) )
                    for i,z in enumerate(hdu1[1].data.Z1):
                        if (hdu1[1].data.ZWARNING[i] == 0) & \
                                (hdu1[1].data.CLASS1[i] == 'ssp_galaxy_glob') & \
                                (hdu2[1].data.ZWARNING[i] == 0) & \
                                (hdu2[1].data.CLASS1[i] == 'ssp_galaxy_glob'):
                            z1.append(z)
                            z2.append(hdu2[1].data.Z1[i])
                            zerr1.append(hdu1[1].data.Z_ERR1[i])
                            zerr2.append(hdu2[1].data.Z_ERR1[i])
        z1 = n.array(z1)
        z2 = n.array(z2)
        zerr1 = n.array(zerr1)
        zerr2 = n.array(zerr2)
        z_diff = z2-z1
        zerr_rms = n.sqrt( (zerr1**2 + zerr2**2) )
        scaled_diff = z_diff / zerr_rms
        hist,binedges = n.histogram(scaled_diff,bins=nbins)
        normhist = hist / float(z1.shape[0])
        bins = n.zeros(nbins)
        for i in range(nbins):
            bins[i] = (binedges[i+1]+binedges[i])/2.
        p.plot(bins, hist, drawstyle='steps-mid', color='black')

        def fit_func(x,a,sigma,mu):
            # Gaussian function to fit to histogram
            return a * n.exp( -((x-mu)**2)/(2*sigma**2) )

        popt,pcov = curve_fit(fit_func, normhist,bins)
        xfit = n.linspace(-6,6,1000)
        yfit = fit_func(xfit, popt[0], popt[1], popt[2])
        p.plot(xfit,yfit,color='cyan')
        p.xlabel(r'$(z_2-z_1)/ \delta z_{rms}$', size=16)
        p.ylabel('Fraction per bin',size=16)
        p.text(3,.01,r'$\sigma_{fit}=1.18$',size=18)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/reobs_errors.pdf')
        p.clf()


    def plate_splits_function(self, plate, mjd, nbins=25, fit=True):
        hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                              'test/bautista/v5_8_guy_split1', '%s' % plate,
                              'v5_8_guy_split1',
                              'redmonster-%s-%s.fits' % (plate,mjd)))
        hdu2 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                              'test/bautista/v5_8_guy_split2', '%s' % plate,
                              'v5_8_guy_split2',
                              'redmonster-%s-%s.fits' % (plate,mjd)))
        hdu3 = fits.open(join(environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'],
                              '%s' % plate, 'spPlate-%s-%s.fits' % (plate,mjd)))
        for i,ebt1 in enumerate(hdu3[5].data.EBOSS_TARGET1):
            if ebt1 & 2 > 0:
                if True: #(hdu1[1].data.ZWARNING[i] == 0) and \
                                #(hdu2[1].data.ZWARNING[i] == 0):
                    self.z1.append(hdu1[1].data.Z1[i])
                    self.zerr1.append(hdu1[1].data.Z_ERR1[i])
                    self.z2.append(hdu2[1].data.Z1[i])
                    self.zerr2.append(hdu2[1].data.Z_ERR1[i])
                    if n.abs(self.z1[-1] - self.z2[-1]) > .01:
                        del(self.z1[-1])
                        del(self.z2[-1])
                        del(self.zerr1[-1])
                        del(self.zerr2[-1])
                    if n.sqrt( (self.zerr1[-1]**2 + self.zerr2[-1]**2) ) == 0:
                        del(self.z1[-1])
                        del(self.z2[-1])
                        del(self.zerr1[-1])
                        del(self.zerr2[-1])


    def plate_splits_errors(self, nbins=25, fit=True, normed=True):
        plates = [7834,7839,7848]
        mjds = [56979,56900,56959]
        self.z1 = []
        self.z2 = []
        self.zerr1 = []
        self.zerr2 = []
        for i,plate in enumerate(plates):
            self.plate_splits_function(plate=plate, mjd=mjds[i],
                                       nbins=nbins, fit=fit)
        self.z1 = n.array(self.z1)
        self.z2 = n.array(self.z2)
        self.zerr1 = n.array(self.zerr1)
        self.zerr2 = n.array(self.zerr2)
        z_diff = self.z2-self.z1
        zerr_rms = n.sqrt( (self.zerr1**2 + self.zerr2**2) )
        scaled_diff = z_diff / zerr_rms
        while True:
            if n.abs(scaled_diff[n.abs(scaled_diff).argmax()]) > 5:
                scaled_diff = n.delete(scaled_diff, n.abs(scaled_diff).argmax())
            else:
                break
        print(n.max(n.abs(scaled_diff)))
        print(scaled_diff.shape)
        hist,binedges = n.histogram(scaled_diff, bins = nbins)
        if normed:
            normhist = hist / float(self.z1.shape[0])
        else:
            normhist = hist
        bins = n.zeros(nbins)
        for i in range(nbins):
            bins[i] = (binedges[i+1]+binedges[i])/2.
        p.plot(bins, normhist, drawstyle='steps-mid', color='black')

        def fit_func(x, a, sigma, mu): # Gaussian function to fit to histogram
            return a * n.exp( -((x-mu)**2)/(2.*sigma**2) )

        if fit:
            popt, pcov = curve_fit(fit_func, bins, normhist)
            xfit = n.linspace(-4,4,1000)
            yfit = fit_func(xfit, popt[0], popt[1], popt[2])
            p.plot(xfit, yfit, color='mediumpurple')
            p.text(2.4,.11, r'$\sigma_{\mathrm{fit}}=$%.2f' % popt[1], size=16)
            p.text(2.4, .10, r'$\mu_{\mathrm{fit}}=$%.2f' % popt[2], size=16)

        p.xlabel(r'$(z_2-z_1)/ (\delta z_1^2+$ $\delta z_2^2)^{1/2}$', size=16)
        p.ylabel('Fraction per bin', size=16)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/reobs_errors.pdf')
        p.clf()


    def sequels_sky_drchi2(self):
        xdata = n.linspace(0,.01,40)
        rm_ydata = []
        idl_ydata = []
        globpath1 = join(environ['REDMONSTER_SPECTRO_REDUX'],self.version, '*')
        for chi2max in xdata:
            print(chi2max)
            total = 0.
            countidl = 0.
            countrm = 0.
            for path in iglob(globpath1):
                plate = basename(path)
                if len(plate) == 4:
                    globpath2 = join(environ['BOSS_SPECTRO_REDUX'],
                                     '%s' % self.version, '%s' % plate,
                                     'spPlate-%s-*.fits' % plate)
                    for file in iglob(globpath2):
                        try:
                            if len(basename(file)) == 23:
                                mjd = basename(file)[13:18]
                                hduplate = fits.open(file)
                                hduidl=fits.open(join(environ['BOSS_SPECTRO_REDUX'],
                                                      '%s' % self.version,
                                                      '%s' % plate,
                                                      '%s' % self.version,
                                                      'spZbest-%s-%s.fits' %
                                                      (plate,mjd)))
                                hdurm = fits.open(
                                        join(environ['REDMONSTER_SPECTRO_REDUX'],
                                             '%s' % self.version, '%s' % plate,
                                             '%s' % self.version,
                                             'redmonster-%s-%s.fits' % (plate,mjd)))
                                for i,zwarn in enumerate(hdurm[1].data.ZWARNING):
                                    if zwarn & 1 > 0:
                                        total += 1.
                                        if hduidl[1].data.RCHI2DIFF[i] > chi2max:
                                            countidl += 1.
                                        if hdurm[1].data.RCHI2DIFF[i] > chi2max:
                                            countrm += 1.
                        except IOError: pass
            rm_ydata.append(countrm/total)
            idl_ydata.append(countidl/total)

        f = p.figure()
        ax = f.add_subplot(1,1,1)
        p.plot(xdata, rm_ydata, drawstyle='steps-mid', color='red',
               label='redmonster')
        p.plot(xdata, idl_ydata, drawstyle='steps-mid', color='blue',
               label='idlspec1d')
        p.xlabel(r'$\Delta\chi_r^2$', size=16)
        p.ylabel(r'Cumulative fraction above threshold', size=16)
        ax.set_yscale('log')
        p.grid(b=True, which='major', color='black', linestyle='-')
        p.grid(b=True, which='minor', color='black', linestyle='--')
        p.legend()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/sky_failure_vs_drchi2.pdf')
        p.clf()


    def sequels_failure_confusion(self):
        galgal = 0
        galstar = 0
        galqso = 0
        starstar = 0
        starqso = 0
        qsoqso = 0
        total = 0
        self.read_redmonster_summary_file()
        openplate = 0
        openmjd = 0
        for i,zwarn in enumerate(self.rm_zwarning):
            if zwarn & 4 > 0:
                total += 1.
                plate = self.rm_plates_summary[i]
                mjd = self.rm_mjds_summary[i]
                fiber = self.rm_fibers_summary[i]
                if openplate != plate or openmjd != mjd:
                    self.read_redmonster_all(plate,mjd)
                    openplate = plate
                    openmjd = mjd
                ind = n.where(self.rm_fibers == fiber)[0][0]
                if self.rm_type[ind][:3] == 'ssp':
                    if self.rm_type2[ind][:3] == 'ssp': galgal += 1
                    elif self.rm_type2[ind][:3] == 'QSO': galqso += 1
                    elif self.rm_type2[ind][:3] == 'CAP': galstar += 1
                if self.rm_type[ind][:3] == 'QSO':
                    if self.rm_type2[ind][:3] == 'ssp': galqso += 1
                    elif self.rm_type2[ind][:3] == 'QSO': qsoqso += 1
                    elif self.rm_type2[ind][:3] == 'CAP': starqso += 1
                if self.rm_type[ind][:3] == 'CAP':
                    if self.rm_type2[ind][:3] == 'ssp': galstar += 1
                    if self.rm_type2[ind][:3] == 'QSO': starqso += 1
                    if self.rm_type2[ind][:3] == 'CAP': starstar += 1
        print('%s galaxy-galaxy confusions of %s, which is %s' % \
                (galgal,total,(galgal/total)*100))
        print('%s galaxy-star confusions of %s, which is %s' % \
                (galstar,total,(galstar/total)*100))
        print('%s galaxy-QSO confusions of %s, which is %s' % \
                (galqso,total,(galqso/total)*100))
        print('%s star-star confusions of %s, which is %s' % \
                (starstar,total,(starstar/total)*100))
        print('%s star-QSO confusions of %s, which is %s' % \
                (starqso,total,(starqso/total)*100))
        print('%s QSO-QSO confusions of %s, which is %s' % \
                (qsoqso,total,(qsoqso/total)*100))


    def rchi2_null_histos(self, nbins=35, reduced=True, normed=True):
        self.read_redmonster_summary_file()
        if reduced:
            rchi2_nulls = self.rm_chi2_null / self.rm_dof
            rchi2_nulls = rchi2_nulls[n.where(rchi2_nulls < 1.8)[0]]
            rchi2_nulls = rchi2_nulls[n.where(rchi2_nulls > 0.6)[0]]
            xdata = n.linspace(.6,1.8,400)
        else:
            rchi2_nulls = rchi2_nulls[n.where(rchi2_nulls < 8000)[0]]
            rchi2_nulls = rchi2_nulls[n.where(rchi2_nulls > 3000)[0]]
            xdata = n.linspace(3000,8000,400)
        # Plot normal histogram
        hist, binedges = n.histogram(rchi2_nulls, bins=nbins)
        normhist = hist / float(rchi2_nulls.shape[0])
        bins = n.zeros(nbins)
        for i in range(nbins):
            bins[i] = (binedges[i+1] + binedges[i]) / 2.
        if normed:
            p.plot(bins, normhist, drawstyle='steps-mid')
            p.ylabel('Fraction per bin')
        else:
            p.plot(bins, hist, drawstyle='steps-mid')
            p.yabel('Number per bin')
        if reduced:
            p.xlabel(r'$\chi_{null,red}^2$',size=16)
            p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/rchi2_null_histo.pdf')
        else:
            p.xlabel(r'$\chi_{null}^2$',size=16)
            p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/chi2_null_histo.pdf')
        p.clf()
        # Plot cumulative histogram
        ydata = []
        for xpoint in xdata:
            ypoint = 0
            for rchi2null in rchi2_nulls:
                if rchi2null < xpoint:
                    ypoint += 1
            ydata.append( ypoint/float(rchi2_nulls.shape[0]) )
        p.plot(xdata,ydata)
        p.ylabel('Cumulative fraction below', size=16)
        p.grid(b=True, which='major', color='black', linestyle='--')
        if reduced:
            p.xlabel(r'$\chi_{null,red}^2$',size=16)
            p.axis([.6,1.8,0,1])
            p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/rchi2_null_cumul_histo.pdf')
        else:
            p.xlabel(r'$\chi_{null}^2$',size=16)
            p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/chi2_null_cumul_histo.pdf')
        p.clf()


    def chi2_null_less_chi2_min(self, nbins=35, normed=True):
        self.read_redmonster_summary_file()
        diffs = self.rm_chi2_null - (self.rm_rchi2s * self.rm_dof)
        diffs = diffs[n.where(diffs < 1000)[0]]
        diffs = diffs[n.where(diffs > -1000)[0]]
        hist, binedges = n.histogram(diffs, bins=nbins)
        normhist = hist / float(diffs.shape[0])
        bins = n.zeros(nbins)
        diffs2 = self.rm_rchi2diff * self.rm_dof
        diffs2 = diffs2[n.where(diffs2 < 1000)[0]]
        diffs2 = diffs2[n.where(diffs2 > -1000)[0]]
        hist2, binedges2 = n.histogram(diffs2,bins=nbins)
        normhist2 = hist2 / float(diffs2.shape[0])
        bins2 = n.zeros(nbins)
        diffs3 = self.rm_sn2_data - (self.rm_rchi2s * self.rm_dof)
        diffs3 = diffs3[n.where(diffs3 < 1000)[0]]
        diffs3 = diffs3[n.where(diffs3 > -1000)[0]]
        hist3, binedges3 = n.histogram(diffs3,bins=nbins)
        normhist3 = hist3 / float(diffs3.shape[0])
        bins3 = n.zeros(nbins)
        for i in range(nbins):
            bins[i] = (binedges[i+1] + binedges[i]) / 2.
            bins2[i] = (binedges2[i+1] + binedges2[i]) / 2.
            bins3[i] = (binedges3[i+1] + binedges3[i]) / 2.
        if normed:
            p.plot(bins, normhist, drawstyle='steps-mid', color='magenta',
                   label=r'$\chi_{\mathrm{null}}^2$')
            p.plot(bins2, normhist2, drawstyle='steps-mid',
                   color='mediumpurple', label=r'$\chi_{\mathrm{fit2}}^2$')
            p.plot(bins3, normhist3, drawstyle='steps-mid', color='cyan',
                   label=r'$\chi_{0}^2$')
            p.ylabel('Fraction per bin')
        else:
            p.plot(bins, hist, drawstyle='steps-mid', color='magenta',
                   label=r'$\chi_{\mathrm{null}}^2$')
            p.plot(bins2, hist2, drawstyle='steps-mid', color='mediumpurple',
                   label=r'$\chi_{\mathrm{fit2}}^2$')
            p.plot(bins3, hist3, drawstyle='steps-mid', color='cyan',
                   label=r'$\chi_{0}^2$')
            p.ylabel('Number per bin')
        if bins[0] < bins2[0]:
            if bins[-1] > bins2[-1]: p.axis([bins[0],bins[-1],0,.25])
            else: p.axis([bins[0],bins2[-1],0,.25])
        else:
            if bins[-1] > bins2[-1]: p.axis([bins2[0],bins[-1],0,.25])
            else: p.axis([bins2[0],bins2[-1],0,.25])
        p.xlabel(r'$\chi^2-\chi_{\mathrm{min}}^2$', size=16)
        p.legend()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/chi2_null_less_chi2_min_histo.pdf')
        p.clf()
        print(bins[0:3], normhist[0:3], bins2[0:3], normhist2[0:3], bins3[0:3],\
                normhist3[0:3], end=' ')


    def sequels_stack_spectra(self):
        pass


    def sequels_example_chi2s(self, plates, mjds, fibers):
        # Create stacked plots of three chi2 vs z curves for three
        # fibers.  plates, mjds, and fibers are lists
        # of plate, mjd, fiberid sets to have chi2 curve plotted
        # 7397 57129 784
        # 7311 57038 465
        # 7305 56991 692
        # This is a comment
        f = p.figure()
        ax1 = f.add_subplot(311)
        specs = spec.Spec(plate=plates[0], mjd=mjds[0], fiberid=[fibers[0]])
        zssp1 = zfinder.ZFinder(fname='ndArch-ssp_galaxy_glob-v000.fits',
                                npoly=4, zmin=-0.01, zmax=1.2)
        zssp1.zchi2(specs.flux, specs.loglambda, specs.ivar, npixstep=2)
        bestzvec = self.chi2_curves_helper(zssp1.zchi2arr, zssp1.zbase)
        p.plot(zssp1.zbase, [max(bestzvec)]*zssp1.zbase.shape[0], '--',
               color='mediumpurple')
        #p.plot(zssp1.zbase, [zssp1.sn2_data[0]]*zssp1.zbase.shape[0], '--',
        #       color='mediumaquamarine')
        #print zssp1.sn2_data[0]
        p.plot(zssp1.zbase, bestzvec, color='black')
        p.text(0.1,3970,'7397-57129-785', fontsize=12)
        p.text(0.89, 3970, r'$\chi_0^2 = $ %.1f' % zssp1.sn2_data[0],
               fontsize=12)
        p.axis([0,1.2,3955,4110])
        ax1.set_yticks([3960,3990,4020,4050,4080,4110])
        ax1.set_xticks([0,0.2,0.4,0.6,0.8,1.0,1.2])
        ax2 = f.add_subplot(312)
        specs = spec.Spec(plate=plates[1], mjd=mjds[1], fiberid=[fibers[1]])
        zssp1 = zfinder.ZFinder(fname='ndArch-ssp_galaxy_glob-v000.fits',
                                npoly=4, zmin=-0.01, zmax=1.2)
        zssp1.zchi2(specs.flux, specs.loglambda, specs.ivar, npixstep=2)
        bestzvec = self.chi2_curves_helper(zssp1.zchi2arr, zssp1.zbase)
        p.plot(zssp1.zbase, [max(bestzvec)]*zssp1.zbase.shape[0], '--',
               color='mediumpurple')
        #p.plot(zssp1.zbase, [zssp1.sn2_data[0]]*zssp1.zbase.shape[0], '--',
                #color='mediumaquamarine')
        #print zssp1.sn2_data[0]
        p.plot(zssp1.zbase, bestzvec, color='black')
        p.text(0.1,4865.80645,'7311-57038-466', fontsize=12)
        p.text(0.89, 4865.80645, r'$\chi_0^2 = $ %.1f' % zssp1.sn2_data[0],
               fontsize=12)
        p.axis([0,1.2,4860,4920])
        ax2.set_yticks([4860,4875,4890,4905,4920])
        ax2.set_xticks([0,0.2,0.4,0.6,0.8,1.0,1.2])
        p.ylabel(r'$\chi^2$', size=16)
        ax3 = f.add_subplot(313)
        specs = spec.Spec(plate=plates[2], mjd=mjds[2], fiberid=[fibers[2]])
        zssp1 = zfinder.ZFinder(fname='ndArch-ssp_galaxy_glob-v000.fits',
                                npoly=4, zmin=-0.01, zmax=1.2)
        zssp1.zchi2(specs.flux, specs.loglambda, specs.ivar, npixstep=2)
        bestzvec = self.chi2_curves_helper(zssp1.zchi2arr, zssp1.zbase)
        p.plot(zssp1.zbase, [max(bestzvec)]*zssp1.zbase.shape[0], '--',
               color='mediumpurple')
        #p.plot(zssp1.zbase, [zssp1.sn2_data[0]]*zssp1.zbase.shape[0], '--',
                #color='mediumaquamarine')
        #print zssp1.sn2_data[0]
        p.plot(zssp1.zbase, bestzvec, color='black')
        p.text(0.1,4001.6451612,'7305-56991-693', fontsize=12)
        p.text(0.89, 4001.6451612, r'$\chi_0^2 = $ %.1f' % zssp1.sn2_data[0],
               fontsize=12)
        p.axis([0,1.2,4000,4017])
        ax3.set_yticks([4000,4004,4008,4012,4016])
        ax3.set_xticks([0,0.2,0.4,0.6,0.8,1.0,1.2])
        p.xlabel(r'$z$', size=16)
        p.subplots_adjust(wspace = .3, hspace = .3)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/example_chi2_vs_z.pdf')
        p.clf()


    def sequels_1poly_vs_4poly_scatters(self):
        hdu1poly = fits.open( join(self.redmonster_spectro_redux + '_poly1',
                                   'redmonsterAll-%s.fits' % self.version) )
        hdu4poly = fits.open( join(self.redmonster_spectro_redux + '_poly4',
                                   'redmonsterAll-%s.fits' % self.version) )
        # Chi2null
        yes1yes4 = []
        yes1no4 = []
        no1yes4 = []
        no1no4 = []
        openplate = None
        openmjd = None
        for i,zwarn in enumerate(hdu1poly[1].data.ZWARNING):
            stderr.write('\r %s of %s' % (i+1, hdu1poly[1].data.ZWARNING.shape[0]))
            plate = hdu1poly[1].data.PLATE[i]
            mjd = hdu1poly[1].data.MJD[i]
            fiberid = hdu1poly[1].data.FIBERID[i]
            if openplate != plate or openmjd != mjd:
                hdu4poly = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, '%s' % plate, self.version, 'redmonster-%s-%s.fits' % (plate,mjd)))
                openplate = plate
                openmjd = mjd
            try:
                fiberind = n.where(hdu4poly[1].data.FIBERID == fiberid)[0][0]
                thesechi2 = (hdu1poly[1].data.CHI2NULL[i], hdu4poly[1].data.CHI2NULL[fiberind])
                if not zwarn & 4:
                    if not hdu4poly[1].data.ZWARNING[fiberind] & 4:
                        yes1yes4.append(thesechi2)
                    else:
                        yes1no4.append(thesechi2)
                else:
                    if not hdu4poly[1].data.ZWARNING[fiberind] & 4:
                        no1yes4.append(thesechi2)
                    else:
                        no1no4.append(thesechi2)
            except IndexError:
                pass
        f = p.figure()
        ax1 = f.add_subplot(311)
        colors = ['black', 'tomato', 'darkturquoise', 'green']
        labels = ['Both', '1 poly', '4 poly', 'Neither']
        chi2list = [yes1yes4, yes1no4, no1yes4, no1no4]
        for i in range(3):
            x = []
            y = []
            for j in range(len(chi2list[i])):
                x.append(chi2list[i][j][0])
                y.append(chi2list[i][j][1])
            if i == 0: p.scatter(x, y, s=1, color=colors[i], label=labels[i], alpha=0.6) # lower alpha for grey points
            else: p.scatter(x, y, s=1, color=colors[i], label=labels[i], alpha=1)
        p.legend(loc=2, prop={'size':8})
        p.plot(n.linspace(0,10000,10000), n.linspace(0,10000,10000), color='black', linestyle='--')
        p.axis([2800,20000,3000,7000])
        ax1.set_yticks([3000,4000,5000,6000,7000])
        p.xlabel(r'$\chi_{\mathrm{null},1}^2$', size=12)
        p.ylabel(r'$\chi_{\mathrm{null},4}^2$', size=12)
        ypoints = []
        xpoints = []
        for pair in yes1yes4:
            xpoints.append(pair[0])
            ypoints.append(pair[1])
        xrms = n.sqrt(n.mean(n.square(xpoints)))
        yrms = n.sqrt(n.mean(n.square(ypoints)))

        # minrchi2
        yes1yes4 = []
        yes1no4 = []
        no1yes4 = []
        no1no4 = []
        for i,zwarn in enumerate(hdu1poly[1].data.ZWARNING):
            stderr.write('\r %s of %s' % (i+1, hdu1poly[1].data.ZWARNING.shape[0]))
            plate = hdu1poly[1].data.PLATE[i]
            mjd = hdu1poly[1].data.MJD[i]
            fiberid = hdu1poly[1].data.FIBERID[i]
            if openplate != plate or openmjd != mjd:
                hdu4poly = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, '%s' % plate, self.version, 'redmonster-%s-%s.fits' % (plate,mjd)))
                openplate = plate
                openmjd = mjd
            try:
                fiberind = n.where(hdu4poly[1].data.FIBERID == fiberid)[0][0]
                thesechi2 = (hdu1poly[1].data.MINRCHI2[i], hdu4poly[1].data.MINRCHI21[fiberind])
                if not zwarn & 4:
                    if not hdu4poly[1].data.ZWARNING[fiberind] & 4:
                        yes1yes4.append(thesechi2)
                    else:
                        yes1no4.append(thesechi2)
                else:
                    if not hdu4poly[1].data.ZWARNING[fiberind] & 4:
                        no1yes4.append(thesechi2)
                    else:
                        no1no4.append(thesechi2)
            except IndexError:
                pass
        f.add_subplot(312)
        colors = ['black', 'tomato', 'darkturquoise', 'green']
        labels = ['Both', '1 poly', '4 poly', 'Neither']
        chi2list = [yes1yes4, yes1no4, no1yes4, no1no4]
        for i in range(3):
            x = []
            y = []
            for j in range(len(chi2list[i])):
                x.append(chi2list[i][j][0])
                y.append(chi2list[i][j][1])
            if i == 0: p.scatter(x, y, s=1, color=colors[i], label=labels[i], alpha=0.6)
            else: p.scatter(x, y, s=1, color=colors[i], label=labels[i], alpha=1)
            p.legend(loc=2, prop={'size':8})
        p.plot(n.linspace(0,2,10000), n.linspace(0,2,10000), color='black', linestyle='--')
        p.axis([0.75,1.4,0.7,1.4])
        p.xlabel(r'$\chi_{\mathrm{r,min},1}^2$',size=12)
        p.ylabel(r'$\chi_{\mathrm{r,min},4}^2$',size=12)

        # rchi2diff
        yes1yes4 = []
        yes1no4 = []
        no1yes4 = []
        no1no4 = []
        for i,zwarn in enumerate(hdu1poly[1].data.ZWARNING):
            stderr.write('\r %s of %s' % (i+1, hdu1poly[1].data.ZWARNING.shape[0]))
            plate = hdu1poly[1].data.PLATE[i]
            mjd = hdu1poly[1].data.MJD[i]
            fiberid = hdu1poly[1].data.FIBERID[i]
            if openplate != plate or openmjd != mjd:
                hdu4poly = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, '%s' % plate, self.version, 'redmonster-%s-%s.fits' % (plate,mjd)))
                openplate = plate
                openmjd = mjd
            try:
                fiberind = n.where(hdu4poly[1].data.FIBERID == fiberid)[0][0]
                thesechi2 = (hdu1poly[1].data.RCHI2DIFF[i], hdu4poly[1].data.RCHI2DIFF[fiberind])
                if not zwarn & 4:
                    if not hdu4poly[1].data.ZWARNING[fiberind] & 4:
                        yes1yes4.append(thesechi2)
                    else:
                        #print '1:  %s %s %s' % (hdu1poly[1].data.PLATE[i],hdu1poly[1].data.MJD[i],hdu1poly[1].data.FIBERID[i])
                        yes1no4.append(thesechi2)
                else:
                    if not hdu4poly[1].data.ZWARNING[fiberind] & 4:
                        no1yes4.append(thesechi2)
                        #print '4:  %s %s %s' % (hdu1poly[1].data.PLATE[i],hdu1poly[1].data.MJD[i],hdu1poly[1].data.FIBERID[i])
                    else:
                        no1no4.append(thesechi2)
            except IndexError:
                pass
        f.add_subplot(313)
        colors = ['black', 'tomato', 'darkturquoise', 'green']
        labels = ['Both', '1 poly', '4 poly', 'Neither']
        chi2list = [yes1yes4, yes1no4, no1yes4, no1no4]
        for i in range(3):
            x = []
            y = []
            for j in range(len(chi2list[i])):
                x.append(chi2list[i][j][0])
                y.append(chi2list[i][j][1])
            if i == 0: p.scatter(x, y, s=1, color=colors[i], label=labels[i], alpha=0.6)
            else: p.scatter(x, y, s=1, color=colors[i], label=labels[i], alpha=1)
            p.legend(loc=2, prop={'size':8})
        p.axis([-0.008,0.05,-0.003,0.05])
        p.xlabel(r'$\Delta\chi_{\mathrm{r},1}^2$',size=12)
        p.ylabel(r'$\Delta\chi_{\mathrm{r},4}^2$',size=12)
        p.plot(n.linspace(-0.1,.1,10000), n.linspace(-0.1,.1,10000), color='black', linestyle='--')

        p.subplots_adjust(hspace = .8)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/1poly_4poly_scatters.pdf')
        p.clf()

        print('X_rms:  %s' % xrms)
        print('Y_rms:  %s' % yrms)

        '''
        import seaborn as sns
        g = (sns.jointplot(n.asarray(x), n.asarray(y), kind='reg').set_axis_labels('x', 'y'))
        g.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/test.pdf')
        '''




# ------------------------------------------------------------------------------
# modified version of zfitter.z_refine2() to create chi2 vs z curves
# for self.sequels_example_chi2s()

    def chi2_curves_helper(self, zchi2, zbase, threshold=23.3, width=15):
        self.zchi2 = zchi2
        self.zbase = zbase
        self.z = n.zeros((zchi2.shape[0],5))
        self.z_err = n.zeros((zchi2.shape[0],5))
        self.minvector = []
        self.zwarning = n.zeros(zchi2.shape[0])
        self.threshold = threshold
        self.width = width
        for ifiber in range(self.zchi2.shape[0]):
            self.minvector.append( (ifiber,) + \
                    n.unravel_index(self.zchi2[ifiber].argmin(),
                                    self.zchi2[ifiber].shape))
            bestzvec = n.zeros( self.zchi2.shape[-1])
            for iz in range(self.zchi2.shape[-1]):
                bestzvec[iz] = n.min( self.zchi2[ifiber,...,iz] )
        return bestzvec

# ------------------------------------------------------------------------------












# S/N per fiber is located in spZbest files in hdu[1].data.SN_MEDIAN .
# You can get just r,i,z bands with x = hdu[1].data.SN_MEDIAN[:,2:] .

# Fiber magnitudes are in spZbest files in hdu[1].data.SPECTROFLUX .
# Units are nanomaggies, convert to magnitudes
# with 22.5 - 2.5 * LOG_10(SPECTROFLUX)

# To see fibers with zwarning != 0, ztype = 'galaxy', and
# boss_target1 = 'cmass', use >>> print n.where( (x.rm_zwarning != 0) &
# (x.rm_type == 'ssp_galaxy_glob') & (x.boss_target1 & 2 == 2) )[0]+1

# Plate 7338 has 6 MJDs, 7340 has 4




# ------------------------------------------------------------------------------


# Below here are re-writes of plotting functions using seaborn



    def sequels_logdv_vs_z_histos_all_sns(self, nbins=12, sns_pal='deep'):
    # Make histograms of log10(dv) in redshift bins for
    # LOWZ and CMASS galaxies
        colors = [
              'tomato','sage','cornflowerblue','sandybrown',
              'mediumpurple','grey'
              ]
        labels = ['0.1<z<0.2','0.2<z<0.3','0.3<z<0.4','0.4<z<0.5']
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')
        f = p.figure()
        '''
            ax1 = f.add_subplot(1,2,1)
            for j,zmin in enumerate(n.linspace(.1,.4,4)):
            zmax = zmin + .1
            errors = n.array([])
            count = 0
            for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            self.get_all_yanny(plate)
            fibers = self.get_okay_lowz()
            fibers = self.redshift_bin_fibers(fibers, zmin, zmax)
            count += len(fibers)
            errors = n.append(errors, self.rm_zerr1[fibers])
            errors = self.dz_to_dv(errors)
            errors = n.log10(errors)
            hist,binedges = n.histogram(errors, bins=nbins)
            bins = n.zeros(nbins)
            for i in xrange(nbins):
            bins[i] = (binedges[i+1]+binedges[i])/2.
            normhist = hist / float(count)
            p.plot(bins,normhist,drawstyle='steps-mid', color=colors[j],
            label=labels[j])
            p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
            p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
            p.title('LOWZ Sample', size=18)
            p.legend()
            p.axis([.55,2,0,.4])
            '''
        ax2 = f.add_subplot(1,1,1)
        labels = ['0.6<z<0.7','0.7<z<0.8','0.8<z<0.9','0.9<z<1.0']
        nbins = 25
        for j,zmin in enumerate(n.linspace(.6,.9,4)):
            #import pdb; pdb.set_trace()
            zmax = zmin + .1
            errors = n.array([])
            zs = n.array([])
            count = 0
            '''
                for plate in self.plates:
                self.read_redmonster(plate)
                #self.read_spPlate(plate)
                #self.read_spZbest(plate)
                #self.get_all_yanny(plate)
                fibers = self.get_okay_cmass()
                fibers = self.redshift_bin_fibers(fibers, zmin, zmax)
                count += len(fibers)
                errors = n.append(errors,self.rm_zerr1[fibers])
                '''
            self.read_redmonster_summary_file()
            for i,z in enumerate(self.rm_z1):
                if (z >= zmin) & (z <= zmax):
                    if (self.rm_type[i] == 'ssp_galaxy_glob') & \
                        (self.rm_zwarning[i] == 0) & (self.rm_zerr1[i] > 0):
                        count += 1
                        errors = n.append(errors,self.rm_zerr1[i])
                        zs = n.append(zs,z)
            #errors.append(self.rm_zerr1[fibers].tolist())
            errors = self.dz_to_dv(zs, errors)
            print(zmin, zmax, n.mean(errors), n.std(errors))
            errors = n.log10(errors)
            hist,binedges = n.histogram(errors, bins=nbins)
            bins = n.zeros(nbins)
            for i in range(nbins):
                bins[i] = (binedges[i+1]+binedges[i])/2.
            normhist = hist / float(count)
            p.plot(bins, normhist, drawstyle='steps-mid', label=labels[j])
        p.minorticks_on()
        p.xlabel(r'$\log_{10} \delta v$ (km s$^{-1}$)', size=14)
        p.ylabel(r'Fraction per bin in $\log_{10} \delta v$', size=14)
        p.title('SEQUELS LRGs')
        p.axis([.5,2.5,0,.25])
        p.legend()
        p.subplots_adjust(wspace = .35)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/dv_vs_z_histos.pdf')
        p.clf()


    def sequels_failure_vs_dchi2_sns(self, drchi2max=.02, npoints=150, sns_pal='muted', rm_line_x=0.005):
    # Makes a plot of SEQUELS LRG failure rate as a function of
    # dchi2 threshold for redmonster and idlspec1d
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')
        rm_data = []
        idl_data = []
        diffs = n.linspace(0,drchi2max,npoints)
        for i,diff in enumerate(diffs):
            print('%s of %s' % (i+1,npoints))
            rm_point, idl_point = self.dchi2_failure_diff_function(diff)
            rm_data.append(rm_point)
            idl_data.append(idl_point)
        f = p.figure()
        ax = f.add_subplot(111)
        p.plot(diffs, rm_data, color=sns.color_palette("RdBu_r", 7)[-1], label='redmonster')
        p.plot(diffs, idl_data, color=sns.color_palette("RdBu_r", 7)[0], label='spectro1d')
        rmcoords01 = (0.01, rm_data[n.abs(n.array(diffs)-0.01).argmin()])
        rmcoords005 = (0.005, rm_data[n.abs(n.array(diffs)-0.005).argmin()])
        idlcoords01 = (0.01, idl_data[n.abs(n.array(diffs)-0.01).argmin()])
        p.plot(n.linspace(0,0.01,1000),[idlcoords01[1]]*1000, color=sns.color_palette("RdBu_r", 7)[0], linestyle='--')
        p.plot([0.01]*1000, n.linspace(0,idlcoords01[1],1000), color=sns.color_palette("RdBu_r", 7)[0], linestyle='--')
        if rm_line_x == 0.01:
            p.plot(n.linspace(0,0.01,1000), [rmcoords01[1]]*1000, color=sns.color_palette("RdBu_r", 7)[-1], linestyle='--')
            p.plot([0.01]*1000, n.linspace(0,rmcoords01[1],1000), color=sns.color_palette("RdBu_r", 7)[-1], linestyle='--')
        else:
            p.plot(n.linspace(0,0.005,1000), [rmcoords005[1]]*1000, color=sns.color_palette("RdBu_r", 7)[-1], linestyle='--')
            p.plot([0.005]*1000, n.linspace(0,rmcoords005[1],1000), color=sns.color_palette("RdBu_r", 7)[-1], linestyle='--')
        p.xlabel(r'$\Delta\chi_{r}^2$ threshold', size=14)
        p.ylabel(r'Cumulative fraction below threshold', size=14)
        #p.grid(b=True, which='major', color='black', linestyle='--')
        p.legend(loc=2, prop={'size':14})
        p.axis([0,.02,0,.7])
        p.tick_params(labelsize=12)
        p.grid(b=True, which='major', color='lightgrey', linestyle='-')
        f.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/drchi2_vs_failure.pdf')
        p.clf()

    def plate_splits_errors_sns(self, nbins=25, fit=True, normed=True, sns_pal='muted'):
        # redshift pdf from splits of extra-deep plates
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')
        '''
        plates = [7834,7839,7848]
        mjds = [56979,56900,56959]
        self.z1 = []
        self.z2 = []
        self.zerr1 = []
        self.zerr2 = []
        for i,plate in enumerate(plates):
            self.plate_splits_function(plate=plate, mjd=mjds[i],
                                       nbins=nbins, fit=fit)
        '''
        c_kms = 299792.458
        directory = '/uufs/astro.utah.edu/common/home/u0814744/compute/scratch/repeatability'
        hdu = fits.open(directory+'/spAll-%s-repeats_lrg.fits' % self.version)

        thing_ids = []
        object_ids1 = []
        object_ids2 = []
        object_ids = {}

        self.z1 = []
        self.z2 = []
        self.zerr1 = []
        self.zerr2 = []

        for thing_id in hdu[1].data.THING_ID:
            if thing_id not in thing_ids:
                thing_ids.append(thing_id)
                w1 = n.where(hdu[1].data.THING_ID == thing_id)[0][0]
                w2 = n.where(hdu[1].data.THING_ID == thing_id)[0][1]
                object_id1 = (hdu[1].data.PLATE[w1], hdu[1].data.MJD[w1], hdu[1].data.FIBERID[w1]-1)
                object_ids1.append(object_id1)
                object_id2 = (hdu[1].data.PLATE[w2], hdu[1].data.MJD[w2], hdu[1].data.FIBERID[w2]-1)
                object_ids2.append(object_id2)
                object_ids[(hdu[1].data.PLATE[w1], hdu[1].data.MJD[w1], hdu[1].data.FIBERID[w1]-1)] = (hdu[1].data.PLATE[w2], hdu[1].data.MJD[w2], hdu[1].data.FIBERID[w2]-1)


        #hdurm = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits'))
        ioerrors = 0
        for i,object_id1 in enumerate(object_ids):
            stderr.write('\r %s of %s ' % (i+1,len(object_ids)))
            try:
                object_id2 = object_ids[object_id1]

                hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_repeats1' % self.version, '%s' % object_id1[0], '%s' % self.version, 'redmonster-%s-%s.fits' % (object_id1[0],object_id1[1])))
                hdu2 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_repeats2' % self.version, '%s' % object_id2[0], '%s' % self.version, 'redmonster-%s-%s.fits' % (object_id2[0],object_id2[1])))
                fiberind1 = n.where(hdu1[1].data.FIBERID == object_id1[2])[0][0]
                fiberind2 = n.where(hdu2[1].data.FIBERID == object_id2[2])[0][0]
                self.z1.append(hdu1[1].data.Z1[fiberind1])
                self.z2.append(hdu2[1].data.Z1[fiberind2])
                self.zerr1.append(hdu1[1].data.Z_ERR1[fiberind1])
                self.zerr2.append(hdu2[1].data.Z_ERR2[fiberind2])

                #dv.append(n.abs(z1-z2)*c_kms/(1+n.min([z1, z2])))
                #drchi2.append(n.min([rchi21, rchi22]))
            except IndexError:
                print("IndexError")
            except IOError:
                ioerrors += 1
                print("IOError! %s %s" % (repr(object_id1), ioerrors))



        self.z1 = n.array(self.z1)
        self.z2 = n.array(self.z2)
        self.zerr1 = n.array(self.zerr1)
        self.zerr2 = n.array(self.zerr2)
        z_diff = self.z2-self.z1
        zerr_rms = n.sqrt( (self.zerr1**2 + self.zerr2**2) )
        scaled_diff = z_diff / zerr_rms
        while True:
            if n.abs(scaled_diff[n.abs(scaled_diff).argmax()]) > 5:
                scaled_diff = n.delete(scaled_diff, n.abs(scaled_diff).argmax())
            else:
                break
        print(n.max(n.abs(scaled_diff)))
        print(scaled_diff.shape)
        hist,binedges = n.histogram(scaled_diff, bins = nbins)
        if normed:
            normhist = hist / float(self.z1.shape[0])
        else:
            normhist = hist
        bins = n.zeros(nbins)
        for i in range(nbins):
            bins[i] = (binedges[i+1]+binedges[i])/2.
        p.plot(bins, normhist, drawstyle='steps-mid', color='black')

        def fit_func(x, a, sigma, mu): # Gaussian function to fit to histogram
            return a * n.exp( -((x-mu)**2)/(2.*sigma**2) )

        if fit:
            popt, pcov = curve_fit(fit_func, bins, normhist)
            xfit = n.linspace(-4,4,1000)
            yfit = fit_func(xfit, popt[0], popt[1], popt[2])
            p.plot(xfit, yfit, color='mediumpurple')
            p.text(.78*(xfit[-1]-xfit[0])+xfit[0], .78*(1.1*max([max(yfit),max(normhist)])), r'$\sigma_{\mathrm{fit}}=$%.2f' % popt[1])
            p.text(.78*(xfit[-1]-xfit[0])+xfit[0], .72*(1.1*max([max(yfit),max(normhist)])), r'$\mu_{\mathrm{fit}}=$%.2f' % popt[2])

        p.xlabel(r'$(z_2-z_1)/ (\delta z_1^2+$ $\delta z_2^2)^{1/2}$', size=14)
        p.ylabel('Fraction per bin', size=14)
        p.tick_params(labelsize=12)
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/reobs_errors.pdf')
        p.clf()


    def logdrchi2_poly_histos_sns(self, nbins=50, sns_pal='muted'):
        # Histograms of log10 delta rchi2 for 1,2,3,4 poly runs
        hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly1' % self.version, 'redmonsterAll-%s.fits' % self.version))
        hdu2 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly2' % self.version, 'redmonsterAll-%s.fits' % self.version))
        hdu3 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly3' % self.version, 'redmonsterAll-%s.fits' % self.version))
        hdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, 'redmonsterAll-%s.fits' % self.version))
        hdulist = [hdu1, hdu2, hdu3, hdu4]
        labels = ['1 poly', '2 poly', '3 poly', '4 poly']
        sns.set_style('white')
        sns_pal = sns.color_palette(["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"])
        sns.set_palette(sns_pal)
        sns.set_context('paper')
        f = p.figure()
        ax = f.add_subplot(111)
        for j,hdu in enumerate(hdulist):
            x = n.delete(hdu[1].data.RCHI2DIFF, n.where(hdu[1].data.RCHI2DIFF == 0)[0])
            x = n.delete(x, n.where(n.isnan(x) == True))
            hist,binedges = n.histogram(n.log10(x), bins=nbins, normed=True)
            bins = n.zeros(nbins)
            for i in range(nbins):
                bins[i] = (binedges[i+1]+binedges[i])/2.
            p.plot(bins, hist, drawstyle='steps-mid', label=labels[j])
        p.plot([n.log10(0.005)]*1000, n.linspace(0,1.2,1000),linestyle='--')
        p.axis([-4,0,0,1.2])
        p.legend(loc=2,prop={'size':14})
        p.xlabel(r'$\log_{10} \Delta \chi^2 / \mathrm{dof}$', size=14)
        p.ylabel('Distribution', size=14)
        p.tick_params(labelsize=12)
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/drchi2_poly_histos.pdf')


    def fiber_poly_differences(self, sns_pal = sns.color_palette("hls", 8)):
        # Find fibers that are successful with 1 poly but not 4 and vice versa, then plot some examples of each
        hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))
        #hdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, 'redmonsterAll-%s.fits' % self.version))
        yes1no4 = []
        no1yes4 = []
        openplate = None
        openmjd = None
        for i,zwarn1 in enumerate(hdu1[1].data.ZWARNING):
            stderr.write('\r %s of %s ' % (i+1, hdu1[1].data.ZWARNING.shape[0]))
            plate = hdu1[1].data.PLATE[i]
            mjd = hdu1[1].data.MJD[i]
            fiberid = hdu1[1].data.FIBERID[i]
            if openplate != plate or openmjd != mjd:
                hdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, '%s' % plate, self.version, 'redmonster-%s-%s.fits' % (plate,mjd)))
                openplate = plate
                openmjd = mjd
            try:
                fiberind = n.where(hdu4[1].data.FIBERID == fiberid)[0][0]
                if not zwarn1 & 4:
                    if hdu4[1].data.ZWARNING[fiberind] & 4 == 4:
                        fiber = (hdu1[1].data.PLATE[i], hdu1[1].data.MJD[i], hdu1[1].data.FIBERID[i])
                        yes1no4.append(fiber)
                        print("1poly success, 4poly failure: plate %s mjd %s fiber %s" % fiber)
                else:
                    if not hdu4[1].data.ZWARNING[fiberind]:
                        fiber = (hdu1[1].data.PLATE[i], hdu1[1].data.MJD[i], hdu1[1].data.FIBERID[i])
                        no1yes4.append(fiber)
                        print("4poly success, 1poly failure: plate %s mjd %s fiber %s" % fiber)
            except IndexError:
                pass
        for i in range(20):
            objid = yes1no4[i]
            print('yes1no4')
            print('plot %s: plate %s mjd %s fiber %s' % (i, objid[0], objid[1], objid[2]))
            hduidl = fits.open( join( environ['BOSS_SPECTRO_REDUX'], self.version, '%s' % objid[0], 'spPlate-%s-%s.fits' % (objid[0],objid[1]) ) )
            hdurm1 = fits.open( join( environ['REDMONSTER_SPECTRO_REDUX'], self.version, '%s' % objid[0], self.version, 'redmonster-%s-%s.fits' % (objid[0],objid[1]) ) )
            hdurm4 = fits.open( join( environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, '%s' % objid[0], self.version,'redmonster-%s-%s.fits'% (objid[0], objid[1]) ) )
            print('1poly z = %s' % hdurm1[1].data.Z1[n.where(hdurm1[1].data.FIBERID == objid[2])[0][0]])
            print('1poly template = %s' % hdurm1[1].data.CLASS1[n.where(hdurm1[1].data.FIBERID == objid[2])[0][0]])
            print('1poly template amplitude = %s' % eval(hdurm1[1].data.THETA1[n.where(hdurm1[1].data.FIBERID == objid[2])[0][0]])[0])
            print('4poly z = %s' % hdurm4[1].data.Z1[n.where(hdurm4[1].data.FIBERID == objid[2])[0][0]])
            print('4poly template = %s' % hdurm4[1].data.CLASS1[n.where(hdurm4[1].data.FIBERID == objid[2])[0][0]])
            print('4poly template amplitude = %s' % eval(hdurm4[1].data.THETA1[n.where(hdurm4[1].data.FIBERID == objid[2])[0][0]])[0])
            print('')
            sns.set_style('white')
            sns.set_palette(sns_pal)
            sns.set_context('paper')
            f = p.figure()
            ax = f.add_subplot(211)
            wave = 10**(hduidl[0].header['COEFF0'] + n.arange(hduidl[0].header['NAXIS1'])*hduidl[0].header['COEFF1'])
            p.plot(wave, convolve(hduidl[0].data[objid[2]], Box1DKernel(5)), color='black', label='Data')
            p.plot(wave, hdurm1[2].data[n.where(hdurm1[1].data.FIBERID == objid[2])[0][0]][0], color=sns_pal[0], label='1 polynomial model')
            p.xlim(wave[0], wave[-1])
            p.ylim(n.sort(hduidl[0].data[objid[2]])[n.round(hduidl[0].data[objid[2]].shape[0]*.05)],
                   n.sort(hduidl[0].data[objid[2]])[n.round(hduidl[0].data[objid[2]].shape[0]*.95)])
            p.legend(loc=2)
            p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
            ax = f.add_subplot(212)
            p.plot(wave, convolve(hduidl[0].data[objid[2]], Box1DKernel(5)), color='black', label='Data')
            p.plot(wave, hdurm4[2].data[n.where(hdurm4[1].data.FIBERID == objid[2])[0][0]][0], color=sns_pal[0], label='4 polynomial model')
            p.xlim(wave[0], wave[-1])
            p.xlim(wave[0], wave[-1])
            p.ylim(n.sort(hduidl[0].data[objid[2]])[n.round(hduidl[0].data[objid[2]].shape[0]*.05)],
                   n.sort(hduidl[0].data[objid[2]])[n.round(hduidl[0].data[objid[2]].shape[0]*.95)])
            p.legend(loc=2)
            p.xlabel(r'Observed wavelength ($\AA$)')
            p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
            p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/yes1no4_%s.pdf' % i)
            p.close()

        for i in range(20):
            objid = no1yes4[i]
            print('no1yes4')
            print('plot %s: plate %s mjd %s fiber %s' % (i, objid[0], objid[1], objid[2]))
            hduidl = fits.open( join( environ['BOSS_SPECTRO_REDUX'], self.version, '%s' % objid[0], 'spPlate-%s-%s.fits' % (objid[0],objid[1]) ) )
            hdurm1 = fits.open( join( environ['REDMONSTER_SPECTRO_REDUX'], self.version, '%s' % objid[0], self.version, 'redmonster-%s-%s.fits' % (objid[0],objid[1]) ) )
            hdurm4 = fits.open( join( environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, '%s' % objid[0], self.version,'redmonster-%s-%s.fits'% (objid[0], objid[1]) ) )
            print('1poly z = %s' % hdurm1[1].data.Z1[n.where(hdurm1[1].data.FIBERID == objid[2])[0][0]])
            print('1poly template = %s' % hdurm1[1].data.CLASS1[n.where(hdurm1[1].data.FIBERID == objid[2])[0][0]])
            print('1poly template amplitude = %s' % eval(hdurm1[1].data.THETA1[n.where(hdurm1[1].data.FIBERID == objid[2])[0][0]])[0])
            print('4poly z = %s' % hdurm4[1].data.Z1[n.where(hdurm4[1].data.FIBERID == objid[2])[0][0]])
            print('4poly template = %s' % hdurm4[1].data.CLASS1[n.where(hdurm4[1].data.FIBERID == objid[2])[0][0]])
            print('4poly template amplitude = %s' % eval(hdurm4[1].data.THETA1[n.where(hdurm4[1].data.FIBERID == objid[2])[0][0]])[0])
            print('')
            sns.set_style('white')
            sns.set_palette(sns_pal)
            sns.set_context('paper')
            f = p.figure()
            ax = f.add_subplot(211)
            wave = 10**(hduidl[0].header['COEFF0'] + n.arange(hduidl[0].header['NAXIS1'])*hduidl[0].header['COEFF1'])
            p.plot(wave, convolve(hduidl[0].data[objid[2]], Box1DKernel(5)), color='black', label='Data')
            p.plot(wave, hdurm1[2].data[n.where(hdurm1[1].data.FIBERID == objid[2])[0][0]][0], color=sns_pal[0], label='1 polynomial model')
            p.xlim(wave[0], wave[-1])
            p.ylim(n.sort(hduidl[0].data[objid[2]])[n.round(hduidl[0].data[objid[2]].shape[0]*.05)],
                   n.sort(hduidl[0].data[objid[2]])[n.round(hduidl[0].data[objid[2]].shape[0]*.95)])
            p.legend(loc=4)
            p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
            ax = f.add_subplot(212)
            p.plot(wave, convolve(hduidl[0].data[objid[2]], Box1DKernel(5)), color='black', label='Data')
            p.plot(wave, hdurm4[2].data[n.where(hdurm4[1].data.FIBERID == objid[2])[0][0]][0], color=sns_pal[0], label='4 polynomial model')
            p.xlim(wave[0], wave[-1])
            p.xlim(wave[0], wave[-1])
            p.ylim(n.sort(hduidl[0].data[objid[2]])[n.round(hduidl[0].data[objid[2]].shape[0]*.05)],
                   n.sort(hduidl[0].data[objid[2]])[n.round(hduidl[0].data[objid[2]].shape[0]*.95)])
            p.legend(loc=4)
            p.xlabel(r'Observed wavelength ($\AA$)')
            p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
            p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/no1yes4_%s.pdf' % i)
            p.close()


    def chi2_compare_poly_sns(self, sns_pal='deep'):
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')
        hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))
        #hdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, 'redmonsterAll-%s.fits' % self.version))
        chi201 = n.array([])
        chi201_yes1no4 = n.array([])
        chi201_no1yes4 = n.array([])
        chi204 = n.array([])
        chi204_yes1no4 = n.array([])
        chi204_no1yes4 = n.array([])
        chi2null1 = n.array([])
        chi2null1_yes1no4 = n.array([])
        chi2null1_no1yes4 = n.array([])
        chi2null4 = n.array([])
        chi2null4_yes1no4 = n.array([])
        chi2null4_no1yes4 = n.array([])
        openplate = None
        openmjd = None
        for i,zwarn in enumerate(hdu1[1].data.ZWARNING):
            stderr.write('\r %s of %s' % (i+1,hdu1[1].data.ZWARNING.shape[0]))
            plate = hdu1[1].data.PLATE[i]
            mjd = hdu1[1].data.MJD[i]
            fiberid = hdu1[1].data.FIBERID[i]
            if openplate != plate or openmjd != mjd:
                hdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, '%s' % plate, self.version, 'redmonster-%s-%s.fits' % (plate,mjd)))
                openplate = plate
                openmjd = mjd
            try:
                fiberind = n.where(hdu4[1].data.FIBERID == fiberid)[0][0]
                if not zwarn & 4:
                    if not hdu4[1].data.ZWARNING[fiberind]:
                        chi201 = n.append(chi201, hdu1[1].data.SN2DATA[i])
                        chi204 = n.append(chi204, hdu4[1].data.SN2DATA[fiberind])
                        chi2null1 = n.append(chi2null1, hdu1[1].data.CHI2NULL[i])
                        chi2null4 = n.append(chi2null4, hdu4[1].data.CHI2NULL[fiberind])
                    else:
                        chi201_yes1no4 = n.append(chi201_yes1no4, hdu1[1].data.SN2DATA[i])
                        chi204_yes1no4 = n.append(chi204_yes1no4, hdu4[1].data.SN2DATA[fiberind])
                        chi2null1_yes1no4 = n.append(chi2null1_yes1no4, hdu1[1].data.CHI2NULL[i])
                        chi2null4_yes1no4 = n.append(chi2null4_yes1no4, hdu4[1].data.CHI2NULL[fiberind])
                else:
                    if not hdu4[1].data.ZWARNING[fiberind]:
                        chi201_no1yes4 = n.append(chi201_no1yes4, hdu1[1].data.SN2DATA[i])
                        chi204_no1yes4 = n.append(chi204_no1yes4, hdu4[1].data.SN2DATA[fiberind])
                        chi2null1_no1yes4 = n.append(chi2null1_no1yes4, hdu1[1].data.CHI2NULL[i])
                        chi2null4_no1yes4 = n.append(chi2null4_no1yes4, hdu4[1].data.CHI2NULL[fiberind])
            except IndexError:
                pass

        f = p.figure()
        ax = f.add_subplot(211)
        p.plot(n.linspace(0,50000,1000),n.linspace(0,50000,1000), color='black', linestyle='--')
        p.scatter( (chi201-chi2null1), (chi204-chi2null4), s=1, color='black', label='Both', alpha=0.6)
        p.scatter( (chi201_yes1no4-chi2null1_yes1no4), (chi204_yes1no4-chi2null4_yes1no4), s=1, color='tomato', label='1 poly')
        p.scatter( (chi201_no1yes4-chi2null1_no1yes4), (chi204_no1yes4-chi2null4_no1yes4), s=1, color='darkturquoise',label='4 poly')
        p.axis([0,50000,0,50000])
        p.legend(loc=4)
        p.xlabel(r'$\chi_{0}^2-\chi_{\mathrm{null},1}^2$')
        p.ylabel(r'$\chi_{0}^2-\chi_{\mathrm{null},4}^2$')

        ax = f.add_subplot(212)
        p.plot(n.linspace(0,1,1000),n.linspace(0,1,1000), color='black', linestyle='--')
        p.scatter( (chi201-chi2null1)/chi201, (chi204-chi2null4)/chi204, s=1, color='black', label='Both', alpha=0.6)
        p.scatter( (chi201_yes1no4-chi2null1_yes1no4)/chi201_yes1no4, (chi204_yes1no4-chi2null4_yes1no4)/chi204_yes1no4, s=1, color='tomato', label='1 poly')
        p.scatter( (chi201_no1yes4-chi2null1_no1yes4)/chi201_no1yes4, (chi204_no1yes4-chi2null4_no1yes4)/chi204_no1yes4, s=1, color='darkturquoise',label='4 poly')
        p.axis([0,1,0,1])
        p.legend(loc=4)
        p.xlabel(r'$\frac{\chi_{0}^2-\chi_{\mathrm{null},1}^2}{\chi_{0}^2}$')
        p.ylabel(r'$\frac{\chi_{0}^2-\chi_{\mathrm{null},4}^2}{\chi_{0}^2}$')

        p.subplots_adjust(hspace = .4)
        p.gcf().subplots_adjust(bottom=.15)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/chi2_compare.pdf')
        p.close()

        # SNS jointplot test
        '''
        f = p.figure()
        ax = f.add_subplot(111)
        sns.jointplot((chi201-chi2null1)/chi201, (chi204-chi2null4)/chi204, kind='reg')
        sns.jointplot((chi201_yes1no4-chi2null1_yes1no4)/chi201_yes1no4, (chi204_yes1no4-chi2null4_yes1no4)/chi204_yes1no4, kind='reg')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/jointplot.pdf')
        '''
        sns.set()
        sns.set_style('white')

        f = p.figure()
        ax = f.add_subplot(111)
        #g = sns.jointplot((chi201_yes1no4-chi2null1_yes1no4)/chi201_yes1no4, (chi204_yes1no4-chi2null4_yes1no4)/chi204_yes1no4, kind="kde", color="k")
        g = sns.JointGrid((chi201_yes1no4-chi2null1_yes1no4)/chi201_yes1no4, (chi204_yes1no4-chi2null4_yes1no4)/chi204_yes1no4, xlim=(0,1), ylim=(0,1))
        g.plot_joint(sns.kdeplot, shade=True, cmap="Greys", n_levels=7)
        g.plot_joint(p.scatter, color='#e74c3c', s=1.5)
        g.plot_marginals(sns.kdeplot, color="black", shade=True)
        g.ax_joint.collections[0].set_alpha(0)
        g.set_axis_labels(r'$\frac{\chi_{0}^2-\chi_{\mathrm{null},1}^2}{\chi_{0}^2}$', r'$\frac{\chi_{0}^2-\chi_{\mathrm{null},4}^2}{\chi_{0}^2}$')
        p.gcf().subplots_adjust(bottom=.15)
        p.gcf().subplots_adjust(left=.15)
        g.fig.suptitle('1 success, 4 failure')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/jointplot1.pdf')
        p.close()

        f = p.figure()
        ax = f.add_subplot(111)
        g = sns.JointGrid((chi201_no1yes4-chi2null1_no1yes4)/chi201_no1yes4, (chi204_no1yes4-chi2null4_no1yes4)/chi204_no1yes4, xlim=(0,1), ylim=(0,1))
        g.plot_joint(sns.kdeplot, shade=True, cmap="Greys", n_levels=10)
        g.plot_joint(p.scatter, color='#e74c3c', s=1.5)
        g.plot_marginals(sns.kdeplot, color="black", shade=True)
        g.ax_joint.collections[0].set_alpha(0)
        g.set_axis_labels(r'$\frac{\chi_{0}^2-\chi_{\mathrm{null},1}^2}{\chi_{0}^2}$', r'$\frac{\chi_{0}^2-\chi_{\mathrm{null},4}^2}{\chi_{0}^2}$')
        p.gcf().subplots_adjust(bottom=.15)
        p.gcf().subplots_adjust(left=.15)
        g.fig.suptitle('1 failure, 4 success')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/jointplot2.pdf')
        p.close()


        # Fit power law to data points and plot on top as well
        def fit_func(x, a, k, b):
            return a*x**k + b
        popt, pcov = curve_fit(fit_func, (chi201-chi2null1)/chi201, (chi204-chi2null4)/chi204)
        print('power law parameters: a=%s, k=%s, b=%s' % (popt[0], popt[1], popt[2]))

        f = p.figure()
        ax = f.add_subplot(111)
        g = sns.JointGrid((chi201-chi2null1)/chi201, (chi204-chi2null4)/chi204, xlim=(0,1), ylim=(0,1))
        g.plot_joint(sns.kdeplot, shade=False, cmap="Purples_d", n_levels=10)
        g.plot_joint(p.scatter, color='black', s=1, alpha=.3)
        p.tick_params(labelsize=12)
        p.grid(b=True, which='major', color='lightgrey', linestyle='-')
        #p.tight_layout()
        g.plot_marginals(sns.kdeplot, color=sns.color_palette('Purples_d')[2], shade=True)
        g.ax_marg_x.xaxis.grid(True)
        g.ax_marg_y.yaxis.grid(True)
        #g.ax_joint.plot(n.linspace(0,1,1000), fit_func(n.linspace(0,1,1000),popt[0], popt[1], popt[2]), linewidth=1, color='k')
        g.ax_joint.plot(n.linspace(0,1,1000), n.linspace(0,1,1000), ':k')
        g.ax_joint.collections[0].set_alpha(0)
        g.set_axis_labels(r'$(\chi_{0}^2-\chi_{\mathrm{null},1}^2)/\chi_{0}^2$', r'$(\chi_{0}^2-\chi_{\mathrm{null},4}^2)/\chi_{0}^2$', size=16)
        p.gcf().subplots_adjust(bottom=.15)
        p.gcf().subplots_adjust(left=.15)
        #g.fig.suptitle('1 failure, 4 success')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/jointplot3.pdf')
        p.close()
        print('median x: %s' % (n.median((chi201-chi2null1)/chi201)))
        print('median y: %s' % (n.median((chi204-chi2null4)/chi204)))

        # compute KDE of the x, y data points
        from scipy.stats.kde import gaussian_kde
        kde_x = gaussian_kde((chi201-chi2null1)/chi201)
        kde_y = gaussian_kde((chi204-chi2null4)/chi204)

        # Fit a gaussian to each of the KDEs
        '''
        def fit_func(x,a,sigma,mu):
            return a * n.exp( -((x-mu)**2)/(2*sigma**2) )
        '''
        import scipy.special as sse
        def fit_func(x, l, s, m):
            return 0.5*l*n.exp(0.5*l*(2*m+l*s*s-2*x))*sse.erfc((m+l*s*s-x)/(n.sqrt(2)*s))

        poptx, pcov = curve_fit(fit_func, n.linspace(0,1,100), kde_x(n.linspace(0,1,100)), p0=(2,.5,.1))
        print(kde_y(n.linspace(0,1,100)))
        popty, pcov = curve_fit(fit_func, n.linspace(0,1,1000), kde_y(n.linspace(0,1,1000)))
        grid = n.zeros((1000,1000))
        kdex = kde_x(n.linspace(0,1,1000))
        kdey = kde_y(n.linspace(0,1,1000))
        for i in range(1000):
            grid[i] = kdex * kdey[i]
        maxcoords = n.unravel_index(grid.argmax(), (1000,1000))
        print(maxcoords)
        print(n.linspace(0,1,1000)[maxcoords[0]], n.linspace(0,1,1000)[maxcoords[1]])

        #plot gaussian fit over histogram of samples from kde, just to check quality of fit
        sns.set_style('white')
        f = p.figure()
        #ax = f.add_subplot(111)
        #p.plot(n.linspace(0,1,1000), kde_x(n.linspace(0,1,1000)), 'r')
        #p.plot(n.linspace(0,1,1000), fit_func(n.linspace(0,1,1000), poptx[0], poptx[1], poptx[2]), 'k')
        #p.hist((chi201-chi2null1)/chi201, normed=1, alpha=.3, bins=25)
        ax = f.add_subplot(111)
        p.plot(n.linspace(0,1,100), kde_y(n.linspace(0,1,100)), 'r')
        p.plot(n.linspace(0,1,1000), fit_func(n.linspace(0,1,1000), popty[0], popty[1], popty[2]), 'k')
        p.hist((chi204-chi2null4)/chi204, normed=1, alpha=.3, bins=20)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/kde_hist.png')
        p.close()

        f = p.figure()
        ax = f.add_subplot(111)
        g = sns.JointGrid((chi201_no1yes4-chi2null1_no1yes4)/chi201_no1yes4, (chi204_no1yes4-chi2null4_no1yes4)/chi204_no1yes4, xlim=(0,1), ylim=(0,1))
        h = sns.JointGrid((chi201_yes1no4-chi2null1_yes1no4)/chi201_yes1no4, (chi204_yes1no4-chi2null4_yes1no4)/chi204_yes1no4, xlim=(0,1), ylim=(0,1))
        g.plot_joint(sns.kdeplot, shade=True, cmap='Blues')
        h.plot_joint(sns.kdeplot, shade=True, cmap='Greens')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/test.pdf')
        p.close()


    def polynomial_and_template_contribution_sns(self, sns_pal='muted'):
        sns.set_style('whitegrid')
        sns.set_palette(sns_pal)
        sns.set_context('paper')

        intmodel1 = [[],[]]
        intpoly1 = [[],[]]
        inttemp1 = [[],[]]
        intmodel4 = [[],[]]
        intpoly4 = [[],[]]
        inttemp4 = [[],[]]
        for path in iglob(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, '*')):
             if len(basename(path)) == 4:
                 plate = basename(path)
                 print(plate)
                 for filepath in iglob(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, '%s' % plate,
                                            self.version, '*')):
                     if len(basename(filepath)) == 26:
                         hduplate = fits.open(join(environ['BOSS_SPECTRO_REDUX'], self.version, '%s' % plate,
                                                   'spPlate-%s-%s.fits' % (plate, basename(filepath)[16:21])))
                         hdu1 = fits.open(filepath)
                         hdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version,
                                               '%s' % plate, self.version, basename(filepath)))
                         wave = 10**(hduplate[0].header['COEFF0'] + n.arange(hduplate[0].header['NAXIS1']) *
                                     hduplate[0].header['COEFF1'])
                         for i in range(hdu1[2].data.shape[0]):
                             if not hdu1[1].data.ZWARNING[i] & 4:
                                 intmodel1[0].append( trapz(hdu1[2].data[i,0], wave) )
                                 temps = read_ndArch(join(environ['REDMONSTER_TEMPLATES_DIR'],
                                                          hdu1[1].data.FNAME1[i]))[0]
                                 this_temp = temps[eval(hdu1[1].data.MINVECTOR1[i])[:-1]][eval(hdu1[1].data.MINVECTOR1[i])[-1]:eval(hdu1[1].data.MINVECTOR1[i])[-1] + hduplate[0].header['NAXIS1']]
                                 inttemp1[0].append( trapz(this_temp * eval(hdu1[1].data.THETA1[i])[0],wave) )
                                 intpoly1[0].append( trapz(poly_array(1, hduplate[0].header['NAXIS1'])[0] *
                                                        eval(hdu1[1].data.THETA1[i])[1], wave) )
                             else:
                                 intmodel1[1].append( trapz(hdu1[2].data[i,0], wave) )
                                 temps = read_ndArch(join(environ['REDMONSTER_TEMPLATES_DIR'],
                                                          hdu1[1].data.FNAME1[i]))[0]
                                 this_temp = temps[eval(hdu1[1].data.MINVECTOR1[i])[:-1]][eval(hdu1[1].data.MINVECTOR1[i])[-1]:eval(hdu1[1].data.MINVECTOR1[i])[-1] + hduplate[0].header['NAXIS1']]
                                 inttemp1[1].append( trapz(this_temp * eval(hdu1[1].data.THETA1[i])[0],wave) )
                                 intpoly1[1].append( trapz(poly_array(1, hduplate[0].header['NAXIS1'])[0] *
                                                        eval(hdu1[1].data.THETA1[i])[1], wave) )
                             if not hdu4[1].data.ZWARNING[i] & 4:
                                 intmodel4[0].append( trapz(hdu4[2].data[i,0],wave) )
                                 temps = read_ndArch(join(environ['REDMONSTER_TEMPLATES_DIR'],
                                                          hdu4[1].data.FNAME1[i]))[0]
                                 this_temp = temps[eval(hdu4[1].data.MINVECTOR1[i])[:-1]][eval(hdu4[1].data.MINVECTOR1[i])[-1]:eval(hdu4[1].data.MINVECTOR1[i])[-1] + hduplate[0].header['NAXIS1']]
                                 inttemp4[0].append( trapz(this_temp * eval(hdu4[1].data.THETA1[i])[0],wave) )
                                 pmat = n.transpose(poly_array(4, hduplate[0].header['NAXIS1']))
                                 intpoly4[0].append( trapz(n.dot(pmat,eval(hdu4[1].data.THETA1[i])[1:]),wave) )
                             else:
                                 intmodel4[1].append( trapz(hdu4[2].data[i,0],wave) )
                                 temps = read_ndArch(join(environ['REDMONSTER_TEMPLATES_DIR'],
                                                          hdu4[1].data.FNAME1[i]))[0]
                                 this_temp = temps[eval(hdu4[1].data.MINVECTOR1[i])[:-1]][eval(hdu4[1].data.MINVECTOR1[i])[-1]:eval(hdu4[1].data.MINVECTOR1[i])[-1] + hduplate[0].header['NAXIS1']]
                                 inttemp4[1].append( trapz(this_temp * eval(hdu4[1].data.THETA1[i])[0],wave) )
                                 pmat = n.transpose(poly_array(4, hduplate[0].header['NAXIS1']))
                                 intpoly4[1].append( trapz(n.dot(pmat,eval(hdu4[1].data.THETA1[i])[1:]),wave) )

        import pdb; pdb.set_trace()
        f = p.figure()
        ax = f.add_subplot(111)
        plt.scatter(n.array(inttemp1[0])/n.array(intmodel1[0]),n.array(intpoly1[0])/n.array(intmodel1[0]),color='black',s=1,alpha=.5,label='ZWARNING = 0')
        plt.scatter(n.array(inttemp1[1])/n.array(intmodel1[1]), n.array(intpoly1[1])/n.array(intmodel1[1]), color='red', s=1, alpha=.5, label='ZWARNING > 0')
        plt.axis([0,2.5,-3,1])
        plt.xlabel(r'$\frac{\int_{\lambda}^{} \theta_{\mathrm{t}} \, \mathrm{d}\lambda}{\int_{\lambda}^{} \theta \, \mathrm{d}\lambda}$')
        plt.ylabel(r'$\frac{\int_{\lambda}^{} \theta_{\mathrm{p}} \, \mathrm{d}\lambda}{\int_{\lambda}^{} \theta \, \mathrm{d}\lambda}$')
        plt.title('Constant polynomial')
        plt.gcf().subplots_adjust(bottom=.2)
        plt.gcf().subplots_adjust(left=.15)
        plt.legend(loc=3)
        plt.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/poly1_contributions.pdf')
        plt.close()

        f = p.figure()
        ax = f.add_subplot(111)
        plt.scatter(n.array(inttemp4[0])/n.array(intmodel4[0]), n.array(intpoly4[0])/n.array(intmodel4[0]), color='black', s=1, alpha=.5, label='ZWARNING = 0')
        plt.scatter(n.array(inttemp4[1])/n.array(intmodel4[1]), n.array(intpoly4[1])/n.array(intmodel4[1]), color='red', s=1, alpha=.5, label='ZWARNING > 0')
        plt.axis([0,2.5,-2,1])
        plt.xlabel(r'$\frac{\int_{\lambda}^{} \theta_{\mathrm{t}} \, \mathrm{d}\lambda}{\int_{\lambda}^{} \theta \, \mathrm{d}\lambda}$')
        plt.ylabel(r'$\frac{\int_{\lambda}^{} \theta_{\mathrm{p}} \, \mathrm{d}\lambda}{\int_{\lambda}^{} \theta \, \mathrm{d}\lambda}$')
        plt.title('Cubic polynomial')
        plt.gcf().subplots_adjust(bottom=.2)
        plt.gcf().subplots_adjust(left=.15)
        plt.legend(loc=3)
        plt.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/poly4_contributions.pdf')
        plt.close()


    def narrow_band_chi2_sns(self, waverange=[3700,4100], sns_pal='muted'):
        sns.set_style('whitegrid')
        sns.set_palette(sns_pal)
        sns.set_context('paper')

        rchi21 = []
        rchi24 = []
        rchi21_yes1no4 = []
        rchi24_yes1no4 = []
        rchi21_no1yes4 = []
        rchi24_no1yes4 = []

        drchi21 = []
        drchi24 = []
        drchi21_yes1no4 = []
        drchi24_yes1no4 = []
        drchi21_no1yes4 = []
        drchi24_no1yes4 = []

        plate = None
        mjd = None
        fiber = None

        hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))
        hdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, 'redmonsterAll-%s.fits' % self.version))
        plotted = False
        nfibers = hdu1[1].data.ZWARNING.shape[0]

        for i,zwarn in enumerate(hdu1[1].data.ZWARNING):
            print('Object %s of %s' % (i+1,nfibers))
            if not (zwarn & 4 and hdu4[1].data.ZWARNING[i] & 4): # only bother with this fiber if at least one has run has !(zwarn & 4)
                if plate != hdu1[1].data.PLATE[i] or mjd != hdu1[1].data.MJD[i]:
                    plate = hdu1[1].data.PLATE[i]
                    mjd = hdu1[1].data.MJD[i]
                    hduidl = fits.open(join(environ['BOSS_SPECTRO_REDUX'], self.version, '%s' % plate, 'spPlate-%s-%s.fits' % (plate,mjd)))
                    wavearr = 10**(hduidl[0].header['COEFF0'] + n.arange(hduidl[0].header['NAXIS1'])*hduidl[0].header['COEFF1'])
                    platehdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, '%s' % plate, self.version, 'redmonster-%s-%s.fits' % (plate,mjd)))
                    platehdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, '%s' % plate, self.version, 'redmonster-%s-%s.fits' % (plate,mjd)))
                fiber = hdu1[1].data.FIBERID[i]
                if not zwarn & 4:
                    this_wave = wavearr / (1 + hdu1[1].data.Z[i])
                else:
                    this_wave = wavearr / (1 + hdu4[1].data.Z[i])
                pix_low = n.abs(this_wave - waverange[0]).argmin()
                pix_high = n.abs(this_wave - waverange[1]).argmin()
                data_slice = hduidl[0].data[fiber][pix_low:pix_high]
                ivar_slice = hduidl[1].data[fiber][pix_low:pix_high]
                model1_slice = platehdu1[2].data[n.where(platehdu1[1].data.FIBERID == fiber)[0][0],0][pix_low:pix_high]
                model4_slice = platehdu4[2].data[n.where(platehdu4[1].data.FIBERID == fiber)[0][0],0][pix_low:pix_high]
                # Repeat for second best model for delta chi2 plot
                this_wave1 = wavearr / (1 + platehdu1[1].data.Z2[n.where(platehdu1[1].data.FIBERID == fiber)[0][0]])
                this_wave4 = wavearr / (1 + platehdu4[1].data.Z2[n.where(platehdu4[1].data.FIBERID == fiber)[0][0]])
                pix_low1 = n.abs(this_wave1 - waverange[0]).argmin()
                pix_high1 = n.abs(this_wave1 - waverange[1]).argmin()
                pix_low4 = n.abs(this_wave4 - waverange[0]).argmin()
                pix_high4 = n.abs(this_wave4 - waverange[1]).argmin()
                data_slice1 = hduidl[0].data[fiber][pix_low1:pix_high1]
                ivar_slice1 = hduidl[1].data[fiber][pix_low1:pix_high1]
                data_slice4 = hduidl[0].data[fiber][pix_low4:pix_high4]
                ivar_slice4 = hduidl[1].data[fiber][pix_low4:pix_high4]
                model1_slice2 = platehdu1[2].data[n.where(platehdu1[1].data.FIBERID == fiber)[0][0],1][pix_low1:pix_high1]
                model4_slice2 = platehdu4[2].data[n.where(platehdu4[1].data.FIBERID == fiber)[0][0],1][pix_low4:pix_high4]
                if not zwarn & 4:
                    if not hdu4[1].data.ZWARNING[i] & 4:
                        rchi21.append(n.sum(((data_slice - model1_slice)**2)*ivar_slice)/data_slice.shape[0])
                        rchi24.append(n.sum(((data_slice - model4_slice)**2)*ivar_slice)/data_slice.shape[0])

                        drchi21.append(n.sum(((data_slice1 - model1_slice2)**2)*ivar_slice1)/data_slice1.shape[0] - rchi21[-1])
                        drchi24.append(n.sum(((data_slice4 - model4_slice2)**2)*ivar_slice4)/data_slice4.shape[0] - rchi24[-1])
                    else:
                        rchi21_yes1no4.append(n.sum(((data_slice - model1_slice)**2)*ivar_slice)/data_slice.shape[0])
                        rchi24_yes1no4.append(n.sum(((data_slice - model4_slice)**2)*ivar_slice)/data_slice.shape[0])

                        drchi21_yes1no4.append(n.sum(((data_slice1 - model1_slice2)**2)*ivar_slice1)/data_slice1.shape[0] - rchi21_yes1no4[-1])
                        drchi24_yes1no4.append(n.sum(((data_slice4 - model4_slice2)**2)*ivar_slice4)/data_slice4.shape[0] - rchi24_yes1no4[-1])
                        if not plotted:
                            if n.random.uniform() < .01:
                                f = p.figure()
                                ax = f.add_subplot(211)
                                p.plot(this_wave[pix_low:pix_high], data_slice, color='black')
                                p.plot(this_wave[pix_low:pix_high], model1_slice, color='cyan')
                                p.title('%s' % (n.sum(((data_slice - model1_slice)**2)*ivar_slice)/data_slice.shape[0]))
                                ax = f.add_subplot(212)
                                p.plot(this_wave[pix_low:pix_high], data_slice, color='black')
                                p.plot(this_wave[pix_low:pix_high], model4_slice, color='cyan')
                                p.title('%s' % (n.sum(((data_slice - model4_slice)**2)*ivar_slice)/data_slice.shape[0]))
                                p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/narrow_test.pdf')
                                p.close()
                                print('Plotted!')
                                time.sleep(2)
                                plotted = True
                else:
                    rchi21_no1yes4.append(n.sum(((data_slice - model1_slice)**2)*ivar_slice)/data_slice.shape[0])
                    rchi24_no1yes4.append(n.sum(((data_slice - model4_slice)**2)*ivar_slice)/data_slice.shape[0])

                    drchi21_no1yes4.append(n.sum(((data_slice1 - model1_slice2)**2)*ivar_slice1)/data_slice1.shape[0] - rchi21_no1yes4[-1])
                    drchi24_no1yes4.append(n.sum(((data_slice4 - model4_slice2)**2)*ivar_slice4)/data_slice4.shape[0] - rchi24_no1yes4[-1])
        f = p.figure()
        ax = f.add_subplot(111)
        p.plot(n.linspace(0.4,1.6,1000), n.linspace(0.4,1.6,1000), '--', color='black')
        p.scatter(rchi21, rchi24, s=1, color='black', label='Both', alpha=0.6)
        p.scatter(rchi21_yes1no4, rchi24_yes1no4, s=1, color='tomato', label='1 poly')
        p.scatter(rchi21_no1yes4, rchi24_no1yes4, s=1, color='darkturquoise', label='4 poly')
        p.axis([0.4,1.6,0.4,1.6])
        p.legend(loc=2)
        p.xlabel(r'$\chi_1^2 / \mathrm{dof}$')
        p.ylabel(r'$\chi_4^2 / \mathrm{dof}$')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/narrow_band_chi2.pdf')
        p.close()

        f = p.figure()
        ax = f.add_subplot(111)
        p.plot(n.linspace(0.4,1.6,1000), n.linspace(0.4,1.6,1000), '--', color='black')
        p.scatter(rchi21_yes1no4, rchi24_yes1no4, s=1, color='tomato', label='1 poly')
        p.scatter(rchi21_no1yes4, rchi24_no1yes4, s=1, color='darkturquoise', label='4 poly')
        p.axis([0.4,1.6,0.4,1.6])
        p.legend(loc=2)
        p.xlabel(r'$\chi_1^2 / \mathrm{dof}$')
        p.ylabel(r'$\chi_4^2 / \mathrm{dof}$')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/narrow_band_chi2_only_failures.pdf')
        p.close()

        f = p.figure()
        ax = f.add_subplot(111)
        p.plot(n.linspace(-1,2,1000), n.linspace(-1,2,1000), '--', color='black')
        p.plot(n.linspace(-1,2,1000), [0]*1000, '--', color='black')
        p.plot([0]*1000, n.linspace(-1,2,1000), '--', color='black')
        p.scatter(drchi21, drchi24, s=1, color='black', label='Both', alpha=0.6)
        p.scatter(drchi21_yes1no4, drchi24_yes1no4, s=1, color='tomato', label='1 poly')
        p.scatter(drchi21_no1yes4, drchi24_no1yes4, s=1, color='darkturquoise', label='4 poly')
        p.axis([-1,2,-1,2])
        p.legend(loc=2)
        p.xlabel(r'$\Delta\chi_1^2 / \mathrm{dof}$')
        p.ylabel(r'$\Delta\chi_4^2 / \mathrm{dof}$')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/narrow_band_dchi2.pdf')

        f = p.figure()
        ax = f.add_subplot(111)
        p.plot(n.linspace(-0.5,1,1000), n.linspace(-0.5,1,1000), '--', color='black')
        p.plot(n.linspace(-0.5,1,1000), [0]*1000,'--', color='black')
        p.plot([0]*1000, n.linspace(-0.5,1,1000),'--', color='black')
        p.scatter(drchi21_yes1no4, drchi24_yes1no4, s=1, color='tomato', label='1 poly')
        p.scatter(drchi21_no1yes4, drchi24_no1yes4, s=1, color='darkturquoise', label='4 poly')
        p.axis([-0.5,1,-0.5,1])
        p.legend(loc=2)
        p.xlabel(r'$\Delta\chi_1^2 / \mathrm{dof}$')
        p.ylabel(r'$\Delta\chi_4^2 / \mathrm{dof}$')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/narrow_band_dchi2_only_failures.pdf')


    def poly_signal_to_noise_histos_sns(self, sns_pal='muted'):
        sns.set_style('whitegrid')
        sns.set_palette(sns_pal)
        sns.set_context('paper')

        yes1no4_r = []
        no1yes4_r = []
        yes1no4_i = []
        no1yes4_i = []
        yes1no4_z = []
        no1yes4_z = []

        hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))
        hdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, 'redmonsterAll-%s.fits' % self.version))

        plate = None
        mjd = None
        nfibers = hdu1[1].data.ZWARNING.shape[0]

        for i,zwarn in enumerate(hdu1[1].data.ZWARNING):
            print('Object %s of %s' % (i+1,nfibers))
            if not (zwarn & 4 and hdu4[1].data.ZWARNING[i] & 4): # only bother with this fiber if at least one has run has !(zwarn & 4)
                if plate != hdu1[1].data.PLATE[i] or mjd != hdu1[1].data.MJD[i]:
                    plate = hdu1[1].data.PLATE[i]
                    mjd = hdu1[1].data.MJD[i]
                    hduidl = fits.open(join(environ['BOSS_SPECTRO_REDUX'], self.version, '%s' % plate, self.version, 'spZbest-%s-%s.fits' % (plate,mjd)))
                    sn_median = hduidl[1].data.SN_MEDIAN[:,2:]

                fiber = hdu1[1].data.FIBERID[i]
                if not zwarn & 4:
                    if hdu4[1].data.ZWARNING[i] & 4:
                        if sn_median[fiber][0] > -1 and sn_median[fiber][0] < 4:
                            yes1no4_r.append(sn_median[fiber][0])
                        if sn_median[fiber][1] > 0 and sn_median[fiber][1] < 6:
                            yes1no4_i.append(sn_median[fiber][1])
                        if sn_median[fiber][2] > 0 and sn_median[fiber][2] < 6:
                            yes1no4_z.append(sn_median[fiber][2])
                else:
                    if not hdu4[1].data.ZWARNING[i] & 4:
                        if sn_median[fiber][0] > -1 and sn_median[fiber][0] < 4:
                            no1yes4_r.append(sn_median[fiber][0])
                        if sn_median[fiber][1] > 0 and sn_median[fiber][1] < 6:
                            no1yes4_i.append(sn_median[fiber][1])
                        if sn_median[fiber][2] > 0 and sn_median[fiber][2] < 6:
                            no1yes4_z.append(sn_median[fiber][2])

        f = p.figure()
        ax = f.add_subplot(311)
        nbins = 25
        hist1, binedges1 = n.histogram(yes1no4_r, bins=nbins, normed=True)
        hist2, binedges2 = n.histogram(no1yes4_r, bins=nbins, normed=True)
        bins1 = n.zeros(nbins)
        bins2 = n.zeros(nbins)
        for i in range(nbins):
            bins1[i] = (binedges1[i+1]+binedges1[i])/2.
            bins2[i] = (binedges2[i+1]+binedges2[i])/2.
        p.plot(bins1, hist1, drawstyle='steps-mid', label='1 poly')
        p.plot(bins2, hist2, drawstyle='steps-mid', label='4 poly')
        lowerx = n.floor( n.min([n.min(bins1), n.min(bins2)]) )
        upperx = n.ceil( n.max([n.max(bins1), n.max(bins2)]) )
        lowery = 0
        uppery = n.around(n.max([n.max(hist1), n.max(hist2)])*1.15,1)
        p.axis([lowerx, upperx, lowery, uppery])
        p.text( (upperx-lowerx)*.03 + lowerx, uppery*.8, '$r$-band', size=8)
        p.legend()
        ax = f.add_subplot(312)
        hist1, binedges1 = n.histogram(yes1no4_i, bins=nbins, normed=True)
        hist2, binedges2 = n.histogram(no1yes4_i, bins=nbins, normed=True)
        bins1 = n.zeros(nbins)
        bins2 = n.zeros(nbins)
        for i in range(nbins):
            bins1[i] = (binedges1[i+1]+binedges1[i])/2.
            bins2[i] = (binedges2[i+1]+binedges2[i])/2.
        p.plot(bins1, hist1, drawstyle='steps-mid', label='1 poly')
        p.plot(bins2, hist2, drawstyle='steps-mid', label='4 poly')
        lowerx = n.floor( n.min([n.min(bins1), n.min(bins2)]) )
        upperx = n.ceil( n.max([n.max(bins1), n.max(bins2)]) )
        lowery = 0
        uppery = n.around(n.max([n.max(hist1), n.max(hist2)])*1.15,1)
        p.axis([lowerx,upperx,lowery,uppery])
        p.text((upperx-lowerx)*.03 + lowerx, uppery*.8,'$i$-band', size=8)
        p.ylabel('Fraction per bin')
        p.legend()
        ax = f.add_subplot(313)
        hist1, binedges1 = n.histogram(yes1no4_z, bins=nbins, normed=True)
        hist2, binedges2 = n.histogram(no1yes4_z, bins=nbins, normed=True)
        bins1 = n.zeros(nbins)
        bins2 = n.zeros(nbins)
        for i in range(nbins):
            bins1[i] = (binedges1[i+1]+binedges1[i])/2.
            bins2[i] = (binedges2[i+1]+binedges2[i])/2.
        p.plot(bins1, hist1, drawstyle='steps-mid', label='1 poly')
        p.plot(bins2, hist2, drawstyle='steps-mid', label='4 poly')
        lowerx = n.floor( n.min([n.min(bins1), n.min(bins2)]) )
        upperx = n.ceil( n.max([n.max(bins1), n.max(bins2)]) )
        lowery = 0
        uppery = n.around(n.max([n.max(hist1), n.max(hist2)])*1.15,1)
        p.axis([lowerx,upperx,lowery,uppery])
        p.text((upperx-lowerx)*.03 + lowerx, uppery*.8,'$z$-band', size=8)
        p.xlabel('Signal to noise ratio')
        p.legend()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/poly_sn_histos.pdf')
        p.close()


    def test_merge_poly_runs(self):
        hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))
        hdu4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly4' % self.version, 'redmonsterAll-%s.fits' % self.version))

        total = 0
        count = 0
        for i,zwarn in enumerate(hdu1[1].data.ZWARNING):
            total += 1.
            if not zwarn & 4:
                count += 1.
            else:
                if not hdu4[1].data.ZWARNING[i] & 4:
                    count += 1
        print(count/total)


    def sequels_sky_drchi2_sns(self, spectro1d=False, nthreshold=50, sns_pal='muted'):
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')

        hdurm = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_sky' % self.version, 'redmonsterAll-%s.fits' % self.version))

        plate = None
        mjd = None

        drchi2_threshold = n.linspace(0,0.01,nthreshold)
        rmfrac = []
        idlfrac = []

        for i,threshold in enumerate(drchi2_threshold):
            stderr.write('\r %s of %s' % (i+1,nthreshold))
            count = 0.
            total = 0.
            countidl = 0.
            totalidl = 0.
            for j,rchi2diff in enumerate(hdurm[1].data.RCHI2DIFF):
                total += 1
                if rchi2diff > threshold:
                    count += 1
                if spectro1d:
                    totalidl += 1
                    if plate != hdurm[1].data.PLATE[j] or mjd != hdurm[1].data.MJD[j]:
                        plate = hdurm[1].data.PLATE[j]
                        mjd = hdurm[1].data.MJD[j]
                        #hduplate = fits.open(join(environ['BOSS_SPECTRO_REDUX'], self.version, '%s' % plate, self.version, 'spZbest-%s-%s.fits' % (plate,mjd)))
                        hduplate = fits.open(join(environ['BOSS_SPECTRO_REDUX'], 'test/bautista/test_dr14', '%s' % plate, 'test_dr14', 'spZbest-%s-%s.fits' % (plate,mjd)))
                    fiber = hdurm[1].data.FIBERID[j]
                    if hduplate[1].data.RCHI2DIFF_NOQSO[fiber] > threshold:
                        countidl += 1
            rmfrac.append(count/total)
            if spectro1d:
                idlfrac.append(countidl/totalidl)

        print(rmfrac)
        f = p.figure()
        ax = f.add_subplot(111)
        if not spectro1d:
            p.plot(drchi2_threshold, rmfrac, color=sns.color_palette(sns_pal)[2], drawstyle='steps-mid')
        else:
            p.plot(drchi2_threshold, rmfrac, drawstyle='steps-mid', color=sns.color_palette(sns_pal)[2], label='redmonster')
            p.plot(drchi2_threshold, idlfrac, drawstyle='steps-mid',color=sns.color_palette(sns_pal)[0], label='spectro1d')
            p.legend(prop={'size':14})
        p.xlabel(r'$\Delta \chi^2/dof$', size=14)
        p.ylabel(r'Cumulative fraction above threshold', size=14)
        ax.set_yscale('log')
        p.tick_params(labelsize=12)
        p.grid(b=True, which='major', color='lightgrey', linestyle='-')
        p.grid(b=True, which='minor', color='lightgrey', linestyle='--')
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/sky_failure_vs_drchi2.pdf')
        p.close()


    def dchi2_dv_repeats(self, sns_pal='muted'):
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')

        c_kms = 299792.458
        directory = '/uufs/astro.utah.edu/common/home/u0814744/compute/scratch/repeatability'
        hdu = fits.open(directory+'/spAll-v5_10_0-repeats_lrg.fits')

        thing_ids = []
        object_ids1 = []
        object_ids2 = []
        object_ids = {}

        dv = []
        drchi2 = []

        for thing_id in hdu[1].data.THING_ID:
            if thing_id not in thing_ids:
                thing_ids.append(thing_id)
                w1 = n.where(hdu[1].data.THING_ID == thing_id)[0][0]
                w2 = n.where(hdu[1].data.THING_ID == thing_id)[0][1]
                object_id1 = (hdu[1].data.PLATE[w1], hdu[1].data.MJD[w1], hdu[1].data.FIBERID[w1]-1)
                object_ids1.append(object_id1)
                object_id2 = (hdu[1].data.PLATE[w2], hdu[1].data.MJD[w2], hdu[1].data.FIBERID[w2]-1)
                object_ids2.append(object_id2)
                object_ids[(hdu[1].data.PLATE[w1], hdu[1].data.MJD[w1], hdu[1].data.FIBERID[w1]-1)] = (hdu[1].data.PLATE[w2], hdu[1].data.MJD[w2], hdu[1].data.FIBERID[w2]-1)


        #hdurm = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits'))
        totalobjs = 0
        for i,object_id1 in enumerate(object_ids):
            stderr.write('\r %s of %s' % (i+1,len(object_ids)))
            try:
                object_id2 = object_ids[object_id1]

                hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_repeats1' % self.version, '%s' % object_id1[0], '%s' % self.version, 'redmonster-%s-%s.fits' % (object_id1[0],object_id1[1])))
                hdu2 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_repeats2' % self.version, '%s' % object_id2[0], '%s' % self.version, 'redmonster-%s-%s.fits' % (object_id2[0],object_id2[1])))
                fiberind1 = n.where(hdu1[1].data.FIBERID == object_id1[2])[0][0]
                fiberind2 = n.where(hdu2[1].data.FIBERID == object_id2[2])[0][0]
                z1 = hdu1[1].data.Z1[fiberind1]
                z2 = hdu2[1].data.Z1[fiberind2]
                rchi21 = hdu1[1].data.RCHI2DIFF[fiberind1]
                rchi22 = hdu2[1].data.RCHI2DIFF[fiberind2]

                dv.append(n.abs(z1-z2)*c_kms/(1+n.min([z1, z2])))
                drchi2.append(n.min([rchi21, rchi22]))
                totalobjs += 1
            except IndexError:
                print("IndexError")
            except IOError:
                ioerrors += 1
                print("IOError! %s %s" % (repr(object_id1), ioerrors))

        dvidl = []
        drchi2idl = []
        for i,object_id1 in enumerate(object_ids):
            stderr.write('\r %s of %s' % (i+1,len(object_ids)))
            try:
                object_id2 = object_ids[object_id1]
                hdu1 = fits.open(join(environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % object_id1[0], '%s' % self.version,
                                      'spZbest-%s-%s.fits' % (object_id1[0],object_id1[1])))
                hdu2 = fits.open(join(environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % object_id2[0], '%s' % self.version,
                                      'spZbest-%s-%s.fits' % (object_id2[0],object_id2[1])))
                z1 = hdu1[1].data.Z_NOQSO[object_id1[2]]
                z2 = hdu2[1].data.Z_NOQSO[object_id2[2]]
                rchi21 = hdu1[1].data.RCHI2DIFF_NOQSO[object_id1[2]]
                rchi22 = hdu2[1].data.RCHI2DIFF_NOQSO[object_id2[2]]

                dvidl.append(n.abs(z1-z2)*c_kms/(1+n.min([z1,z2])))
                drchi2idl.append(n.min([rchi21,rchi22]))
            except IndexError:
                print("IndexError")
            except IOError as e:
                print("IOError")


        print("Total objects: %s" % len(dv)*2)
        confobjs = 0
        cataobjs = 0
        confobjs01 = 0
        cataobjs01 = 0
        confobjs002 = 0
        cataobjs002 = 0
        for i,chi2 in enumerate(drchi2):
            if chi2 > 0.005:
                confobjs += 1.
                if dv[i] > 1000:
                    cataobjs += 1.
            if chi2 > 0.01:
                confobjs01 += 1.
                if dv[i] > 1000:
                    cataobjs01 += 1.
            if chi2 > 0.0015:
                confobjs002 += 1.
                if dv[i] > 1000:
                    cataobjs002 += 1.

        print("Total objects: %s" % (totalobjs))
        print("Redmonster catastrophic failures at 0.002: %s of %s -- %s percent" % (cataobjs002, confobjs002*2, cataobjs002/(confobjs002*2)))
        print("Redmonster catastrophic failures at 0.005: %s of %s -- %s percent" % (cataobjs, confobjs, cataobjs/(confobjs*2)))
        print("Redmonster catastrophic failures at 0.01: %s of %s -- %s percent" % (cataobjs01, confobjs01*2, cataobjs01/(confobjs01*2)))

        print("Total objects: %s" % len(dvidl)*2)
        import pdb; pdb.set_trace()
        confobjs = 0
        cataobjs = 0
        confobjs01 = 0
        cataobjs01 = 0
        for i,chi2 in enumerate(drchi2idl):
            if chi2 > 0.005:
                confobjs += 1.
                if dvidl[i] > 1000:
                    cataobjs += 1.
            if chi2 > 0.01:
                confobjs01 += 1.
                if dvidl[i] > 1000:
                    cataobjs01 += 1
        print("Spectro1d catastrophic failures at 0.005: %s of %s -- %s percent" % (cataobjs, confobjs*2, cataobjs/(confobjs*2)))
        print("Spectro1d catastrophic failures at 0.01: %s of %s -- %s percent" % (cataobjs01, confobjs01*2, cataobjs01/(confobjs01*2)))

        f = p.figure()
        ax = f.add_subplot(111)
        p.scatter(drchi2, dv, alpha=0.4, color='black', s=2)
        p.axis([1e-6, 1, 0.1, 1e6])
        #ax.ylim(0.1, 1e6)
        #ax.xlim(1e-6, 1)
        ax.set_xscale('log')
        ax.set_yscale('log')
        #ylim = p.ylim()
        #p.plot( [1e-2, 1e-2], 1e6, 'b--', lw=2)
        p.plot( [1e-2]*1000, n.linspace(0.1,1e6,1000), '--', color=sns.color_palette('muted')[0], lw=1.5)
        #p.plot( [5e-3, 5e-3], 1e6, 'r--', lw=2)
        p.plot([5e-3]*1000, n.linspace(0.1,1e6,1000), '--', color=sns.color_palette('muted')[2], lw=1.5)
        #p.plot( [1e-6, 1], [1000, 1000], 'm--', lw=2)
        p.plot(n.linspace(1e-6,1,1000), [1000]*1000, '--', color=sns.color_palette('muted')[3], lw=1.5)
        p.xlabel(r'$\Delta \chi^2/dof$', size=14)
        p.ylabel(r'$\Delta v$ (km/s)', size=14)
        p.tick_params(labelsize=12)
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/repeat_dchi2_dv.pdf')


    def make_n_of_z_table(self):
        hdu = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))
        hduidldr13 = fits.open(join(environ['BOSS_SPECTRO_REDUX'], 'v5_8_0', 'spAll-v5_8_0.fits'))
        hduidldr14 = fits.open('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14/spAll-test_dr14.fits')

        spectra = {
                   'total':0.,
                   'poor':0,
                   'stellar':0,
                   '0.0<z<0.5':0,
                   '0.5<z<0.6':0,
                   '0.6<z<0.7':0,
                   '0.7<z<0.8':0,
                   '0.8<z<0.9':0,
                   '0.9<z<1.0':0,
                   '1.0<z<1.1':0,
                   '1.1<z<1.2':0,
                    'z>1.2':0
                    }

        spectraidldr13 = {
                   'total':0.,
                   'poor':0,
                   'stellar':0,
                   '0.0<z<0.5':0,
                   '0.5<z<0.6':0,
                   '0.6<z<0.7':0,
                   '0.7<z<0.8':0,
                   '0.8<z<0.9':0,
                   '0.9<z<1.0':0,
                   '1.0<z<1.1':0,
                   '1.1<z<1.2':0,
                    'z>1.2':0
                    }
        spectraidldr14 = {
                   'total':0.,
                   'poor':0,
                   'stellar':0,
                   '0.0<z<0.5':0,
                   '0.5<z<0.6':0,
                   '0.6<z<0.7':0,
                   '0.7<z<0.8':0,
                   '0.8<z<0.9':0,
                   '0.9<z<1.0':0,
                   '1.0<z<1.1':0,
                   '1.1<z<1.2':0,
                    'z>1.2':0
                    }


        for i,zwarn in enumerate(hdu[1].data.ZWARNING):
            spectra['total'] += 1.
            if zwarn & 4:
                spectra['poor'] += 1.
            else:
                if hdu[1].data.CLASS[i] == 'CAP':
                    spectra['stellar'] += 1.
                else:
                    if hdu[1].data.Z[i] > 0.0 and hdu[1].data.Z[i] < 0.5:
                        spectra['0.0<z<0.5'] += 1.
                    elif hdu[1].data.Z[i] > 0.5 and hdu[1].data.Z[i] < 0.6:
                        spectra['0.5<z<0.6'] += 1.
                    elif hdu[1].data.Z[i] > 0.6 and hdu[1].data.Z[i] < 0.7:
                        spectra['0.6<z<0.7'] += 1.
                    elif hdu[1].data.Z[i] > 0.7 and hdu[1].data.Z[i] < 0.8:
                        spectra['0.7<z<0.8'] += 1.
                    elif hdu[1].data.Z[i] > 0.8 and hdu[1].data.Z[i] < 0.9:
                        spectra['0.8<z<0.9'] += 1.
                    elif hdu[1].data.Z[i] > 0.9 and hdu[1].data.Z[i] < 1.0:
                        spectra['0.9<z<1.0'] += 1.
                    elif hdu[1].data.Z[i] > 1.0 and hdu[1].data.Z[i] < 1.1:
                        spectra['1.0<z<1.1'] += 1.
                    elif hdu[1].data.Z[i] > 1.1 and hdu[1].data.Z[i] < 1.2:
                        spectra['1.1<z<1.2'] += 1.
                    elif hdu[1].data.Z[i] > 1.2:
                        spectra['z>1.2'] += 1.

        for i,ebt1 in enumerate(hduidldr13[1].data.EBOSS_TARGET1):
            if ebt1 & 2:
                if hduidldr13[1].data.SPECPRIMARY[i] > 0:
                    spectraidldr13['total'] += 1.
                    if hduidldr13[1].data.RCHI2DIFF_NOQSO[i] < 0.01:
                        spectraidldr13['poor'] += 1.
                    else:
                        if hduidldr13[1].data.CLASS_NOQSO[i] == 'STAR':
                            spectraidldr13['stellar'] += 1.
                        else:
                            if hduidldr13[1].data.Z_NOQSO[i] > 0.0 and hduidldr13[1].data.Z_NOQSO[i] < 0.5:
                                spectraidldr13['0.0<z<0.5'] += 1.
                            elif hduidldr13[1].data.Z_NOQSO[i] > 0.5 and hduidldr13[1].data.Z_NOQSO[i] < 0.6:
                                spectraidldr13['0.5<z<0.6'] += 1.
                            elif hduidldr13[1].data.Z_NOQSO[i] > 0.6 and hduidldr13[1].data.Z_NOQSO[i] < 0.7:
                                spectraidldr13['0.6<z<0.7'] += 1.
                            elif hduidldr13[1].data.Z_NOQSO[i] > 0.7 and hduidldr13[1].data.Z_NOQSO[i] < 0.8:
                                spectraidldr13['0.7<z<0.8'] += 1.
                            elif hduidldr13[1].data.Z_NOQSO[i] > 0.8 and hduidldr13[1].data.Z_NOQSO[i] < 0.9:
                                spectraidldr13['0.8<z<0.9'] += 1.
                            elif hduidldr13[1].data.Z_NOQSO[i] > 0.9 and hduidldr13[1].data.Z_NOQSO[i] < 1.0:
                                spectraidldr13['0.9<z<1.0'] += 1.
                            elif hduidldr13[1].data.Z_NOQSO[i] > 1.0 and hduidldr13[1].data.Z_NOQSO[i] < 1.1:
                                spectraidldr13['1.0<z<1.1'] += 1.
                            elif hduidldr13[1].data.Z_NOQSO[i] > 1.1 and hduidldr13[1].data.Z_NOQSO[i] < 1.2:
                                spectraidldr13['1.1<z<1.2'] += 1.
                            elif hduidldr13[1].data.Z_NOQSO[i] > 1.2:
                                spectraidldr13['z>1.2'] += 1.

        for i,ebt1 in enumerate(hduidldr14[1].data.EBOSS_TARGET1):
            if ebt1 & 2:
                if hduidldr14[1].data.SPECPRIMARY[i] > 0:
                    spectraidldr14['total'] += 1.
                    if hduidldr14[1].data.RCHI2DIFF_NOQSO[i] < 0.01:
                        spectraidldr14['poor'] += 1.
                    else:
                        if hduidldr14[1].data.CLASS_NOQSO[i] == 'STAR':
                            spectraidldr14['stellar'] += 1.
                        else:
                            if hduidldr14[1].data.Z_NOQSO[i] > 0.0 and hduidldr14[1].data.Z_NOQSO[i] < 0.5:
                                spectraidldr14['0.0<z<0.5'] += 1.
                            elif hduidldr14[1].data.Z_NOQSO[i] > 0.5 and hduidldr14[1].data.Z_NOQSO[i] < 0.6:
                                spectraidldr14['0.5<z<0.6'] += 1.
                            elif hduidldr14[1].data.Z_NOQSO[i] > 0.6 and hduidldr14[1].data.Z_NOQSO[i] < 0.7:
                                spectraidldr14['0.6<z<0.7'] += 1.
                            elif hduidldr14[1].data.Z_NOQSO[i] > 0.7 and hduidldr14[1].data.Z_NOQSO[i] < 0.8:
                                spectraidldr14['0.7<z<0.8'] += 1.
                            elif hduidldr14[1].data.Z_NOQSO[i] > 0.8 and hduidldr14[1].data.Z_NOQSO[i] < 0.9:
                                spectraidldr14['0.8<z<0.9'] += 1.
                            elif hduidldr14[1].data.Z_NOQSO[i] > 0.9 and hduidldr14[1].data.Z_NOQSO[i] < 1.0:
                                spectraidldr14['0.9<z<1.0'] += 1.
                            elif hduidldr14[1].data.Z_NOQSO[i] > 1.0 and hduidldr14[1].data.Z_NOQSO[i] < 1.1:
                                spectraidldr14['1.0<z<1.1'] += 1.
                            elif hduidldr14[1].data.Z_NOQSO[i] > 1.1 and hduidldr14[1].data.Z_NOQSO[i] < 1.2:
                                spectraidldr14['1.1<z<1.2'] += 1.
                            elif hduidldr14[1].data.Z_NOQSO[i] > 1.2:
                                spectraidldr14['z>1.2'] += 1.

        print('REDMONSTER %s' % self.version)
        for entry in spectra:
            print('Fraction %s: %s' % (entry, spectra[entry]/spectra['total']))
            print('N(%s): %s' % (entry, (spectra[entry]/spectra['total'])*60))
        print('Total tracers: %s' % ((spectra['0.6<z<0.7']/spectra['total'])*60 + (spectra['0.7<z<0.8']/spectra['total'])*60 + (spectra['0.8<z<0.9']/spectra['total'])*60 + (spectra['0.9<z<1.0']/spectra['total'])*60))
        print('')
        print('IDL DR13')
        for entry in spectraidldr13:
            print('Fraction %s: %s' % (entry, spectraidldr13[entry]/spectraidldr13['total']))
            print('N(%s): %s' % (entry, (spectraidldr13[entry]/spectraidldr13['total'])*60))
        print('Total tracers: %s' % ((spectraidldr13['0.6<z<0.7']/spectraidldr13['total'])*60 + (spectraidldr13['0.7<z<0.8']/spectraidldr13['total'])*60 + (spectraidldr13['0.8<z<0.9']/spectraidldr13['total'])*60 + (spectraidldr13['0.9<z<1.0']/spectraidldr13['total'])*60))
        print('')
        print('IDL DR14')
        for entry in spectraidldr14:
            try:
                print('Fraction %s: %s' % (entry, spectraidldr14[entry]/spectraidldr14['total']))
                print('N(%s): %s' % (entry, (spectraidldr14[entry]/spectraidldr14['total'])*60))
            except ZeroDivisionError:
                print('%s has no objects at all' % entry)
        print('Total tracers: %s' % ((spectraidldr14['0.6<z<0.7']/spectraidldr14['total'])*60 + (spectraidldr14['0.7<z<0.8']/spectraidldr14['total'])*60 + (spectraidldr14['0.8<z<0.9']/spectraidldr14['total'])*60 + (spectraidldr14['0.9<z<1.0']/spectraidldr14['total'])*60))




    def failure_vs_fiberid(self, sns_pal='muted'):
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')

        hdu = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))

        totals = {}
        counts = {}
        for i in [2*x for x in range(500)]:
            totals[i] = 1
            counts[i] = 0

        for i,zwarn in enumerate(hdu[1].data.ZWARNING):
            fiber = hdu[1].data.FIBERID[i]
            if fiber % 2 == 0:
                totals[fiber] += 1.
                if zwarn > 0:
                    counts[fiber] += 1.

        p.plot(n.array([2*x for x in range(500)])+1, n.array(list(counts.values()))/n.array(list(totals.values())), color=sns.color_palette("Set2", 10)[1], drawstyle='steps-mid')
        p.plot(n.array([2*x for x in range(500)])+1, convolve(n.array(list(counts.values()))/n.array(list(totals.values())), Box1DKernel(5)), color='black', drawstyle='steps-mid')
        #sp.axes([1,1000, 0, n.max( n.array(counts.values())/n.array(totals.values()) )*1.2])
        p.xlabel(r'Fiber number', size=14)
        p.ylabel(r'Failure rate', size=14)
        p.tick_params(labelsize=12)
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/failure_vs_fiberid.pdf')
        p.close()

        total, count = 0., 0.
        for i in [2*x for x in range(250)]:
            total += totals[i]
            count += counts[i]
        print(count/total)
        total, count = 0., 0.
        for i in [2*x+500 for x in range(250)]:
            total += totals[i]
            count += counts[i]
        print(count/total)



    def failure_rate_on_plate(self, nbins=40, sns_pal='muted'):
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')

        hdu = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly1' % self.version, 'redmonsterAll-%s.fits' % self.version))

        p.close()

        xfocal, yfocal = [], []
        xzwarn, yzwarn = [], []
        plate = None
        mjd = None
        for i,zwarn in enumerate(hdu[1].data.ZWARNING):
            stderr.write('\r %s of %s' % (i+1,hdu[1].data.ZWARNING.shape[0]))
            fiberid = hdu[1].data.FIBERID[i]
            if plate != hdu[1].data.PLATE[i] or mjd != hdu[1].data.MJD[i]:
                plate = hdu[1].data.PLATE[i]
                mjd = hdu[1].data.MJD[i]
                hduidl = fits.open(join(environ['BOSS_SPECTRO_REDUX'], self.version, '%s' % plate, 'spPlate-%s-%s.fits' % (plate, mjd)))
                #hduidl = fits.open('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14/%s/spPlate-%s-%s.fits' % (plate,plate,mjd))
            xfocal.append(hduidl[5].data.XFOCAL[fiberid])
            yfocal.append(hduidl[5].data.YFOCAL[fiberid])
            if zwarn > 0:
                xzwarn.append(hduidl[5].data.XFOCAL[fiberid])
                yzwarn.append(hduidl[5].data.YFOCAL[fiberid])
        xfocal = n.array(xfocal)
        yfocal = n.array(yfocal)
        xzwarn = n.array(xzwarn)
        yzwarn = n.array(yzwarn)

        totals, x_edges, y_edges, image = p.hist2d(xfocal, yfocal, bins=nbins, norm=LogNorm())
        totals = n.rot90(totals)
        p.close()
        failures, xbinedges, ybinedges, image = p.hist2d(xzwarn, yzwarn, bins=[x_edges,y_edges], norm=LogNorm())
        failures = n.rot90(failures)
        p.close()

        '''
        xbins = n.zeros(xbinedges.shape[0]-1)
        ybins = n.zeros(ybinedges.shape[0]-1)
        for i in xrange(xbinedges.shape[0]-1):
            xbins[i] = (xbinedges[i+1] + xbinedges[i])/2.
            ybins[i] = (ybinedges[i+1] + ybinedges[i])/2.
        '''

        hist = failures / totals
        p.imshow(hist, interpolation='nearest', origin='lower', extent=[xbinedges[0], xbinedges[-1], ybinedges[0], ybinedges[-1]], cmap='cool')
        cbar = p.colorbar()
        cbar.set_label('Failure rate', size=14)
        cbar.ax.tick_params(labelsize=12)
        p.clim(0,0.25)
        p.tick_params(labelsize=12)
        p.xlabel('XFOCAL', size=14)
        p.ylabel('YFOCAL', size=14)
        f = p.gcf()
        #f.subplots_adjust(bottom=0.2)
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/failure_vs_plate.pdf')
        p.close()

        faildict = {}

        xbins = n.zeros(len(xbinedges)-1)
        ybins = n.zeros(len(ybinedges)-1)
        for i in range(xbins.shape[0]):
            xbins[i] = (xbinedges[i+1] + xbinedges[i])/2.
            ybins[i] = (ybinedges[i+1] + ybinedges[i])/2.

        for i,x in enumerate(xbins):
            for j,y in enumerate(ybins):
                dist = n.floor(n.sqrt(x**2 + y**2))
                if dist <= 300:
                    if dist in faildict:
                        faildict[dist][0] += 1.
                        faildict[dist][1] += hist[i,j]
                    else:
                        faildict[dist] = [1, hist[i,j]]

        fail = []
        dist = []
        for key in faildict:
            dist.append(key/300.)
            fail.append( faildict[key][1]/faildict[key][0])

        f = p.figure()
        f.add_subplot(111)
        p.plot(dist, convolve(fail,Box1DKernel(5)), drawstyle='steps-mid')
        p.xlabel(r'$r/R_\mathrm{plate}$', size=14)
        p.ylabel(r'Failure rate', size=14)
        p.tick_params(labelsize=12)
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/failure_vs_dist.pdf')
        p.close()

        inner_50_fails = []
        for i,this_dist in enumerate(dist):
            if this_dist <= n.sqrt(.5):
                inner_50_fails.append(fail[i])
        inner_50_fails = n.array(inner_50_fails)
        print 'Average failure rate of inner 50 percent of fibers is %s' % (n.sum(inner_50_fails)/inner_50_fails.shape[0])


    def failure_vs_sn_sns(self,sn_max=5,nbins=20):
        # Makes plot of eBOSS LRG target failure rate (zwarning > 0)
        # vs median S/N in r-, i-, and z-bands
        f = p.figure()
        ax = f.add_subplot(1,1,1)
        total = 0
        bad_fibers = []
        bad_r_sn = []
        bad_i_sn = []
        bad_z_sn = []
        r_sn = []
        i_sn = []
        z_sn = []
        rmax = 0
        imax = 0
        zmax = 0
        globpath = join( self.redmonster_spectro_redux,'*')
        openplate = 0
        openmjd = 0
        hdurm = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))
        self.rm_zwarning = hdurm[1].data.ZWARNING
        for i,fiber in enumerate(hdurm[1].data.FIBERID):
            plate = hdurm[1].data.PLATE[i]
            mjd = hdurm[1].data.MJD[i]
            stderr.write('\r %s of %s' % (i+1,hdurm[1].data.FIBERID.shape[0]))
            if (openplate != plate) or (openmjd != mjd):
                hduidl = fits.open(join('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14', '%s' % plate, 'test_dr14', 'spZbest-%s-%s.fits' % (plate, mjd)))
                openplate = plate
                openmjd = mjd
                self.sn_median = hduidl[1].data.SN_MEDIAN[:,2:]
            if (self.sn_median[fiber,0] <= sn_max):
                total += 1
                r_sn.append(self.sn_median[fiber,0])
                if self.sn_median[fiber,0] > rmax:
                    rmax = self.sn_median[fiber,0]
                i_sn.append(self.sn_median[fiber,1])
                if self.sn_median[fiber,1] > imax:
                    imax = self.sn_median[fiber,1]
                z_sn.append(self.sn_median[fiber,2])
                if self.sn_median[fiber,2] > zmax:
                    zmax = self.sn_median[fiber,2]
                if (self.rm_zwarning[i] & 4):
                    bad_fibers.append(fiber)
                    bad_r_sn.append(self.sn_median[fiber,0])
                    bad_i_sn.append(self.sn_median[fiber,1])
                    bad_z_sn.append(self.sn_median[fiber,2])
        nbinsarr = n.linspace(0,sn_max,nbins+1)
        rtotal,rbinedges = n.histogram(r_sn,bins=nbinsarr)
        itotal,ibinedges = n.histogram(i_sn,bins=nbinsarr)
        ztotal,zbinedges = n.histogram(z_sn,bins=nbinsarr)
        rhist,rbinedges = n.histogram(bad_r_sn,bins=nbinsarr)
        ihist,ibinedges = n.histogram(bad_i_sn,bins=nbinsarr)
        zhist,zbinedges = n.histogram(bad_z_sn,bins=nbinsarr)
        rbins = n.zeros(nbins)
        ibins = n.zeros(nbins)
        zbins = n.zeros(nbins)
        for i in range(nbins):
            rbins[i] = (rbinedges[i+1]+rbinedges[i])/2.
            ibins[i] = (ibinedges[i+1]+ibinedges[i])/2.
            zbins[i] = (zbinedges[i+1]+zbinedges[i])/2.
        rhist = rhist / rtotal.astype(float)
        ihist = ihist / itotal.astype(float)
        zhist = zhist / ztotal.astype(float)
        for i in range(nbins):
            if i != 0 and i != (nbins-1):
                if isnan(rhist[i]):
                    try:
                        rhist[i] = (rhist[i-1] + rhist[i+1]) / 2.
                    except:
                        rhist[i] = 0
                if isnan(ihist[i]):
                    try:
                        ihist[i] = (ihist[i-1] + ihist[i+1]) / 2.
                    except:
                        ihist[i] = 0
                if isnan(zhist[i]):
                    try:
                        zhist[i] = (zhist[i-1] + zhist[i+1]) / 2.
                    except:
                        zhist[i] = 0

        p.plot(rbins,rhist,color=sns.color_palette("hls", 8)[4],label='r-band', drawstyle='steps-mid')
        p.plot(ibins,ihist,color=sns.color_palette("hls", 8)[5],label='i-band', drawstyle='steps-mid')
        p.plot(zbins,zhist,color=sns.color_palette("hls", 8)[6],label='z-band', drawstyle='steps-mid')
        ax.set_yscale('log')
        p.xlabel(r'Median S/N per 69 km s$^{-1}$ coadded pixel',size=14)
        p.ylabel(r'eBOSS LRG target failure rate', size=14)
        print(rbins)
        print(rhist)
        print(rtotal)
        print(total)
        print(rmax)
        print(imax)
        print(zmax)
        p.legend(prop={'size':14})
        p.tick_params(labelsize=12)
        p.grid(b=True, which='major', color='lightgrey', linestyle='-')
        p.grid(b=True, which='minor', color='lightgrey', linestyle='--')
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/failure_vs_sn.pdf')
        p.close()


    def failure_vs_imag_sns(self,imin=18,imax=24,nbins=25):
        # Makes plot of SEQUELS LRG failure rate (zwarning > 0)
        # vs i-band magnitude
        f = p.figure()
        ax = f.add_subplot(1,1,1)
        total = 0
        bad_i_mag = []
        i_mag = []
        openplate = 0
        openmjd = 0
        hdurm = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))
        self.rm_zwarning = hdurm[1].data.ZWARNING
        for i,fiber in enumerate(hdurm[1].data.FIBERID):
            plate = hdurm[1].data.PLATE[i]
            mjd = hdurm[1].data.MJD[i]
            stderr.write('\r %s of %s' % (i+1,hdurm[1].data.FIBERID.shape[0]))
            if (openplate != plate) and (openmjd != mjd):
                hduzbest = fits.open(join('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14', '%s' % plate, 'test_dr14', 'spZbest-%s-%s.fits' % (plate, mjd)))
                hduspplate = fits.open(join('/uufs/chpc.utah.edu/common/home/sdss00/ebosswork/eboss/spectro/redux/test/bautista/test_dr14/', '%s' % plate, 'spPlate-%s-%s.fits' % (plate,mjd)))
                self.spectroflux = 22.5 - 2.5*n.log10(hduzbest[1].data.SPECTROFLUX)
                openplate = plate
                openmjd = mjd
            if (self.spectroflux[fiber,3] <= imax):
                total += 1.
                i_mag.append(self.spectroflux[fiber,3])
                if (self.rm_zwarning[i] & 4 > 0):
                    bad_i_mag.append(self.spectroflux[fiber,3])
        nbinsarr = n.linspace(imin,imax,nbins+1)
        itotal,ibinedges = n.histogram(i_mag,bins=nbinsarr)
        ihist,ibinedges = n.histogram(bad_i_mag,bins=nbinsarr)
        ibins = n.zeros(nbins)
        for i in range(nbins):
            ibins[i] = (ibinedges[i+1]+ibinedges[i])/2.
        ihist = ihist / itotal.astype(float)
        for i in range(nbins):
            if i != 0 and i != (nbins-1):
                if isnan(ihist[i]):
                    try:
                        ihist[i] = (ihist[i-1] + ihist[i+1]) / 2.
                    except:
                        ihist[i] = 0
        p.plot(ibins,ihist,drawstyle='steps-mid',label='i-band')
        p.axis([imin,imax,.01,1])
        ax.set_yscale('log')
        p.axvline(21.8,linestyle='--',color='k')
        p.xlabel(r'$i$-band magnitude',size=14)
        p.ylabel(r'Failure rate', size=14)
        #print rbins
        #print rhist
        #print rtotal
        #p.legend()
        p.tick_params(labelsize=12)
        p.grid(b=True, which='major', color='lightgrey', linestyle='-')
        p.grid(b=True, which='minor', color='lightgrey', linestyle='--')
        f.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/failure_vs_imag.pdf')
        p.close


    def zerr_reductions_scatter(self, sns_pal='muted'):
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')

        hdu591 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], 'v5_9_1', 'redmonsterAll-v5_9_1.fits'))

        fibercount = 0
        dz591 = []
        dz5100 = []

        openplate = None
        openmjd = None

        while fibercount < 1000:
            print(fibercount)
            ind591 = n.random.randint(0,hdu591[1].data.ZWARNING.shape[0])
            plate, mjd, fiberid = hdu591[1].data.PLATE[ind591], hdu591[1].data.MJD[ind591], hdu591[1].data.FIBERID[ind591]
            if not hdu591[1].data.ZWARNING[ind591] & 4:
                if exists(join(environ['REDMONSTER_SPECTRO_REDUX'], 'v5_10_0', '%s' % plate, 'v5_10_0', 'redmonster-%s-%s-%s.fits' %
                               (plate, mjd,fiberid))):
                    hdu5100 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], 'v5_10_0', '%s' % plate, 'v5_10_0',
                                             'redmonster-%s-%s-%s.fits' % (plate,mjd,fiberid)))
                    if not hdu5100[1].data.ZWARNING[0] & 4:
                        fibercount += 1
                        dz591.append(hdu591[1].data.Z_ERR[ind591])
                        dz5100.append(hdu5100[1].data.Z_ERR1[0])

        f = p.figure()
        ax = f.add_subplot(111)
        p.plot(n.linspace(0,0.01,100), n.linspace(0,0.01,100), color=sns.color_palette("Set2", 10)[1], linestyle='--', alpha=0.8)
        p.scatter(dz591, dz5100, alpha=0.8, color='black', s=2)
        p.xlabel(r'$\delta z$ (v5_9_1)', size=14)
        p.ylabel(r'$\delta z$ (v5_10_0)', size=14)
        print(max(dz591))
        print(max(dz5100))
        print(min(dz591))
        print(min(dz5100))
        p.axis([min(dz591), max(dz591), min(dz5100), max(dz591)])
        p.tick_params(labelsize=12)
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/zerr_reductions.pdf')


    def completeness_purity_contours(self, sns_pal='muted', chi2min=0, chi2max=0.01, nchi2=50):
        sns.set_style('white')
        sns.set_palette(sns_pal)
        sns.set_context('paper')

        chi2s = n.linspace(chi2min, chi2max, nchi2)

        # Calculate completeness as function of drchi2
        hdurm = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits' % self.version))
        hduidl = fits.open(join(environ['BOSS_SPECTRO_REDUX'], self.version, 'spAll-%s.fits' % self.version))

        compidl = []
        comprm = []
        indices = []

        for i,chi2 in enumerate(chi2s):
            total = 0.
            count = 0.
            for j,rmchi2 in enumerate(hdurm[1].data.RCHI2DIFF):
                stderr.write('\r completeness %s of %s, %s of %s' % (i+1,chi2s.shape[0],j+1,hdurm[1].data.RCHI2DIFF.shape[0]))
                total += 1
                if rmchi2 >= chi2:
                    count += 1
            if total != 0:
                if chi2 != 0.005:
                    comprm.append(count/total)
                else:
                    rmcstar = count/total
            else:
                comprm.append(1.)
            total = 0.
            count = 0.
            if i == 0:
                for j,idlchi2 in enumerate(hduidl[1].data.RCHI2DIFF_NOQSO):
                    stderr.write('\r completeness %s of %s, %s of %s' % (i+1,chi2s.shape[0],j+1,hduidl[1].data.RCHI2DIFF_NOQSO.shape[0]))
                    if hduidl[1].data.EBOSS_TARGET1[j] & 2:
                        indices.append(j)
                        total += 1.
                        if idlchi2 >= chi2:
                            count += 1.
            else:
                for j,ind in enumerate(indices):
                    stderr.write('\r completeness %s of %s, %s of %s' % (i+1,chi2s.shape[0],j+1,len(indices)))
                    total += 1.
                    if hduidl[1].data.RCHI2DIFF_NOQSO[ind] >= chi2:
                        count += 1
            if total != 0:
                if chi2 != 0.01:
                    compidl.append(count/total)
                else:
                    idlcstar = count/total
            else:
                compidl.append(1.)

        # Calculate purity as a function of drchi2
        c_kms = 299792.458
        directory = '/uufs/astro.utah.edu/common/home/u0814744/compute/scratch/repeatability'
        hdu = fits.open(directory+'/spAll-v5_10_0-repeats_lrg.fits')

        thing_ids = []
        object_ids1 = []
        object_ids2 = []
        object_ids = {}

        dv = []
        drchi2 = []

        for thing_id in hdu[1].data.THING_ID:
            if thing_id not in thing_ids:
                thing_ids.append(thing_id)
                w1 = n.where(hdu[1].data.THING_ID == thing_id)[0][0]
                w2 = n.where(hdu[1].data.THING_ID == thing_id)[0][1]
                object_id1 = (hdu[1].data.PLATE[w1], hdu[1].data.MJD[w1], hdu[1].data.FIBERID[w1]-1)
                object_ids1.append(object_id1)
                object_id2 = (hdu[1].data.PLATE[w2], hdu[1].data.MJD[w2], hdu[1].data.FIBERID[w2]-1)
                object_ids2.append(object_id2)
                object_ids[(hdu[1].data.PLATE[w1], hdu[1].data.MJD[w1], hdu[1].data.FIBERID[w1]-1)] = (hdu[1].data.PLATE[w2], hdu[1].data.MJD[w2], hdu[1].data.FIBERID[w2]-1)


        #hdurm = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], self.version, 'redmonsterAll-%s.fits'))
        totalobjs = 0
        for i,object_id1 in enumerate(object_ids):
            stderr.write('\r %s of %s' % (i+1,len(object_ids)))
            try:
                object_id2 = object_ids[object_id1]

                hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_repeats1' % self.version, '%s' % object_id1[0], '%s' % self.version, 'redmonster-%s-%s.fits' % (object_id1[0],object_id1[1])))
                hdu2 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s_repeats2' % self.version, '%s' % object_id2[0], '%s' % self.version, 'redmonster-%s-%s.fits' % (object_id2[0],object_id2[1])))
                fiberind1 = n.where(hdu1[1].data.FIBERID == object_id1[2])[0][0]
                fiberind2 = n.where(hdu2[1].data.FIBERID == object_id2[2])[0][0]
                z1 = hdu1[1].data.Z1[fiberind1]
                z2 = hdu2[1].data.Z1[fiberind2]
                rchi21 = hdu1[1].data.RCHI2DIFF[fiberind1]
                rchi22 = hdu2[1].data.RCHI2DIFF[fiberind2]

                dv.append(n.abs(z1-z2)*c_kms/(1+n.min([z1, z2])))
                drchi2.append(n.min([rchi21, rchi22]))
                totalobjs += 1
            except IndexError:
                print("IndexError")
            except IOError:
                ioerrors += 1
                print("IOError! %s %s" % (repr(object_id1), ioerrors))

        dvidl = []
        drchi2idl = []
        for i,object_id1 in enumerate(object_ids):
            stderr.write('\r %s of %s' % (i+1,len(object_ids)))
            try:
                object_id2 = object_ids[object_id1]
                hdu1 = fits.open(join(environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % object_id1[0], '%s' % self.version, 'spZbest-%s-%s.fits' % (object_id1[0],object_id1[1])))
                hdu2 = fits.open(join(environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % object_id2[0], '%s' % self.version, 'spZbest-%s-%s.fits' % (object_id2[0],object_id2[1])))
                z1 = hdu1[1].data.Z_NOQSO[object_id1[2]]
                z2 = hdu2[1].data.Z_NOQSO[object_id2[2]]
                rchi21 = hdu1[1].data.RCHI2DIFF_NOQSO[object_id1[2]]
                rchi22 = hdu2[1].data.RCHI2DIFF_NOQSO[object_id2[2]]

                dvidl.append(n.abs(z1-z2)*c_kms/(1+n.min([z1,z2])))
                drchi2idl.append(n.min([rchi21,rchi22]))
            except IndexError:
                print("IndexError")
            except IOError as e:
                print("IOError")

        catarm = []
        cataidl = []
        for i,chi2 in enumerate(chi2s):
            stderr.write('\r catastrophic failures %s of %s' % (i+1,chi2s.shape[0]))
            total = 0.
            count = 0.
            for j,rmchi2 in enumerate(drchi2):
                if rmchi2 >= chi2:
                    total += 1.
                    if dv[j] > 1000:
                        count += 1.
            if total != 0:
                if chi2 != 0.005:
                    catarm.append(1-count/(total*2))
                else:
                    rmpstar = 1-count/(total*2)
            else:
                catarm.append(0.)
            total = 0.
            count = 0.
            for j,idlchi2 in enumerate(drchi2idl):
                if idlchi2 >= chi2:
                    total += 1.
                    if dvidl[j] > 1000:
                        count += 1.
            if total != 0:
                if chi2 != 0.01:
                    cataidl.append(1-count/(total*2))
                else:
                    idlpstar = 1-count/(total*2)
            else:
                cataidl.append(0.)

        f = p.figure()
        ax = f.add_subplot(111)
        p.scatter(comprm, catarm, marker='o', c=n.delete(chi2s,12), label='redmonster')
        cbar = p.colorbar()
        cbar.set_label(r'$\Delta\chi^2/$dof', size=14)
        cbar.ax.tick_params(labelsize=12)
        p.scatter(compidl, cataidl, marker='D', c=n.delete(chi2s,24), label='spectro1d')
        p.scatter(rmcstar, rmpstar, marker='*', color='red')
        p.scatter(idlcstar, idlpstar, marker='*', color='red')
        p.xlabel(r'Completeness', size=14)
        p.ylabel(r'Purity', size=14)
        p.legend()
        p.grid(b=True, which='major', color='lightgrey', linestyle='-')
        p.tick_params(labelsize=12)
        p.tight_layout()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/comp_pur_contour.pdf')
