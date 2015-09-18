import numpy as n
from astropy.io import fits
from os.path import join, basename
from os import environ
from glob import iglob
from redmonster.sandbox import yanny as y
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
from math import isnan
from astropy.convolution import convolve, Box1DKernel
from scipy.optimize import curve_fit

class verify_rm:
    
    def __init__(self,version='v5_8_0',plates=[3686,3687,3804,3805,3853,3855,3856,3860],mjds={3686:55268,3687:55269,3804:55267,3805:55269,3853:55268,3855:55268,3856:55269,3860:55269}):
        self.version = version
        self.plates = plates
        self.mjds = mjds
        self.redmonster_spectro_redux = join( environ['REDMONSTER_SPECTRO_REDUX'], '%s' % self.version)
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
        '''
        for i,plate in enumerate(self.plates):
            if i == 0: self.compare_redshifts(plate,self.fibers3686,self.zperson3686)
            elif i == 1: self.compare_redshifts(plate,self.fibers3687,self.zperson3687)
            elif i == 2: self.compare_redshifts(plate,self.fibers3804,self.zperson3804)
            elif i == 3: self.compare_redshifts(plate,self.fibers3805,self.zperson3805)
            elif i == 4: self.compare_redshifts(plate,self.fibers3853,self.zperson3853)
            elif i == 5: self.compare_redshifts(plate,self.fibers3855,self.zperson3855)
            elif i == 6: self.compare_redshifts(plate,self.fibers3856,self.zperson3856)
            elif i == 7: self.compare_redshifts(plate,self.fibers3860,self.zperson3860)
        '''


    def yanny_to_arrays(self):
        # Convert yanny file to arrays
        # Read yanny file
        x = y.yanny(filename='/uufs/astro.utah.edu/common/home/u0814744/boss/spInspect_alltest_bolton.par.txt', np=True)
        #x = y.yanny(filename='/Users/boltonlab3/boss/spInspect_alltest_bolton.par.txt', np=True)
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
        redmonsterpath = join( self.redmonster_spectro_redux, '%s' % plate, '%s' % self.version, 'redmonster-%s-%s.fits' % (plate,self.mjds[plate]) )
        hdu = fits.open(redmonsterpath)
        self.rm_z1 = hdu[1].data.Z1
        self.rm_zerr1 = hdu[1].data.Z_ERR1
        self.rm_fibers = hdu[1].data.FIBERID + 1 # +1 here because rm fibers are 0-based and idlspec2d are 1-based
        self.rm_type = hdu[1].data.CLASS
        self.rm_zwarning = hdu[1].data.ZWARNING
    
    
    def read_spPlate(self,plate):
        # Read in the spPlate file for a given plate
        spPlatepath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, 'spPlate-%s-%s.fits' % (plate, self.mjds[plate]) )
        hdu = fits.open(spPlatepath)
        self.boss_target1 = hdu[5].data.BOSS_TARGET1
    
    
    def read_spZbest(self,plate):
        # Read in the spZbest file for a given plate
        spZbestpath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, '%s' % self.version, 'spZbest-%s-%s.fits' % (plate, self.mjds[plate]) )
        hdu = fits.open(spZbestpath)
        self.sn_median = hdu[1].data.SN_MEDIAN[:,2:]
        self.spectroflux = 22.5 - 2.5*n.log10(hdu[1].data.SPECTROFLUX) # In i-band, note conversion from nanomaggies to magnitudes
    

    def get_cmass(self):
        # Return (0-based) indices of CMASS targets
        return n.where( self.boss_target1 & 2 == 2 )[0].tolist()
    

    def get_lowz(self):
        # Return (0-based indices) of LOWZ targets
        return n.where( self.boss_target1 & 1 == 1 )[0].tolist()
    
    
    def get_okay_cmass(self):
        # Return (0-based) indices of CMASS targets that have the yanny comment 'v5_4_9 ok' and imag <= 21.5
        # self.get_fibers() and self.get_comments() need to have already been called on this plate for this method to work properly
        okay_fibers = (n.asarray(self.vifibers)[n.where(n.asarray(self.comments) == 'v5_4_9 ok')[0].tolist()]-1).tolist() # -1 due to fibers being 1-based and python using 0-based
        return n.asarray(okay_fibers)[n.where( (self.boss_target1[okay_fibers] & 2 == 2) & (self.spectroflux[okay_fibers][:,3] <= 21.5) )[0].tolist()].tolist()
    
    
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
        print count
    

    def cmass_completeness(self):
        # Prints percent of all CMASS targets with rm_zwarning == 0
        vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_cmass()
            vals.append( float(len(n.where( self.rm_zwarning[fibers] == 0 )[0].tolist())) / float(len(fibers)) )
        avg = n.sum(vals) / float(len(vals))
        print avg
                       

    def lowz_completeness(self):
        # Prints percent of all LOWZ targets with rm_zwarning == 0
        vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_lowz()
            vals.append( float(len(n.where( self.rm_zwarning[fibers] == 0 )[0].tolist())) / float(len(fibers)) )
        avg = n.sum(vals) / float(len(vals))
        print avg


    def cmass_galaxy_completeness(self):
        # Prints percent of all CMASS targets that have rm_warning == 0 and were classified as 'ssp_em_galaxy'
        #vals = []
        count = 0
        total = 0
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_cmass()
            for fiber in fibers:
                if self.rm_type[fiber] == 'ssp_em_galaxy':
                    total += 1
                    if self.rm_zwarning[fiber] == 0:
                        count += 1
            #vals.append( float(len(n.where((self.rm_zwarning[fibers] == 0) & (self.rm_type[fibers] == 'ssp_em_galaxy'))[0].tolist())) / float(len(n.where(self.rm_type[fibers] == 'ssp_em_galaxy')[0].tolist())) )
            #vals.append( float(len(n.where((self.rm_zwarning[fibers] == 0) & (self.rm_type[fibers] == 'ssp_em_galaxy'))[0].tolist())) / float(len(fibers)) )
        #avg = n.sum(vals) / float(len(vals))
        avg = float(count) / float(total)
        print count
        print total
        print avg


    def lowz_galaxy_completeness(self):
        # Prints percent of all LOWZ targets that have rm_zwarning == 0 and were classified as 'ssp_em_galaxy'
        #vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_lowz()
            for fiber in fibers:
                if self.rm_type[fiber] == 'ssp_em_galaxy':
                    total += 1
                    if self.rm_zwarning[fiber] == 0:
                        count += 1
            #vals.append( float(len(n.where((self.rm_zwarning[fibers] == 0) & (self.rm_type[fibers] == 'ssp_em_galaxy'))[0].tolist())) / float(len(n.where(self.rm_type[fibers] == 'ssp_em_galaxy')[0].tolist())) )
            #vals.append( float(len(n.where((self.rm_zwarning[fibers] == 0) & (self.rm_type[fibers] == 'ssp_em_galaxy'))[0].tolist())) / float(len(fibers)) )
        #avg = n.sum(vals) / float(len(vals))
        avg = float(count) / float(total)
        print count
        print total
        print avg

    def count_okay_cmass_fibers(self):
        # Prints number of CMASS targets with yanny comment 'v5_4_9 ok' and imag <= 21.5
        count = 0
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            self.get_all_yanny(plate)
            count += len(self.get_okay_cmass())
        print count
            
    def cmass_okay_completeness(self):
        # Prints fraction of CMASS targets having yanny comment 'v5_4_9 ok' and imag <= 21.5 that have rm_zwarning == 0
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
        print '%s out of %s' % (count,total)
        print float(count) / float(total)

    def cmass_okay_galaxy_completeness(self):
        # Prints fraction of targets classified by RM as 'ssp_em_galaxies' in the subset of CMASS targets having yanny comment 'v5_4_9 ok'
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
            count += len( n.where( self.rm_type[fibers] == 'ssp_em_galaxy')[0].tolist() )
        print '%s out of %s' % (count,total)
        print float(count) / float(total)


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
        # Make histograms of log10(dv) in redshift bins for LOWZ and CMASS galaxies
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
            for i in xrange(nbins):
                bins[i] = (binedges[i+1]+binedges[i])/2.
            normhist = hist / float(count)
            p.plot(bins,normhist,drawstyle='steps-mid', color=colors[j], label=labels[j])
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
            for i in xrange(nbins):
                bins[i] = (binedges[i+1]+binedges[i])/2.
            normhist = hist / float(count)
            p.plot(bins,normhist,drawstyle='steps-mid', color=colors[j], label=labels[j])
        p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        p.title('CMASS Sample', size=18)
        p.axis([.9,2.4,0,.3])
        p.legend()
        p.subplots_adjust(wspace = .35)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/dv_histo_cmass.pdf')


    def identify_catastrophic_failures(self):
        # Identify fibers in 'okay CMASS' sample with zwarning == 0 and abs(z_rm - z_person) > .005
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
                if (fiber+1) in self.vifibers: # +1 to go from python indexing to boss fiber conventions
                    total += 1
                    vi_index = n.where( n.asarray(self.vifibers) == (fiber+1) )[0][0]
                    if self.rm_zwarning[fiber] == 0:
                        if self.rm_type[fiber] == 'ssp_em_galaxy':
                            if n.abs(self.rm_z1[fiber] - self.zperson[vi_index]) >= 0.005:
                                self.bad_plates.append(plate)
                                self.bad_fibers.append(fiber)
                                self.bad_rm_z.append(self.rm_z1[fiber])
                                self.bad_zperson.append(self.zperson[vi_index])
                                self.bad_type.append(self.rm_type[fiber])
                                count_bad += 1
        print '%s catastrophic failures out of %s fibers, or %s PERCENT (not fraction!) of the total' % (count_bad,total,(count_bad/float(total))*100)
        for i,fiber in enumerate(self.bad_fibers):
            print 'Plate %s, fiber %s, redmonster z = %s, zperson = %s' % (self.bad_plates[i],fiber,self.bad_rm_z[i], self.bad_zperson[i])


    def identify_unclear_impurities(self):
        # Identify fibers that have zwarning == 0 but no confident visual redshift
        pass


    def identify_recoverable_incompleteness(self):
        # Identify fibers with confident visual redshift and 'galaxy' classification but have zwarning != 0 or rm_type == 'star'
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
                if (fiber+1) in self.vifibers: # +1 to go from python indexing to boss fiber conventions
                    total += 1
                    vi_index = n.where( n.asarray(self.vifibers) == (fiber+1) )[0][0]
                    if (self.rm_zwarning[fiber] != 0) | (self.rm_type[fiber] != 'ssp_em_galaxy'):
                        if (self.zperson[vi_index] != -9) & (self.vitype[vi_index] == 4):
                            self.recoverable_fibers.append(fiber)
                            self.recoverable_plates.append(plate)
                            self.recoverable_rm_z.append(self.rm_z1[fiber])
                            self.recoverable_rm_type.append(self.rm_type[fiber])
                            self.recoverable_zperson.append(self.zperson[vi_index])
                            count_recoverable += 1
        print '%s recoverable failures out of %s fibers, or %s PERCENT (not fraction!) of the total' % (count_recoverable,total,(count_recoverable/float(total))*100)
        for i,fiber in enumerate(self.recoverable_fibers):
            print 'Plate %s, fiber %s, redmonster z = %s, redmonster class = %s, zperson = %s' % (self.recoverable_plates[i],fiber,self.recoverable_rm_z[i], self.recoverable_rm_type[i], self.recoverable_zperson[i])
        big_diff_num = len( n.where( n.abs(n.asarray(self.recoverable_rm_z)-n.asarray(self.recoverable_zperson)) >= .005 )[0] )
        f = p.figure()
        ax1 = f.add_subplot(1,1,1)
        p.plot(self.recoverable_rm_z,self.recoverable_zperson, 'k.')
        p.plot(n.linspace(0,1,1000),n.linspace(0,1,1000),'red')
        p.axis([-.5,3,0,1])
        p.xlabel(r'$z_{redmonster}$',size=16)
        p.ylabel(r'$z_{visual}$',size=16)
        p.title('Objects with "recoverable" redshifts', size=18)
        p.text(1.25, .2, '%s out of %s fibers with confident visual' % (count_recoverable,total), fontsize=10)
        p.text(1.25,.15, 'redshift and called "galaxy" but have', size=10)
        p.text(1.25, .1, 'zwarning > 0 or class != "galaxy". Of', size=10)
        p.text(1.25,.05, 'these, %s have $\delta z > 0.005$.' % (big_diff_num), size=10)
        #p.savefig('recov.pdf')


    def cmass_failure_vs_sn(self,sn_max=7,nbins=29):
        # Makes plot of CMASS failure rate (zwarning > 0) vs median S/N in r-, i-, and z-bands
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
        for i in xrange(nbins):
            rbins[i] = (rbinedges[i+1]+rbinedges[i])/2.
            ibins[i] = (ibinedges[i+1]+ibinedges[i])/2.
            zbins[i] = (zbinedges[i+1]+zbinedges[i])/2.
        rhist = rhist / map(float,rtotal)
        ihist = ihist / map(float,itotal)
        zhist = zhist / map(float,ztotal)
        for i in xrange(nbins):
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


# --------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------- FULL SEQUELS METHODS ONLY BELOW THIS LINE ------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------


    def read_redmonster_all(self,plate, mjd):
        # Read a redmonster file in the context of looking at SEQUELS LRG dataset
        redmonsterpath = join( self.redmonster_spectro_redux, '%s' % plate, '%s' % self.version, 'redmonster-%s-%s.fits' % (plate,mjd) )
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
        summary_path = join( self.redmonster_spectro_redux, 'redmonsterAll-%s.fits' % self.version )
        hdu = fits.open(summary_path)
        self.rm_z1 = hdu[1].data.Z
        self.rm_zerr1 = hdu[1].data.Z_ERR
        #self.rm_fibers = hdu[1].data.FIBERID + 1 # +1 here because rm fibers are 0-based and idlspec2d are 1-based
        self.rm_type = hdu[1].data.CLASS
        self.rm_zwarning = hdu[1].data.ZWARNING
        self.rm_fibers_summary = hdu[1].data.FIBERID
        self.rm_plates_summary = hdu[1].data.PLATE
        self.rm_mjds_summary = hdu[1].data.MJD
        self.rm_rchi2s = hdu[1].data.MINRCHI2
        self.rm_dof = hdu[1].data.DOF
        self.rm_rchi2diff = hdu[1].data.RCHI2DIFF
        self.rm_chi2_null = hdu[1].data.CHI2NULL


    def read_spPlate_all(self,plate, mjd=None):
        # Read in the spPlate file for a given plate in the context of the entire DR10 dataset
        if mjd is not None:
            hdu = fits.open( join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, 'spPlate-%s-%s.fits' % (plate,mjd) ) )
            try: self.eboss_target0 = hdu[5].data.EBOSS_TARGET0
            except: pass
            try: self.eboss_target1 = hdu[5].data.EBOSS_TARGET1
            except: pass
        else:
            globpath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, 'spPlate-%s-*.fits' % plate )
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
            hdu = fits.open(join(environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, '%s' % self.version, 'spZbest-%s-%s.fits' % (plate,mjd) ))
            self.sn_median = hdu[1].data.SN_MEDIAN[:,2:]
            self.spectroflux = 22.5 - 2.5*n.log10(hdu[1].data.SPECTROFLUX) # In i-band, note conversion from nanomaggies to magnitudes
            self.idl_rchi2s = hdu[1].data.RCHI2
            self.idl_dof = hdu[1].data.DOF
            self.idl_rchi2diff = hdu[1].data.RCHI2DIFF_NOQSO
            #self.modelmag = hdu[1].data.MODELMAG[:,2:]
            #self.extinction = hdu[1].data.EXTINCTION[:,2:]
        else:
            globpath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, '%s' % self.version, 'spZbest-%s-*.fits' % plate )
            spZbestpaths = []
            for spZbestpath in iglob(globpath):
                spZbestpaths.append(spZbestpath)
            spZbestpaths.sort()
            hdu = fits.open(spZbestpaths[0])
            self.sn_median = hdu[1].data.SN_MEDIAN[:,2:]
            self.spectroflux = 22.5 - 2.5*n.log10(hdu[1].data.SPECTROFLUX) # In i-band, note conversion from nanomaggies to magnitudes
            self.idl_dof = hdu[1].data.DOF
            self.idl_rchi2diff = hdu[1].data.RCHI2DIFF_NOQSO
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
        print count
        print total
        print avg


    def sequels_galaxy_completeness_all(self):
        # Prints percent of all DR10 CMASS targets that have rm_warning == 0 and were classified as 'ssp_em_galaxy'
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
        #        if (self.rm_zwarning[fiber] == 0) & (self.rm_type[fiber] == 'ssp_em_galaxy'):
        #            count += 1
        self.read_redmonster_summary_file()
        for i,zwarn in enumerate(self.rm_zwarning):
            #if zwarn == 0:
            if (zwarn & 4) == 0:
                total += 1
                if self.rm_type[i] == 'ssp_em_galaxy':
                    count += 1
        avg = float(count) / float(total)
        print count
        print total
        print avg


    def sequels_logdv_vs_z_histos_all(self, nbins=12):
        # Make histograms of log10(dv) in redshift bins for LOWZ and CMASS galaxies
        colors = ['tomato','sage','cornflowerblue','sandybrown','mediumpurple','grey']
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
            p.plot(bins,normhist,drawstyle='steps-mid', color=colors[j], label=labels[j])
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
                    if (self.rm_type[i] == 'ssp_em_galaxy') & (self.rm_zwarning[i] == 0) & (self.rm_zerr1[i] > 0):
                        count += 1
                        errors = n.append(errors,self.rm_zerr1[i])
                        zs = n.append(zs,z)
            #errors.append(self.rm_zerr1[fibers].tolist())
            errors = self.dz_to_dv(zs, errors)
            errors = n.log10(errors)
            hist,binedges = n.histogram(errors, bins=nbins)
            bins = n.zeros(nbins)
            for i in xrange(nbins):
                bins[i] = (binedges[i+1]+binedges[i])/2.
            normhist = hist / float(count)
            p.plot(bins,normhist,drawstyle='steps-mid', color=colors[j], label=labels[j])
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
        # Makes plot of SEQUELS LRG target failure rate (zwarning > 0) vs median S/N in r-, i-, and z-bands
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
            print '%s-%s-%s' % (plate,fiber,mjd)
            if (openplate != plate) and (openmjd != mjd):
                self.read_spZbest_all(plate,mjd)
                openplate = plate
                openmjd = mjd
            if (self.sn_median[fiber,0] <= sn_max):
                total += 1
                r_sn.append(self.sn_median[fiber,0])
                if self.sn_median[fiber,0] > rmax: rmax = self.sn_median[fiber,0]
                i_sn.append(self.sn_median[fiber,1])
                if self.sn_median[fiber,1] > imax: imax = self.sn_median[fiber,1]
                z_sn.append(self.sn_median[fiber,2])
                if self.sn_median[fiber,2] > zmax: zmax = self.sn_median[fiber,2]
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
        for i in xrange(nbins):
            rbins[i] = (rbinedges[i+1]+rbinedges[i])/2.
            ibins[i] = (ibinedges[i+1]+ibinedges[i])/2.
            zbins[i] = (zbinedges[i+1]+zbinedges[i])/2.
        rhist = rhist / map(float,rtotal)
        ihist = ihist / map(float,itotal)
        zhist = zhist / map(float,ztotal)
        for i in xrange(nbins):
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
        print rbins
        print rhist
        print rtotal
        print total
        print rmax
        print imax
        print zmax
        p.legend()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/failure_vs_sn.pdf')
        p.clf()


    def sequels_logdv_vs_sn_histos_all(self, nbins=25):
        # Make histograms of log10(dv) in S/N bins in bands r,i,z for SEQUELS LRG targets
        colors = ['tomato','sage','cornflowerblue','sandybrown','mediumpurple','grey'] #['purple', 'cyan', 'blue', 'lime', 'red', 'black']
        labels = ['1.0<S/N<1.5','1.5<S/N<2.0','2.0<S/N<2.5','2.5<S/N<3.0','3.0<S/N<3.5','3.5<S/N<4.0','4.0<S/N<4.5']
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
                print '%s-%s-%s' % (plate,fiber,mjd)
                if (openplate != plate) and (openmjd != mjd):
                    self.read_spZbest_all(plate,mjd)
                    self.read_spPlate_all(plate,mjd)
                    openplate = plate
                    openmjd = mjd
                if (self.rm_zwarning[i] == 0) & (self.rm_zerr1[i] > 0):
                    if (self.sn_median[fiber][0] >= sn_min) & (self.sn_median[fiber][0] <= sn_max):
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
        errors1 = n.log10(errors1)
        hist1,binedges1 = n.histogram(errors1, bins=nbins)
        bins1 = n.zeros(nbins)
        for i in xrange(nbins):
            bins1[i] = (binedges1[i+1]+binedges1[i])/2.
        normhist1 = hist1 / float(count1)
        p.plot(bins1,normhist1,drawstyle='steps-mid', color=colors[0], label=labels[0])
        errors2 = self.dz_to_dv(z2,errors2)
        errors2 = n.log10(errors2)
        hist2,binedges2 = n.histogram(errors2, bins=nbins)
        bins2 = n.zeros(nbins)
        for i in xrange(nbins):
            bins2[i] = (binedges2[i+1]+binedges2[i])/2.
        normhist2 = hist2 / float(count2)
        p.plot(bins2,normhist2,drawstyle='steps-mid', color=colors[1], label=labels[1])
        errors3 = self.dz_to_dv(z3,errors3)
        errors3 = n.log10(errors3)
        hist3,binedges3 = n.histogram(errors3, bins=nbins)
        bins3 = n.zeros(nbins)
        for i in xrange(nbins):
            bins3[i] = (binedges3[i+1]+binedges3[i])/2.
        normhist3 = hist3 / float(count3)
        p.plot(bins3,normhist3,drawstyle='steps-mid', color=colors[2], label=labels[2])
        errors4 = self.dz_to_dv(z4,errors4)
        errors4 = n.log10(errors4)
        hist4,binedges4 = n.histogram(errors4, bins=nbins)
        bins4 = n.zeros(nbins)
        for i in xrange(nbins):
            bins4[i] = (binedges4[i+1]+binedges4[i])/2.
        normhist4 = hist4 / float(count4)
        p.plot(bins4,normhist4,drawstyle='steps-mid', color=colors[3], label=labels[3])
        errors5 = self.dz_to_dv(z5,errors5)
        errors5 = n.log10(errors5)
        hist5,binedges5 = n.histogram(errors5, bins=nbins)
        bins5 = n.zeros(nbins)
        for i in xrange(nbins):
            bins5[i] = (binedges5[i+1]+binedges5[i])/2.
        normhist5 = hist5 / float(count5)
        #p.plot(bins5,normhist5,drawstyle='steps-mid', color=colors[4], label=labels[4])
        errors6 = self.dz_to_dv(z6,errors6)
        errors6 = n.log10(errors6)
        hist6,binedges6 = n.histogram(errors6, bins=nbins)
        bins6 = n.zeros(nbins)
        for i in xrange(nbins):
            bins6[i] = (binedges6[i+1]+binedges6[i])/2.
        normhist6 = hist6 / float(count6)
        #p.plot(bins6,normhist6,drawstyle='steps-mid', color=colors[5], label=labels[5])
        p.text(0.8, 0.2, '$r$-band', fontsize=12)
        #p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        #p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        #p.title('r-band', size=18)
        p.axis([.5,2.5,0,.25])
        p.legend(prop={'size':6})
        print count1
        print count2
        print count3
        print count4
        print count5
        print count6
        print (count1+count2+count3+count4+count5+count6)/float(self.rm_fibers_summary.shape[0])

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
                print '%s-%s-%s' % (plate,fiber,mjd)
                if (openplate != plate) and (openmjd != mjd):
                    self.read_spZbest_all(plate,mjd)
                    self.read_spPlate_all(plate,mjd)
                    openplate = plate
                    openmjd = mjd
                if (self.rm_zwarning[i] == 0) & (self.rm_zerr1[i] > 0):
                    if (self.sn_median[fiber][1] >= sn_min) & (self.sn_median[fiber][1] <= sn_max):
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
        errors1 = n.log10(errors1)
        hist1,binedges1 = n.histogram(errors1, bins=nbins)
        bins1 = n.zeros(nbins)
        for i in xrange(nbins):
            bins1[i] = (binedges1[i+1]+binedges1[i])/2.
        normhist1 = hist1 / float(count1)
        p.plot(bins1,normhist1,drawstyle='steps-mid', color=colors[0], label=labels[0])
        errors2 = self.dz_to_dv(z2,errors2)
        errors2 = n.log10(errors2)
        hist2,binedges2 = n.histogram(errors2, bins=nbins)
        bins2 = n.zeros(nbins)
        for i in xrange(nbins):
            bins2[i] = (binedges2[i+1]+binedges2[i])/2.
        normhist2 = hist2 / float(count2)
        p.plot(bins2,normhist2,drawstyle='steps-mid', color=colors[1], label=labels[1])
        errors3 = self.dz_to_dv(z3,errors3)
        errors3 = n.log10(errors3)
        hist3,binedges3 = n.histogram(errors3, bins=nbins)
        bins3 = n.zeros(nbins)
        for i in xrange(nbins):
            bins3[i] = (binedges3[i+1]+binedges3[i])/2.
        normhist3 = hist3 / float(count3)
        p.plot(bins3,normhist3,drawstyle='steps-mid', color=colors[2], label=labels[2])
        errors4 = self.dz_to_dv(z4,errors4)
        errors4 = n.log10(errors4)
        hist4,binedges4 = n.histogram(errors4, bins=nbins)
        bins4 = n.zeros(nbins)
        for i in xrange(nbins):
            bins4[i] = (binedges4[i+1]+binedges4[i])/2.
        normhist4 = hist4 / float(count4)
        p.plot(bins4,normhist4,drawstyle='steps-mid', color=colors[3], label=labels[3])
        errors5 = self.dz_to_dv(z5,errors5)
        errors5 = n.log10(errors5)
        hist5,binedges5 = n.histogram(errors5, bins=nbins)
        bins5 = n.zeros(nbins)
        for i in xrange(nbins):
            bins5[i] = (binedges5[i+1]+binedges5[i])/2.
        normhist5 = hist5 / float(count5)
        p.plot(bins5,normhist5,drawstyle='steps-mid', color=colors[4], label=labels[4])
        errors6 = self.dz_to_dv(z6,errors6)
        errors6 = n.log10(errors6)
        hist6,binedges6 = n.histogram(errors6, bins=nbins)
        bins6 = n.zeros(nbins)
        for i in xrange(nbins):
            bins6[i] = (binedges6[i+1]+binedges6[i])/2.
        normhist6 = hist6 / float(count6)
        p.plot(bins6,normhist6,drawstyle='steps-mid', color=colors[5], label=labels[5])
        p.text(0.8, 0.2, '$i$-band', fontsize=12)
        #p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        #p.title('r-band', size=18)
        p.axis([.5,2.5,0,.25])
        p.legend(prop={'size':6})
        print count1
        print count2
        print count3
        print count4
        print count5
        print count6
        print (count1+count2+count3+count4+count5+count6)/float(self.rm_fibers_summary.shape[0])

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
                print '%s-%s-%s' % (plate,fiber,mjd)
                if (openplate != plate) and (openmjd != mjd):
                    self.read_spZbest_all(plate,mjd)
                    self.read_spPlate_all(plate,mjd)
                    openplate = plate
                    openmjd = mjd
                if (self.rm_zwarning[i] == 0) & (self.rm_zerr1[i] > 0):
                    if (self.sn_median[fiber][2] >= sn_min) & (self.sn_median[fiber][2] <= sn_max):
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
        errors1 = n.log10(errors1)
        hist1,binedges1 = n.histogram(errors1, bins=nbins)
        bins1 = n.zeros(nbins)
        for i in xrange(nbins):
            bins1[i] = (binedges1[i+1]+binedges1[i])/2.
        normhist1 = hist1 / float(count1)
        p.plot(bins1,normhist1,drawstyle='steps-mid', color=colors[0], label=labels[1])
        errors2 = self.dz_to_dv(z2,errors2)
        errors2 = n.log10(errors2)
        hist2,binedges2 = n.histogram(errors2, bins=nbins)
        bins2 = n.zeros(nbins)
        for i in xrange(nbins):
            bins2[i] = (binedges2[i+1]+binedges2[i])/2.
        normhist2 = hist2 / float(count2)
        p.plot(bins2,normhist2,drawstyle='steps-mid', color=colors[1], label=labels[2])
        errors3 = self.dz_to_dv(z3,errors3)
        errors3 = n.log10(errors3)
        hist3,binedges3 = n.histogram(errors3, bins=nbins)
        bins3 = n.zeros(nbins)
        for i in xrange(nbins):
            bins3[i] = (binedges3[i+1]+binedges3[i])/2.
        normhist3 = hist3 / float(count3)
        p.plot(bins3,normhist3,drawstyle='steps-mid', color=colors[2], label=labels[3])
        errors4 = self.dz_to_dv(z4,errors4)
        errors4 = n.log10(errors4)
        hist4,binedges4 = n.histogram(errors4, bins=nbins)
        bins4 = n.zeros(nbins)
        for i in xrange(nbins):
            bins4[i] = (binedges4[i+1]+binedges4[i])/2.
        normhist4 = hist4 / float(count4)
        p.plot(bins4,normhist4,drawstyle='steps-mid', color=colors[3], label=labels[4])
        errors5 = self.dz_to_dv(z5,errors5)
        errors5 = n.log10(errors5)
        hist5,binedges5 = n.histogram(errors5, bins=nbins)
        bins5 = n.zeros(nbins)
        for i in xrange(nbins):
            bins5[i] = (binedges5[i+1]+binedges5[i])/2.
        normhist5 = hist5 / float(count5)
        p.plot(bins5,normhist5,drawstyle='steps-mid', color=colors[4], label=labels[5])
        errors6 = self.dz_to_dv(z6,errors6)
        errors6 = n.log10(errors6)
        hist6,binedges6 = n.histogram(errors6, bins=nbins)
        bins6 = n.zeros(nbins)
        for i in xrange(nbins):
            bins6[i] = (binedges6[i+1]+binedges6[i])/2.
        normhist6 = hist6 / float(count6)
        p.plot(bins6,normhist6,drawstyle='steps-mid', color=colors[5], label=labels[6])
        p.text(0.8, 0.28, '$z$-band', fontsize=12)
        p.xlabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        #p.ylabel(r'Fraction per bin in $\log_{10} \delta$v', size=16)
        p.axis([.5,2.5,0,.35])
        p.legend(prop={'size':6})
        p.subplots_adjust(hspace = .5)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/dv_vs_sn_histos.pdf')
        p.clf()
        print count1
        print count2
        print count3
        print count4
        print count5
        print count6
        print (count1+count2+count3+count4+count5+count6)/float(self.rm_fibers_summary.shape[0])


    def sequels_failure_vs_imag_all(self,imin=18,imax=24,nbins=21):
        # Makes plot of SEQUELS LRG failure rate (zwarning > 0) vs i-band magnitude
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
            print '%s-%s-%s' % (plate,mjd,fiber)
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
        for i in xrange(nbins):
            ibins[i] = (ibinedges[i+1]+ibinedges[i])/2.
        ihist = ihist / map(float,itotal)
        for i in xrange(nbins):
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


    def cmass_logdv_vs_z_scatter_all(self,nobjs=100000):
    # Makes a scatterplot nobjs CMASS targets of redshift vs log(dv)
        self.read_redmonster_summary_file()
        errors = []
        zs = []
        for i in xrange(nobjs):
            if (self.rm_zwarning[i] == 0) & (self.rm_type[i] == 'ssp_em_galaxy') & (self.rm_zerr1[i] != -1):
                errors.append(self.rm_zerr1[i])
                zs.append(self.rm_z1[i])
        errors = self.dz_to_dv(n.asarray(errors))
        logerrs = n.log10(errors)
        p.scatter(zs,logerrs, marker='.')
        p.axhline(2.48,linestyle='--',color='k')
        p.xlabel('Redshift',size=16)
        p.ylabel(r'$\log_{10} \delta$v (km s$^{-1}$)', size=16)
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/dv_vs_z_scatter.pdf')
        p.clf()


    def sequels_chi2_histos(self,nbins=50, rchi2=True):
        # Makes histogram of SEQUELS chi2 values for redmonster and idlspec1d (in chi2 or rchi2)
        rm_rchi2s = []
        idl_rchi2s = []
        openplate = 0
        openmjd = 0
        total = 0
        self.read_redmonster_summary_file()
        for i,fiber in enumerate(self.rm_fibers_summary):
            plate = self.rm_plates_summary[i]
            mjd = self.rm_mjds_summary[i]
            print '%s-%s-%s' % (plate,mjd,fiber)
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
                    idl_rchi2s.append(self.idl_rchi2s[fiber] * self.idl_dof[fiber])
                #rm_rchi2s.append(self.rm_rchi2diff[i])
                #idl_rchi2s.append(self.idl_rchi2diff[fiber])
        rmhist,rmbinedges = n.histogram(rm_rchi2s,nbins)
        rmbins = n.zeros(nbins)
        for i in xrange(nbins):
            rmbins[i] = (rmbinedges[i+1]+rmbinedges[i])/2.
        rmhist = rmhist / float(total)
        idlhist, idlbinedges = n.histogram(idl_rchi2s,nbins)
        idlbins = n.zeros(nbins)
        for i in xrange(nbins):
            idlbins[i] = (idlbinedges[i+1]+idlbinedges[i])/2.
        idlhist = idlhist / float(total)
        p.plot(rmbins, rmhist, color='red', drawstyle='steps-mid', label='redmonster')
        p.plot(idlbins, idlhist, color='blue', drawstyle='steps-mid', label='idlspec1d')
        p.xlabel(r'$\chi_r^2$', size=16)
        p.ylabel(r'Fraction per bin', size=16)
        p.legend()
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/rchi2_histos.pdf')
        p.clf()


    def sequels_drchi2_histos(self, drchi2max=.02, nbins=50, rchi2=True):
        # Makes histogram of SEQUELS delta-rchi2 (or delta-chi2) values for redmonster and idlspec1d
        rm_drchi2s = []
        idl_drchi2s = []
        openplate = 0
        openmjd = 0
        total = 0
        self.read_redmonster_summary_file()
        for i,fiber in enumerate(self.rm_fibers_summary):
            plate = self.rm_plates_summary[i]
            mjd = self.rm_mjds_summary[i]
            print '%s-%s-%s' % (plate,mjd,fiber)
            if (openplate != plate) and (openmjd != mjd):
                self.read_spZbest_all(plate,mjd)
                self.read_spPlate_all(plate,mjd)
                openplate = plate
                openmjd = mjd
            if (self.rm_rchi2diff[i] < drchi2max) and (self.idl_rchi2diff[fiber] < drchi2max):
                total += 1
                if rchi2:
                    rm_drchi2s.append( self.rm_rchi2diff[i] )
                    idl_drchi2s.append( self.idl_rchi2diff[fiber] )
                else:
                    rm_drchi2s.append( self.rm_rchi2diff[i] * self.rm_dof[i] )
                    idl_drchi2s.append( self.idl_rchi2diff[fiber] * self.idl_dof[fiber] )
        rmhist, rmbinedges = n.histogram(rm_drchi2s,nbins)
        rmbins = n.zeros(nbins)
        for i in xrange(nbins):
            rmbins[i] = (rmbinedges[i+1]+rmbinedges[i])/2.
        rmhist = rmhist / float(total)
        idlhist, idlbinedges = n.histogram(idl_drchi2s,nbins)
        idlbins = n.zeros(nbins)
        for i in xrange(nbins):
            idlbins[i] = (idlbinedges[i+1]+idlbinedges[i])/2.
        idlhist = idlhist / float(total)
        p.plot(rmbins, rmhist, color='red', drawstyle='steps-mid', label='redmonster')
        p.plot(idlbins, idlhist, color='blue', drawstyle='steps-mid', label='idlspec1d')
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
                self.read_spZbest_all(plate,mjd)
                self.read_spPlate_all(plate,mjd)
                openplate = plate
                openmjd = mjd
            #if (self.rm_rchi2diff[i] < drchi2max) and (self.idl_rchi2diff[fiber] < drchi2max):
            total += 1
            if self.rm_rchi2diff[i] < diff: rm_failures += 1.
            if self.idl_rchi2diff[fiber] < diff: idl_failures += 1.
        return (rm_failures/total), (idl_failures/total)


    def sequels_failure_vs_dchi2(self, drchi2max=.02, npoints=150):
        # Makes a plot of SEQUELS LRG failure rate as a function of dchi2 threshold for redmonster and idlspec1d
        rm_data = []
        idl_data = []
        diffs = n.linspace(0,drchi2max,npoints)
        for i,diff in enumerate(diffs):
            print '%s of %s' % (i+1,npoints)
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
        # Makes a histogram of (z2-z1)/sqrt(dz1**2 + dz2**2) with best fit Gaussian overplotted for all SEQUELS LRG targets with repeat observations
        globpath = join( self.redmonster_spectro_redux,'*')
        z1 = []
        z2 = []
        zerr1 = []
        zerr2 = []
        for path in iglob(globpath):
            plate = basename(path)
            if plate != 'redmonsterAll-%s.fits' % self.version:
                print plate
                mjds = []
                mjdglobpath = join( self.redmonster_spectro_redux, '%s' % plate, '%s' % self.version, 'redmonster-%s-*.fits' % plate)
                for mjdpath in iglob(mjdglobpath):
                    mjd = basename(mjdpath)[16:21]
                    if mjd not in mjds:
                        mjds.append(mjd)
                if len(mjds) > 1:
                    print 'Plate %s has multiple MJDs' % plate
                    hdu1 = fits.open( join( self.redmonster_spectro_redux, plate, self.version, 'redmonster-%s-%s.fits' % (plate,mjds[0]) ) )
                    hdu2 = fits.open( join( self.redmonster_spectro_redux, plate, self.version, 'redmonster-%s-%s.fits' % (plate,mjds[1]) ) )
                    for i,z in enumerate(hdu1[1].data.Z1):
                        if (hdu1[1].data.ZWARNING[i] == 0) & (hdu1[1].data.CLASS1[i] == 'ssp_em_galaxy') & (hdu2[1].data.ZWARNING[i] == 0) & (hdu2[1].data.CLASS1[i] == 'ssp_em_galaxy'):
                            z1.append(z)
                            z2.append(hdu2[1].data.Z1[i])
                            zerr1.append(hdu1[1].data.Z_ERR1[i])
                            zerr2.append(hdu2[1].data.Z_ERR1[i])
        z1 = n.array(z1)
        z2 = n.array(z2)
        zerr1 = n.array(zerr1)
        zerr2 = n.array(zerr2)
        z_diff = z2-z1
        zerr_rms = n.sqrt( (zerr1**2 + zerr2**2)/2. ) # In original paper, this was n.sqrt( (zerr1**2 + zerr2**2) )
        scaled_diff = z_diff / zerr_rms
        hist,binedges = n.histogram(scaled_diff,bins=nbins)
        normhist = hist / float(z1.shape[0])
        bins = n.zeros(nbins)
        for i in xrange(nbins):
            bins[i] = (binedges[i+1]+binedges[i])/2.
        import pdb; pdb.set_trace()
        p.plot(bins, hist, drawstyle='steps-mid', color='black')
        
        def fit_func(x,a,sigma,mu): # Gaussian function to fit to histogram
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
        hdu1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], 'test/bautista/v5_8_guy_split1', '%s' % plate, 'v5_8_guy_split1', 'redmonster-%s-%s.fits' % (plate,mjd)))
        hdu2 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], 'test/bautista/v5_8_guy_split2', '%s' % plate, 'v5_8_guy_split2', 'redmonster-%s-%s.fits' % (plate,mjd)))
        hdu3 = fits.open(join(environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'], '%s' % plate, 'spPlate-%s-%s.fits' % (plate,mjd)))
        for i,ebt1 in enumerate(hdu3[5].data.EBOSS_TARGET1):
            if ebt1 & 2 > 0:
                if True: #(hdu1[1].data.ZWARNING[i] == 0) and (hdu2[1].data.ZWARNING[i] == 0):
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
            self.plate_splits_function(plate=plate, mjd=mjds[i], nbins=nbins, fit=fit)
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
        print n.max(n.abs(scaled_diff))
        print scaled_diff.shape
        hist,binedges = n.histogram(scaled_diff, bins = nbins)
        if normed:
            normhist = hist / float(self.z1.shape[0])
        else:
            normhist = hist
        bins = n.zeros(nbins)
        for i in xrange(nbins):
            bins[i] = (binedges[i+1]+binedges[i])/2.
        p.plot(bins, normhist, drawstyle='steps-mid', color='black')

        def fit_func(x, a, sigma, mu): # Gaussian function to fit to histogram
            return a * n.exp( -((x-mu)**2)/(2.*sigma**2) )

        if fit:
            popt, pcov = curve_fit(fit_func, bins, normhist)
            xfit = n.linspace(-4,4,1000)
            yfit = fit_func(xfit, popt[0], popt[1], popt[2])
            p.plot(xfit, yfit, color='cyan')
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
        globpath1 = join(environ['REDMONSTER_SPECTRO_REDUX'],'v5_8_0_bad', '*')
        for chi2max in xdata:
            print chi2max
            total = 0.
            countidl = 0.
            countrm = 0.
            for path in iglob(globpath1):
                plate = basename(path)
                if len(plate) == 4:
                    globpath2 = join(environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, 'spPlate-%s-*.fits' % plate)
                    for file in iglob(globpath2):
                        if len(basename(file)) == 23:
                            mjd = basename(file)[13:18]
                            hduplate = fits.open(file)
                            hduidl = fits.open(join(environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, '%s' % self.version, 'spZbest-%s-%s.fits' % (plate,mjd)))
                            hdurm = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, '%s' % self.version, 'redmonster-%s-%s.fits' % (plate,mjd)))
                            #for i,ebt0 in enumerate(hduplate[5].data.EBOSS_TARGET0):
                                #if ebt0 == 0:
                            for i,zwarn in enumerate(hdurm[1].data.ZWARNING):
                                if zwarn & 1 > 0:
                                    total += 1.
                                    if hduidl[1].data.RCHI2DIFF[i] > chi2max: countidl += 1.
                                    if hdurm[1].data.RCHI2DIFF[i] > chi2max: countrm += 1.
            rm_ydata.append(countrm/total)
            idl_ydata.append(countidl/total)

        f = p.figure()
        ax = f.add_subplot(1,1,1)
        p.plot(xdata, rm_ydata, drawstyle='steps-mid', color='red', label='redmonster')
        p.plot(xdata, idl_ydata, drawstyle='steps-mid', color='blue', label='idlspec1d')
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
        print '%s galaxy-galaxy confusions of %s, which is %s' % (galgal,total,(galgal/total)*100)
        print '%s galaxy-star confusions of %s, which is %s' % (galstar,total,(galstar/total)*100)
        print '%s galaxy-QSO confusions of %s, which is %s' % (galqso,total,(galqso/total)*100)
        print '%s star-star confusions of %s, which is %s' % (starstar,total,(starstar/total)*100)
        print '%s star-QSO confusions of %s, which is %s' % (starqso,total,(starqso/total)*100)
        print '%s QSO-QSO confusions of %s, which is %s' % (qsoqso,total,(qsoqso/total)*100)


    def rchi2_null_histo(self, nbins=25):
        self.read_redmonster_summary_file()
        rchi2_nulls = self.rm_chi2_null / self.rm_dof
        hist, binedges = n.histogram(rchi2_nulls, bins=nbins)
        bins = n.zeros(nbins)
        for i in xrange(nbins):
            bins[i] = (binedges[i+1] + binedges[i]) / 2.
        p.plot(bins, hist, drawstyle='steps-mid')
        p.savefig('/uufs/astro.utah.edu/common/home/u0814744/boss/rchi2_null_histo.pdf')













# S/N per fiber is located in spZbest files in hdu[1].data.SN_MEDIAN .  You can get just r,i,z bands with x = hdu[1].data.SN_MEDIAN[:,2:] .

# Fiber magnitudes are in spZbest files in hdu[1].data.SPECTROFLUX . Units are nanomaggies, convert to magnitudes
# with 22.5 - 2.5 * LOG_10(SPECTROFLUX)

# To see fibers with zwarning != 0, ztype = 'galaxy', and boss_target1 = 'cmass', use >>> print n.where( (x.rm_zwarning != 0) & (x.rm_type == 'ssp_em_galaxy') & (x.boss_target1 & 2 == 2) )[0]+1

# Plate 7338 has 6 MJDs, 7340 has 4





















