import numpy as n
from astropy.io import fits
from os.path import join
from os import environ
from glob import iglob
from redmonster.sandbox import yanny as y
import matplotlib.pyplot as p

class verify_rm:
    
    def __init__(self,version='v5_7_0',plates=[3686,3687,3804,3805,3853,3855,3856,3860],mjds={3686:55268,3687:55269,3804:55267,3805:55269,3853:55268,3855:55268,3856:55269,3860:55269}):
        self.version = version
        self.plates = plates
        self.mjds = mjds
        self.redmonster_spectro_redux = join( environ['REDMONSTER_SPECTRO_REDUX'], '%s' % self.version)
        self.vifibers = None
        self.zperson = None
        self.zpipe = None
        self.vitype = None
        self.comments = None
        self.yanny_to_arrays()
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
        vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_cmass()
            #vals.append( float(len(n.where((self.rm_zwarning[fibers] == 0) & (self.rm_type[fibers] == 'ssp_em_galaxy'))[0].tolist())) / float(len(n.where(self.rm_type[fibers] == 'ssp_em_galaxy')[0].tolist())) )
            vals.append( float(len(n.where((self.rm_zwarning[fibers] == 0) & (self.rm_type[fibers] == 'ssp_em_galaxy'))[0].tolist())) / float(len(fibers)) )
        avg = n.sum(vals) / float(len(vals))
        print avg


    def lowz_galaxy_completeness(self):
        # Prints percent of all LOWZ targets that have rm_zwarning == 0 and were classified as 'ssp_em_galaxy'
        vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            self.read_spZbest(plate)
            fibers = self.get_lowz()
            #vals.append( float(len(n.where((self.rm_zwarning[fibers] == 0) & (self.rm_type[fibers] == 'ssp_em_galaxy'))[0].tolist())) / float(len(n.where(self.rm_type[fibers] == 'ssp_em_galaxy')[0].tolist())) )
            vals.append( float(len(n.where((self.rm_zwarning[fibers] == 0) & (self.rm_type[fibers] == 'ssp_em_galaxy'))[0].tolist())) / float(len(fibers)) )
        avg = n.sum(vals) / float(len(vals))
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
        # Prints fraction of CMASS targets having yanny comment 'v5_4_9 ok' that have rm_zwarning == 0
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


    def dz_to_dv(self, dz):
        # Convert redshift error dz to velocity error dv
        c_kms = 299792.458 # speed of light in km s^-1
        return dz * c_kms


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
                        if n.abs(self.rm_z1[fiber] - self.zperson[vi_index]) >= 0.005:
                            self.bad_plates.append(plate)
                            self.bad_fibers.append(fiber)
                            self.bad_rm_z.append(self.rm_z1[fiber])
                            self.bad_zperson.append(self.zperson[vi_index])
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
        p.savefig('recov.pdf')


    def cmass_failure_vs_sn(self,nbins=25):
        # Makes plot of CMASS failure rate (zwarning > 0) vs median S/N in r-, i-, and z-bands
        f = p.figure()
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
                if self.sn_median[fiber,2] <= 4:
                    total += 1.
                    r_sn.append(self.sn_median[fiber,0])
                    i_sn.append(self.sn_median[fiber,1])
                    z_sn.append(self.sn_median[fiber,2])
                    if (self.rm_zwarning[fiber] > 0):
                        bad_fibers.append(fiber)
                        bad_r_sn.append(self.sn_median[fiber,0])
                        bad_i_sn.append(self.sn_median[fiber,1])
                        bad_z_sn.append(self.sn_median[fiber,2])
        rtotal,junk = n.histogram(r_sn,bins=nbins)
        itotal,junk = n.histogram(i_sn,bins=nbins)
        ztotal,junk = n.histogram(z_sn,bins=nbins)
        rhist,rbinedges = n.histogram(bad_r_sn,bins=nbins)
        ihist,ibinedges = n.histogram(bad_i_sn,bins=nbins)
        zhist,zbinedges = n.histogram(bad_z_sn,bins=nbins)
        rbins = n.zeros(nbins)
        ibins = n.zeros(nbins)
        zbins = n.zeros(nbins)
        #import pdb; pdb.set_trace()
        for i in xrange(nbins):
            rbins[i] = (rbinedges[i+1]+rbinedges[i])/2.
            ibins[i] = (rbinedges[i+1]+rbinedges[i])/2.
            zbins[i] = (rbinedges[i+1]+rbinedges[i])/2.
        rhist = rhist / map(float,rtotal)
        ihist = ihist / map(float,itotal)
        zhist = zhist / map(float,ztotal)
        #p.plot(rbins,rhist,color='purple',label='r-band')
        #p.plot(ibins,ihist,color='blue',label='i-band')
        #p.plot(zbins,zhist,color='cyan',label='z-band')
        p.plot(rbins,rtotal,label='r')
        p.plot(ibins,itotal,label='i')
        p.plot(zbins,ztotal,label='z')
        print zbins
        print ztotal
        print n.sort(z_sn)
        p.legend()
        p.savefig('failure_vs_sn.pdf')












'''
    def compare_redshifts(self,plate,visual_fibers,visual_z): #visual_fibers is list of fiber numbers visually inspected for a given plate, and visual_z is list of visually determined redshifts for those fibers
        self.read_redmonster(plate)
        self.read_spPlate(plate)
        for i,fiber in enumerate(visual_fibers):
            rm_ind = n.where(self.rm_fibers == fiber)[0]
            try:
                rm_ind = rm_ind[0]
                self.rm_z.append(self.rm_z1[rm_ind])
                self.rm_class.append(self.rm_type[rm_ind])
                self.vis_z.append(visual_z[i])
                self.rm_zwarning.append(self.rm_zwarning[rm_ind])
            except:
                pass
'''

# S/N per fiber is located in spZbest files in hdu[1].data.SN_MEDIAN .  You can get just r,i,z bands with x = hdu[1].data.SN_MEDIAN[:,2:] .

# Fiber magnitudes are in spZbest files in hdu[1].data.SPECTROFLUX . Units are nanomaggies, convert to magnitudes
# with 22.5 - 2.5 * LOG_10(SPECTROFLUX)

# To see fibers with zwarning != 0, ztype = 'galaxy', and boss_target1 = 'cmass', use >>> print n.where( (x.rm_zwarning != 0) & (x.rm_type == 'ssp_em_galaxy') & (x.boss_target1 & 2 == 2) )[0]+1





















