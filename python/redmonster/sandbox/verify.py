import numpy as n
from astropy.io import fits
from os.path import join
from os import environ
from glob import iglob
from redmonster.sandbox import yanny as y

class verify_rm:
    
    def __init__(self,version='v5_7_0',plates=[3686,3687,3804,3805,3853,3855,3856,3860],mjds={3686:55268,3687:55269,3804:55267,3805:55269,3853:55268,3855:55268,3856:55269,3860:55269}):
        self.version = version
        self.plates = plates
        self.mjds = mjds
        self.redmonster_spectro_redux = join( environ['REDMONSTER_SPECTRO_REDUX'], '%s' % self.version)
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
        # Read yanny file
        x = y.yanny(filename='/uufs/astro.utah.edu/common/home/u0814744/boss/spInspect_alltest_bolton.par.txt', np=True)
        #x = y.yanny(filename='/Users/boltonlab3/boss/spInspect_alltest_bolton.par.txt', np=True)
        # Get fibers, zpipe, zperson for each plate
        args = n.where(x['BOSSOBJECT']['plate'] == 3686)[0]
        self.fibers3686 = []
        self.zpipe3686 = []
        self.zperson3686 = []
        self.type3686 = []
        for i in args:
            self.fibers3686.append( x['BOSSOBJECT'][i][2])
            self.zpipe3686.append( x['BOSSOBJECT'][i][5])
            self.zperson3686.append( x['BOSSOBJECT'][i][6])
            self.type3686.append( x['BOSSOBJECT'][i][7])
        args = n.where(x['BOSSOBJECT']['plate'] == 3687)[0]
        self.fibers3687 = []
        self.zpipe3687 = []
        self.zperson3687 = []
        self.type3687 = []
        for i in args:
            self.fibers3687.append( x['BOSSOBJECT'][i][2])
            self.zpipe3687.append( x['BOSSOBJECT'][i][5])
            self.zperson3687.append( x['BOSSOBJECT'][i][6])
            self.type3687.append( x['BOSSOBJECT'][i][7])
        args = n.where(x['BOSSOBJECT']['plate'] == 3804)[0]
        self.fibers3804 = []
        self.zpipe3804 = []
        self.zperson3804 = []
        self.type3804 = []
        for i in args:
            self.fibers3804.append( x['BOSSOBJECT'][i][2])
            self.zpipe3804.append( x['BOSSOBJECT'][i][5])
            self.zperson3804.append( x['BOSSOBJECT'][i][6])
            self.type3804.append( x['BOSSOBJECT'][i][7])
        args = n.where(x['BOSSOBJECT']['plate'] == 3805)[0]
        self.fibers3805 = []
        self.zpipe3805 = []
        self.zperson3805 = []
        self.type3805 = []
        for i in args:
            self.fibers3805.append( x['BOSSOBJECT'][i][2])
            self.zpipe3805.append( x['BOSSOBJECT'][i][5])
            self.zperson3805.append( x['BOSSOBJECT'][i][6])
            self.type3805.append( x['BOSSOBJECT'][i][7])
        args = n.where(x['BOSSOBJECT']['plate'] == 3853)[0]
        self.fibers3853 = []
        self.zpipe3853 = []
        self.zperson3853 = []
        self.type3853 = []
        for i in args:
            self.fibers3853.append( x['BOSSOBJECT'][i][2])
            self.zpipe3853.append( x['BOSSOBJECT'][i][5])
            self.zperson3853.append( x['BOSSOBJECT'][i][6])
            self.type3853.append( x['BOSSOBJECT'][i][7])
        args = n.where(x['BOSSOBJECT']['plate'] == 3855)[0]
        self.fibers3855 = []
        self.zpipe3855 = []
        self.zperson3855 = []
        self.type3855 = []
        for i in args:
            self.fibers3855.append( x['BOSSOBJECT'][i][2])
            self.zpipe3855.append( x['BOSSOBJECT'][i][5])
            self.zperson3855.append( x['BOSSOBJECT'][i][6])
            self.type3855.append( x['BOSSOBJECT'][i][7])
        args = n.where(x['BOSSOBJECT']['plate'] == 3856)[0]
        self.fibers3856 = []
        self.zpipe3856 = []
        self.zperson3856 = []
        self.type3856 = []
        for i in args:
            self.fibers3856.append( x['BOSSOBJECT'][i][2])
            self.zpipe3856.append( x['BOSSOBJECT'][i][5])
            self.zperson3856.append( x['BOSSOBJECT'][i][6])
            self.type3856.append( x['BOSSOBJECT'][i][7])
        args = n.where(x['BOSSOBJECT']['plate'] == 3860)[0]
        self.fibers3860 = []
        self.zpipe3860 = []
        self.zperson3860 = []
        self.type3860 = []
        for i in args:
            self.fibers3860.append( x['BOSSOBJECT'][i][2])
            self.zpipe3860.append( x['BOSSOBJECT'][i][5])
            self.zperson3860.append( x['BOSSOBJECT'][i][6])
            self.type3860.append( x['BOSSOBJECT'][i][7])

    def read_redmonster(self,plate):
        redmonsterpath = join( self.redmonster_spectro_redux, '%s' % plate, '%s' % self.version, 'redmonster-%s-%s.fits' % (plate,self.mjds[plate]) )
        #for path in iglob(platepath):
        #    hdu = fits.open(path)
        hdu = fits.open(redmonsterpath)
        self.rm_z1 = hdu[1].data.Z1
        self.rm_zerr1 = hdu[1].data.Z_ERR1
        self.rm_fibers = hdu[1].data.FIBERID + 1 # +1 here because rm fibers are 0-based and idlspec2d are 1-based
        self.rm_type = hdu[1].data.CLASS
        self.rm_zwarning = hdu[1].data.ZWARNING
    
    def read_spPlate(self,plate):
        spPlatepath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, 'spPlate-%s-%s.fits' % (plate, self.mjds[plate]) )
        hdu = fits.open(spPlatepath)
        self.boss_target1 = hdu[5].data.BOSS_TARGET1
    
    def read_spZbest(self,plate):
        spZbestpath = join( environ['BOSS_SPECTRO_REDUX'], '%s' % self.version, '%s' % plate, '%s' % self.version, 'spZall-%s-%s.fits' % (plate, self.mjds[plate]) )
        hdu = fits.open(spZbestpath)
        self.sn_median = hdu[1].data.SN_MEDIAN[:,2:]

    def get_cmass(self):
        return n.where( self.boss_target1 & 2 == 2 )[0].tolist()

    def get_lowz(self):
        return n.where( self.boss_target1 & 1 == 1 )[0].tolist()

    def cmass_completeness(self):
        vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            fibers = self.get_cmass()
            vals.append( float(len(n.where( self.rm_zwarning[fibers] == 0 )[0].tolist())) / float(len(fibers)) )
        avg = n.sum(vals) / float(len(vals))
        return avg

    def lowz_completeness(self):
        vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            fibers = self.get_lowz()
            vals.append( float(len(n.where( self.rm_zwarning[fibers] == 0 )[0].tolist())) / float(len(fibers)) )
        avg = n.sum(vals) / float(len(vals))
        return avg

    def cmass_galaxy_completeness(self):
        vals = []
        for plate in self.plates:
            self.read_redmonster(plate)
            self.read_spPlate(plate)
            fibers = self.get_cmass()
            #fibers = n.where(self.rm_type[fibers] == 'ssp_em_galaxy')[0].tolist() # SHOULD THIS USE AS THE TOTAL NUMBER OF FIBERS THOSE CLASSIFIED AS GALAXIES BY RM OR BY VISUAL INSPECTION?
            #vals.append( float(len(n.where( (self.rm_zwarning[fibers] == 0) )[0].tolist())) / float(len(fibers)) )
            vals.append( float(len(n.where((self.rm_zwarning[fibers] == 0) & (self.rm_type[fibers] == 'ssp_em_galaxy'))[0].tolist())) / float(len(n.where(self.rm_type[fibers] == 'ssp_em_galaxy')[0].tolist())) )
        avg = n.sum(vals) / float(len(vals))
        return avg




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

# S/N per fiber is located in spZbest files in hdu[1].data.SN_MEDIAN.  You can get just r,i,z bands with x = hdu[1].data.SN_MEDIAN[:,2:] .





















