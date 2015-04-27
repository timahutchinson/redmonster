import numpy as n
from astropy.io import fits
from os.path import join
from os import environ
from glob import iglob
from redmonster.sandbox import yanny as y

class verify_rm:

    def __init__(self):
        self.plates = [3686,3687,3804,3805,3853,3855,3856,3860]
        self.redmonster_dir = '/uufs/chpc.utah.edu/common/home/bolton_data0/redmonster/v5_5_12/'
        self.yanny_to_arrays()
        self.rm_z = []
        self.vis_z = []
        for i,plate in enumerate(plates):
            if i == 0: self.compare_redshifts(plate,self.fibers3686,self.zperson3686)
            elif i == 1: self.compare_redshifts(plate,self.fibers3687,self.zperson3687)
            elif i == 2: self.compare_redshifts(plate,self.fibers3804,self.zperson3804)
            elif i == 3: self.compare_redshifts(plate,self.fibers3805,self.zperson3805)
            elif i == 4: self.compare_redshifts(plate,self.fibers3853,self.zperson3853)
            elif i == 5: self.compare_redshifts(plate,self.fibers3855,self.zperson3855)
            elif i == 6: self.compare_redshifts(plate,self.fibers3856,self.zperson3856)
            elif i == 7: self.compare_redshifts(plate,self.fibers3860,self.zperson3860)


    def yanny_to_arrays(self):
        # Read yanny file
        x = y.yanny(filename='/uufs/astro.utah.edu/common/home/u0814744/boss/spInspect_alltest_bolton.par.txt', np=True)
        # Get fibers, zpipe, zperson for each plate
        args = n.where(x['BOSSOBJECT']['plate'] == 3686)[0]
        self.fibers3686 = []
        self.zpipe3686 = []
        self.zperson3686 = []
        for i in args:
            self.fibers3686.append( x['BOSSOBJECT'][i][2])
            self.zpipe3686.append( x['BOSSOBJECT'][i][5])
            self.zperson3686.append( x['BOSSOBJECT'][i][6])
        args = n.where(x['BOSSOBJECT']['plate'] == 3687)[0]
        self.fibers3687 = []
        self.zpipe3687 = []
        self.zperson3687 = []
        for i in args:
            self.fibers3687.append( x['BOSSOBJECT'][i][2])
            self.zpipe3687.append( x['BOSSOBJECT'][i][5])
            self.zperson3687.append( x['BOSSOBJECT'][i][6])
        args = n.where(x['BOSSOBJECT']['plate'] == 3804)[0]
        self.fibers3804 = []
        self.zpipe3804 = []
        self.zperson3804 = []
        for i in args:
            self.fibers3804.append( x['BOSSOBJECT'][i][2])
            self.zpipe3804.append( x['BOSSOBJECT'][i][5])
            self.zperson3804.append( x['BOSSOBJECT'][i][6])
        args = n.where(x['BOSSOBJECT']['plate'] == 3805)[0]
        self.fibers3805 = []
        self.zpipe3805 = []
        self.zperson3805 = []
        for i in args:
            self.fibers3805.append( x['BOSSOBJECT'][i][2])
            self.zpipe3805.append( x['BOSSOBJECT'][i][5])
            self.zperson3805.append( x['BOSSOBJECT'][i][6])
        args = n.where(x['BOSSOBJECT']['plate'] == 3853)[0]
        self.fibers3853 = []
        self.zpipe3853 = []
        self.zperson3853 = []
        for i in args:
            self.fibers3853.append( x['BOSSOBJECT'][i][2])
            self.zpipe3853.append( x['BOSSOBJECT'][i][5])
            self.zperson3853.append( x['BOSSOBJECT'][i][6])
        args = n.where(x['BOSSOBJECT']['plate'] == 3855)[0]
        self.fibers3855 = []
        self.zpipe3855 = []
        self.zperson3855 = []
        for i in args:
            self.fibers3855.append( x['BOSSOBJECT'][i][2])
            self.zpipe3855.append( x['BOSSOBJECT'][i][5])
            self.zperson3855.append( x['BOSSOBJECT'][i][6])
        args = n.where(x['BOSSOBJECT']['plate'] == 3856)[0]
        self.fibers3856 = []
        self.zpipe3856 = []
        self.zperson3856 = []
        for i in args:
            self.fibers3856.append( x['BOSSOBJECT'][i][2])
            self.zpipe3856.append( x['BOSSOBJECT'][i][5])
            self.zperson3856.append( x['BOSSOBJECT'][i][6])
        args = n.where(x['BOSSOBJECT']['plate'] == 3860)[0]
        self.fibers3860 = []
        self.zpipe3860 = []
        self.zperson3860 = []
        for i in args:
            self.fibers3860.append( x['BOSSOBJECT'][i][2])
            self.zpipe3860.append( x['BOSSOBJECT'][i][5])
            self.zperson3860.append( x['BOSSOBJECT'][i][6])

    def read_redmonster(self,plate):
        platepath = join( self.redmonster_dir, plate, 'v5_5_12', 'redmonster*' )
        for path in iglob(platepath):
            hdu = fits.open(path)
        self.rm_z1 = hdu[1].data.Z1
        self.rm_zerr1 = hdu[1].data.Z_ERR1
        self.rm_fibers = hdu[1].data.Z2 + 1 # +1 here because rm fibers are 0-based and idlspec2d are 1-based
        self.rm_type = hdu[1].data.CLASS

    def compare_redshifts(self,plate,visual_fibers,visual_z):
        self.read_redmonster(plate)
        for i,fiber in enumerate(visual_fibers):
            rm_ind = n.where(self.rm_fibers == fiber)[0]
            if rm_ind != []:
                rm_ind = rm_ind[0]
                self.rm_z.append(self.rm_z1[rm_ind])
                self.vis_z.append(visual_z[i])




















