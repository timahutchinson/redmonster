# Write output files after running entirety of redmonster
#
# Tim Hutchinson, University of Utah, August 2014
#
# Edited to reflect changes made to zfitter and zpicker, July 2015
# t.hutchinson@utah.edu

import numpy as n
from astropy.io import fits
from os import environ, makedirs, getcwd, remove
from os.path import exists, join, basename
from astropy.io import fits
from time import gmtime, strftime
from glob import iglob
import re

class Write_Redmonster:
    '''
        Class to write output file at the end of running redmonster.
        
        The zpick argument is simply the entire object created by running
        redmonster.physics.zpicker.py . The self.dest argument is a string
        containing the path in which to save the output file.
        
        If no dest argument is given, or the path does not exist, then the write_rm() method
        will default to writing in $REDMONSTER_SPECTRO_REDUX/$RUN2D/pppp/$RUN1D/ .  If the
        necessary environmental variables are also not specified, it will write in
        directory in which it is being run.
        
        The default behavior is to not clobber any older version of the output file in the
        given directory.  Setting clobber=True will overwrite old versions of the output file.
        '''
    def __init__(self, zpick, dest=None, clobber=True):
        self.clobber = clobber
        self.zpick = zpick
        #if dest and exists(dest): self.dest = dest
        if dest is not None:
            if exists(dest):
                self.dest = dest
            else:
                try:
                    makedirs(dest)
                    self.dest = dest
                except:
                    self.dest = None
        else:
            bsr = environ['REDMONSTER_SPECTRO_REDUX']
            run2d = environ['RUN2D']
            run1d = environ['RUN1D']
            if bsr and run2d and run1d:
                testpath = join(bsr, run2d, '%s' % zpick.plate, run1d)
                if exists(testpath):
                    self.dest = testpath
                else:
                    try:
                        makedirs(testpath)
                        self.dest = testpath
                    except:
                        self.dest = None
            else: self.dest = None

    '''
    def create_hdulist(self):
        # Get old header, append new stuff
        hdr = self.zpick.hdr
        hdr.extend([('VERS_RM','v0.-1','Version of redmonster used'),('DATE_RM',strftime("%Y-%m-%d_%H:%M:%S", gmtime()),'Time of redmonster completion'), ('NFIBERS', len(self.zpick.z), 'Number of fibers'), ('NZ', len(self.zpick.z[0]), 'Number of redshifts retained')])
        prihdu = fits.PrimaryHDU(header=self.zpick.hdr)
        # Columns for 1st BIN table
        colslist = []
        colslist.append( fits.Column(name='FIBERID', format='J', array=self.zpick.fiberid) )
        colslist.append( fits.Column(name='DOF', format='J', array=self.zpick.dof) )
        if hasattr(self.zpick, 'boss_target1'):
            colslist.append( fits.Column(name='BOSS_TARGET1', format='J', array=self.zpick.boss_target1) )
        elif hasattr(self.zpick, 'eboss_target1'):
            colslist.append( fits.Column(name='EBOSS_TARGET1', format='J', array=self.zpick.eboss_target1) )
        for i in xrange(len(self.zpick.z[0])):
            zlist = []
            zerrlist = []
            classlist = []
            subclasslist = []
            minvectorlist = []
            npolylist = []
            fnamelist = []
            npixsteplist = []
            minrchi2list = []
            for j in xrange(len(self.zpick.z)):
                zlist.append( self.zpick.z[j][i] )
                zerrlist.append( self.zpick.z_err[j][i] )
                classlist.append( self.zpick.type[j][i] )
                subclasslist.append( repr(self.zpick.subtype[j][i]) )
                fnamelist.append( self.zpick.fname[j][i] )
                minvectorlist.append( repr(self.zpick.minvector[j][i]) )
                npolylist.append( self.zpick.npoly[j][i] )
                npixsteplist.append( self.zpick.npoly[j][i] )
                minrchi2list.append( self.zpick.minrchi2[j][i] )
            colslist.append( fits.Column(name='Z%s' % (i+1), format='E', array=zlist) )
            colslist.append( fits.Column(name='Z_ERR%s' % (i+1), format='E', array=zerrlist) )
            colslist.append( fits.Column(name='CLASS%s' % (i+1), format='%iA' % max(map(len,classlist)), array=classlist) )
            colslist.append( fits.Column(name='SUBCLASS%s' % (i+1), format='%iA' % max(map(len,subclasslist)), array=subclasslist) )
            colslist.append( fits.Column(name='FNAME%s' % (i+1), format='%iA' % max(map(len,fnamelist)), array=fnamelist) )
            colslist.append( fits.Column(name='MINVECTOR%s' % (i+1), format='%iA' % max(map(len,minvectorlist)), array=minvectorlist) )
            colslist.append( fits.Column(name='MINRCHI2%s' % (i+1), format='E', array=minrchi2list) )
            colslist.append( fits.Column(name='NPOLY%s' % (i+1), format='J', array=npolylist) )
            colslist.append( fits.Column(name='NPIXSTEP%s' % (i+1), format='J', array=npixsteplist) )
        colslist.append( fits.Column(name='ZWARNING', format='J', array=self.zpick.zwarning) )
        colslist.append( fits.Column(name='RCHI2DIFF', format='E', array=self.zpick.rchi2diff) )
        cols = fits.ColDefs(colslist)
        tbhdu = fits.BinTableHDU.from_columns(cols) #tbhdu = fits.new_table(cols)
        # ImageHDU of models
        sechdu = fits.ImageHDU(data=self.zpick.models)
        self.thdulist = fits.HDUList([prihdu, tbhdu, sechdu]) #self.thdulist = fits.HDUList([prihdu, tbhdu])
'''
    
    def create_hdulist(self):
        # Get old header, append new stuff
        hdr = self.zpick.hdr
        hdr.extend([('VERS_RM','v0.-1','Version of redmonster used'),('DATE_RM',strftime("%Y-%m-%d_%H:%M:%S", gmtime()),'Time of redmonster completion'), ('NFIBERS', self.zpick.z.shape[0], 'Number of fibers')])
        prihdu = fits.PrimaryHDU(header=self.zpick.hdr)
        # Columns for 1st BIN table
        col1 = fits.Column(name='Z1', format='E', array=self.zpick.z[:,0])
        col2 = fits.Column(name='Z2', format='E', array=self.zpick.z[:,1])
        '''
            col2_2 = fits.Column(name='Z3', format='E', array=self.zpick.z[:,2])
            col2_3 = fits.Column(name='Z4', format='E', array=self.zpick.z[:,3])
            col2_4 = fits.Column(name='Z5', format='E', array=self.zpick.z[:,4])
            '''
        col3 = fits.Column(name='Z_ERR1', format='E', array=self.zpick.z_err[:,0])
        col4 = fits.Column(name='Z_ERR2', format='E', array=self.zpick.z_err[:,1])
        '''
            col4_2 = fits.Column(name='Z_ERR3', format='E', array=self.zpick.z_err[:,2])
            col4_3 = fits.Column(name='Z_ERR4', format='E', array=self.zpick.z_err[:,3])
            col4_4 = fits.Column(name='Z_ERR5', format='E', array=self.zpick.z_err[:,4])
            '''
        classx = n.array(map(repr,self.zpick.type))
        maxlen = max(map(len,classx))
        col5 = fits.Column(name='CLASS', format='%iA'%maxlen, array=self.zpick.type)
        # Change dictionary values of subclass to strings to be written into fits file.  eval('dictstring') will turn them back into dictionaries later.
        if type(self.zpick.subtype[0]) is not dict: # Skip this and write directly if these are already dicts
            subclass = n.array(map(repr,self.zpick.subtype))
        else:
            subclass = n.array(self.zpick.subtype)
        maxlen = max(map(len,subclass))
        col6 = fits.Column(name='SUBCLASS', format='%iA'%maxlen, array=subclass)
        col7 = fits.Column(name='FIBERID', format='J', array=self.zpick.fiberid)
        if type(self.zpick.minvector[0]) is not str:
            minvec = n.array(map(repr,self.zpick.minvector)) # Change tuples of minvector to strings to be written into fits file. eval('minvector') will turn them back into tuples later.
        else:
            minvec = n.array(self.zpick.minvector)
        maxlen = max(map(len,minvec))
        col8 = fits.Column(name='MINVECTOR', format='%iA'%maxlen, array=minvec)
        col9 = fits.Column(name='ZWARNING', format='E', array=map(int,self.zpick.zwarning))
        col10 = fits.Column(name='DOF', format='E', array=self.zpick.dof)
        col11 = fits.Column(name='NPOLY', format='E', array=self.zpick.npoly)
        fname = n.array(map(repr,self.zpick.fname))
        maxlen = max(map(len,fname))
        col12 = fits.Column(name='FNAME', format='%iA'%maxlen, array=fname)
        col13 = fits.Column(name='NPIXSTEP', format='E', array=self.zpick.npixstep)
        col14 = fits.Column(name='RCHI2DIFF', format='E', array=self.zpick.chi2diff)
        #col14_2 = fits.Column(name='NUMZ', formate='J', array=self.zpick.num_z)
        try:
            col15 = fits.Column(name='BOSS_TARGET1', format='E', array=self.zpick.boss_target1)
        except:
            pass
        try:
            col16 = fits.Column(name='EBOSS_TARGET1', format='E', array=self.zpick.eboss_target1)
        except:
            pass
        if col15:
            if col16:
                cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16]) # TEMPORARY!!!
            else:
                cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15])
        else:
            if col16:
                cols = fits.ColDefs([col1, col2, col2_2, col2_3, col2_4, col3, col4, col4_2, col4_3, col4_4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col16])
            else:
                cols = fits.ColDefs([col1, col2, col2_2, col2_3, col2_4, col3, col4, col4_2, col4_3, col4_4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14])
        tbhdu = fits.BinTableHDU.from_columns(cols) #tbhdu = fits.new_table(cols)
        # ImageHDU of models
        sechdu = fits.ImageHDU(data=self.zpick.models)
        self.thdulist = fits.HDUList([prihdu, tbhdu, sechdu]) #self.thdulist = fits.HDUList([prihdu, tbhdu])

    
    
    def write_fiber(self):
        self.create_hdulist()
        if self.clobber:
            if self.dest is not None:
                self.thdulist.writeto(join(self.dest, '%s' % 'redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0])), clobber=self.clobber)
                print 'Writing redmonster file to %s' % join(self.dest, '%s' % 'redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0]))
            else:
                self.thdulist.writeto('redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0]), clobber=self.clobber)
                print 'Writing redmonster file to %s' % join( getcwd(), 'redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd) )
        else:
            if self.dest is not None:
                if exists(join(self.dest, '%s' % 'redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0]))):
                    self.thdulist.writeto(join(self.dest, '%s' % 'redmonster-%s-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0], strftime("%Y-%m-%d_%H:%M:%S", gmtime()))))
                    print 'Writing redmonster file to %s' % join(self.dest, '%s' % 'redmonster-%s-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0], strftime("%Y-%m-%d_%H:%M:%S", gmtime())))
                else:
                    self.thdulist.writeto(join(self.dest, '%s' % 'redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0])))
                    print 'Writing redmonster file to %s' % join(self.dest, '%s' % 'redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0]))
            else:
                if exists('redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0])):
                    self.thdulist.writeto('redmonster-%s-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0], strftime("%Y-%m-%d_%H:%M:%S", gmtime())))
                    print 'Writing redmonster file to %s' % join( getcwd(), 'redmonster-%s-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0], strftime("%Y-%m-%d_%H:%M:%S", gmtime())))
                else:
                    self.thdulist.writeto('redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0]))
                    print 'Writing redmonster file to %s' % join( getcwd(), 'redmonster-%s-%s-%03d.fits' % (self.zpick.plate, self.zpick.mjd, self.zpick.fiberid[0]))



    def write_plate(self):
        self.create_hdulist()
        
        if self.clobber:
            if self.dest is not None:
                self.thdulist.writeto(join(self.dest, '%s' % 'redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd)), clobber=self.clobber)
                print 'Writing redmonster file to %s' % join(self.dest, '%s' % 'redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd))
            else:
                self.thdulist.writeto('redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd), clobber=self.clobber)
                print 'Writing redmonster file to %s' % join( getcwd(), 'redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd) )
        else:
            if self.dest is not None:
                if exists(join(self.dest, '%s' % 'redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd))):
                    self.thdulist.writeto(join(self.dest, '%s' % 'redmonster-%s-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd, strftime("%Y-%m-%d_%H:%M:%S", gmtime()))))
                    print 'Writing redmonster file to %s' % join(self.dest, '%s' % 'redmonster-%s-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd, strftime("%Y-%m-%d_%H:%M:%S", gmtime())))
                else:
                    self.thdulist.writeto(join(self.dest, '%s' % 'redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd)))
                    print 'Writing redmonster file to %s' % join(self.dest, '%s' % 'redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd))
            else:
                if exists('redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd)):
                    self.thdulist.writeto('redmonster-%s-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd, strftime("%Y-%m-%d_%H:%M:%S", gmtime())))
                    print 'Writing redmonster file to %s' % join( getcwd(), 'redmonster-%s-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd, strftime("%Y-%m-%d_%H:%M:%S", gmtime())))
                else:
                    self.thdulist.writeto('redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd))
                    print 'Writing redmonster file to %s' % join( getcwd(), 'redmonster-%s-%s.fits' % (self.zpick.plate, self.zpick.mjd))


# Combine individual fiber fits files into a single plate file, or combine all plate files into an spAll-like file
# To combine fiber files, create object for a given plate, mjd and call method merge_fibers()
# To create spAll-like file, instantiate with no plate, mjd and call methond merge_plates()
#
# Tim Hutchinson, University of Utah, November 2014
# t.hutchinson@utah.edu

class Merge_Redmonster:
    
    def __init__(self, plate=None, mjd=None, temp=None):
        self.plate = plate
        self.mjd = mjd
        self.temp = temp
    
    
    def merge_fibers(self):
        self.filepaths = []
        self.type = []
        self.subtype = []
        self.fiberid = []
        self.minvector = []
        self.minrchi2 = []
        self.zwarning = []
        self.dof = []
        self.npoly = []
        self.fname = []
        self.npixstep = []
        self.boss_target1 = []
        self.chi2diff = []
        self.models = None
        self.hdr = None
        
        try: topdir = environ['REDMONSTER_SPECTRO_REDUX']
        except: topdir = None
        try: run2d = environ['RUN2D']
        except: run2d = None
        try: run1d = environ['RUN1D']
        except: run1d = None
        fiberdir = join(topdir, run2d, '%s' % self.plate, run1d, 'redmonster-%s-%s-*.fits' % (self.plate, self.mjd)) if topdir and run2d and run1d else None
        
        if fiberdir:
            for path in iglob(fiberdir):
                self.filepaths.append(path)
                fiberfile = basename(path)
                self.fiberid.append( int(fiberfile[22:25]) )
            self.z = n.zeros( (len(self.fiberid),5) )
            self.z_err = n.zeros( self.z.shape )
            try: self.hdr = fits.open( join( environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'], '%s' % self.plate, 'spPlate-%s-%s.fits' % (self.plate,self.mjd) ) )[0].header
            except: self.hdr = fits.Header()
            npix = fits.open( join( environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'], '%s' % self.plate, 'spPlate-%s-%s.fits' % (self.plate,self.mjd) ) )[0].data.shape[1]
            self.models = n.zeros( (self.z.shape[0],npix) )
            self.filepaths.sort()
            self.fiberid.sort()
            for i, path in enumerate(self.filepaths):
                hdu = fits.open(path)
                self.z[i,0] = hdu[1].data.Z1[0]
                self.z[i,1] = hdu[1].data.Z2[0]
                self.z_err[i,0] = hdu[1].data.Z_ERR1[0]
                self.z_err[i,1] = hdu[1].data.Z_ERR2[0]
                self.type.append(hdu[1].data.CLASS1[0])
                self.subtype.append(hdu[1].data.SUBCLASS1[0])
                self.minvector.append(hdu[1].data.MINVECTOR1[0])
                self.zwarning.append(hdu[1].data.ZWARNING[0])
                self.dof.append(hdu[1].data.DOF[0])
                self.npoly.append(hdu[1].data.NPOLY1[0])
                self.fname.append(hdu[1].data.FNAME1[0])
                self.npixstep.append(hdu[1].data.NPIXSTEP1[0])
                self.chi2diff.append(hdu[1].data.RCHI2DIFF[0])
                try:
                    self.boss_target1.append(hdu[1].data.BOSS_TARGET1[0])
                except:
                    pass
                try:
                    self.eboss_target1.append(hdu[1].data.EBOSS_TARGET1[0])
                except:
                    pass
                self.models[i] = hdu[2].data[0]
                #remove(path)
            output = Write_Redmonster(self, clobber=True)
            output.write_plate()


    def merge_plates(self):
        self.type = []
        self.subtype = []
        self.fiberid = []
        self.minvector = []
        self.zwarning = []
        self.dof = []
        self.npoly = []
        self.fname = []
        self.npixstep = []
        self.chi2diff = []
        self.boss_target1 = []
        self.eboss_target1 = []
        self.plates = []
        self.models = n.zeros((1,1))
        self.hdr = fits.Header()
        
        try: topdir = environ['REDMONSTER_SPECTRO_REDUX']
        except: topdir = None
        try: run2d = environ['RUN2D']
        except: run2d = None
        try: run1d = environ['RUN1D']
        except: run1d = None
        platedir = join( topdir, run2d, '*') if topdir and run2d else None
        
        if platedir:
            for path in iglob(platedir):
                self.plates.append( basename(path) )
            self.plates.sort()
            for listitem in self.plates:
                if listitem[-5:] == '.fits': self.plates.remove(listitem)
            self.fiberid = self.plates
            for plate in self.plates:
                print 'Merging plate %s' % plate
                mjds = []
                try:
                    for x in iglob( join( topdir, run2d, '%s' % plate, run1d, 'redmonster-%s-*.fits' % plate) ):
                        if basename(x)[16:21] not in mjds: mjds.append(basename(x)[16:21])
                #if mjds is not basename(x)[16:21]: mjds.append(basename(x)[16:21])
                #else: mjds.append( basename(x)[16:21] )
                except: mjds = None
                if mjds is not [] and mjds is not None:
                    for mjd in mjds:
                        filepath = join( topdir, run2d, str(plate), run1d, 'redmonster-%s-%s.fits' % (plate, mjd))
                        #npix = fits.open( join( environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'], '%s' % plate, 'spPlate-%s-%s.fits' % (plate, mjd) ) )[0].data.shape[1]
                        if exists(filepath):
                            hdu = fits.open(filepath)
                            self.type += hdu[1].data.CLASS.tolist()
                            self.subtype += hdu[1].data.SUBCLASS.tolist()
                            self.minvector += hdu[1].data.MINVECTOR.tolist()
                            self.zwarning += hdu[1].data.ZWARNING.tolist()
                            self.dof += hdu[1].data.DOF.tolist()
                            self.npoly += hdu[1].data.NPOLY.tolist()
                            self.fname += hdu[1].data.FNAME.tolist()
                            self.npixstep += hdu[1].data.NPIXSTEP.tolist()
                            self.chi2diff += hdu[1].data.CHI2DIFF.tolist()
                            try: self.z1 = n.append(self.z1, hdu[1].data.Z1)
                            except: self.z1 = hdu[1].data.Z1
                            try: self.z_err1 = n.append(self.z_err1, hdu[1].data.Z_ERR1)
                            except: self.z_err1 = hdu[1].data.Z_ERR1
                            try: self.z2 = n.append(self.z2, hdu[1].data.Z2)
                            except: self.z2 = hdu[1].data.Z2
                            try: self.z_err2 = n.append(self.z_err2, hdu[1].data.Z_ERR2)
                            except: self.z_err2 = hdu[1].data.Z_ERR2
        self.z = n.zeros( (self.z1.shape[0],2) )
        self.z_err = n.zeros( self.z.shape )
        self.z[:,0] = self.z1
        self.z[:,1] = self.z2
        self.z_err[:,0] = self.z_err1
        self.z_err[:,1] = self.z_err2
        
        output = Write_Redmonster(self)
        output.create_hdulist()
        output.thdulist.writeto( join( topdir, run2d, 'redmonster-all-%s.fits' % run2d), clobber=True)


    def merge_chi2(self):
        
        try: topdir = environ['REDMONSTER_SPECTRO_REDUX']
        except: topdir = None
        try: run2d = environ['RUN2D']
        except: run2d = None
        try: run1d = environ['RUN1D']
        except: run1d = None
        chi2path = join( topdir, run2d, '%s' % self.plate, run1d, 'chi2arr-%s-%s-%s-*.fits' % (self.temp, self.plate, self.mjd) ) if topdir and run2d and run1d else None
        
        fiberid = []
        paths = []
        
        if chi2path:
            for file in iglob(chi2path):
                paths.append( file )
                m = re.search( 'chi2arr-%s-%s-%s-(\d+).fits' % (self.temp, self.plate, self.mjd), basename(file) )
                if m.group(1): fiberid.append( int(m.group(1)) )
            fiberid.sort()
            paths.sort()
            
            for i,path in enumerate(paths):
                chi2arr = fits.open(path)[0].data
                try:
                    chi2arrs
                except NameError:
                    chi2arrs = n.zeros( (len(fiberid),) + chi2arr.shape[1:] )
                    chi2arrs[i] = chi2arr
                else:
                    chi2arrs[i] = chi2arr
                remove(path)

            prihdu = fits.PrimaryHDU(chi2arrs)
            col1 = fits.Column(name='FIBERID', format='J', array=fiberid)
            cols = fits.ColDefs([col1])
            tbhdu = fits.BinTableHDU.from_columns(cols)
            thdulist = fits.HDUList([prihdu,tbhdu])
            thdulist.writeto( join( topdir, run2d, '%s' % self.plate, run1d, 'chi2arr-%s-%s-%s.fits' % (self.temp, self.plate, self.mjd) ), clobber=True)
