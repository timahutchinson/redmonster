# Top level script to run redmonster as a whole.
#
# Required arguments:
# plate (int): 4-digit plate number
# mjd (int): 4-digit mjd for given plate
# templates (list of strings): List of template names, as strings, to be run, up to a max of 5.  Template files must be in $REDMONSTER_DIR/templates/
#
# Optional arguments:
# fiberid (list of ints): List of fibers, as ints, to run.  If not specified, defaults to all fibers.
# zmin, zmax (list): Lists of minimum and maximum redshift for each template.  zmin and zmax must either be specified for
#                    all templates or not at all.  If not given, all templates will be run over entire possible range.
# npoly (list of ints): List of number of polynomial terms for each template.  If not specified, defaults to 4 for all templates.
#               If given, len(npoly) must equal len(templates)
# npixstep (list of ints): Size of pixel steps to use when doing cross-correlation.  If not given, defaults to 1 for all templates.
#                  If given, len(npixstep) must equal len(templates)
# dest (string): Full path to directory in which output fill will be written.  If not specified, defaults to
#                $REDMONSTER_SPECTRO_REDUX/$RUN2D/pppp/$RUN1D/ , where pppp is the 4 digit plate id.
# clobber (Boolean): Default behavior is to overwrite older output files for same plate/mjd.  Setting to false will cause new
#                    version to be written.
#
# Tim Hutchinson, University of Utah, October 2014
# t.hutchinson@utah.edu

import numpy as n
from astropy.io import fits
from redmonster.datamgr import spec, io
from redmonster.physics import zfinder, zfitter, zpicker
from redmonster.physics import misc
from time import gmtime, strftime
from os.path import exists
from ConfigParser import SafeConfigParser
import sys

class Zfind:

    def __init__(self, inifile=None, dest=None, clobber=False):
        self.inifile = inifile
        self.dest = dest
        self.clobber = clobber
        if self.inifile: self.set_templates_from_inifile()

    def set_templates_from_inifile(self):
        self.labels = []
        self.templates = []
        self.zmin = []
        self.zmax = []
        self.npoly = []
        self.npixstep = []
        
        if exists(self.inifile):
            self.option = SafeConfigParser()
            self.option.optionxform = unicode
            r = self.option.read(self.inifile)
            if len(r) == 1:
                for section in self.option.sections():
                    print section
                    print "template=%r" % self.option.get(section,'template')
                    print "zmin=%r" % self.option.get(section,'zmin')
                    
                    self.labels.append(section)
                    self.templates.append(self.option.get(section,'template'))
                    self.zmin.append(self.option.getfloat(section,'zmin'))
                    self.zmax.append(self.option.getfloat(section,'zmax'))
                    self.npoly.append(self.option.getint(section,'npoly'))
                    self.npixstep.append(self.option.getint(section,'npixstep'))
            self.set_templates()
        else: print "WARNING: %r does not exist" % self.inifile
                        
    def set_templates(self, templates=None, zmin=None, zmax=None, npoly=None, npixstep=None):
        if templates: self.templates = templates
        if zmin: self.zmin = zmin
        if zmax: self.zmax = zmax
        if npoly: self.npoly = npoly
        if npixstep: self.npixstep = npixstep
        
        if type(self.templates) is str:
            try: self.templates = [self.templates]
            except:
                print 'Templates not a list and unable to convert to list!'
                sys.exit(1)
        if type(self.templates) is list: self.templates = map(str, self.templates)
        if self.zmin is not None:
            if type(self.zmin) is not list:
                try:
                    self.zmin = [self.zmin]
                except:
                    try:
                        self.zmin = self.zmin.tolist()
                    except:
                        print 'Can\'t convert zmin to list - defaulting to full zrange!'
                        self.zmin = None
                        self.zmax = None
            if type(self.zmin) is list:
                if len(self.zmin) != len(self.templates):
                    print 'Length of zmin doesn\'t match length of templates - defaulting to full zrange!'
                    self.zmin = None
                    self.zmax = None
                if self.zmax is None:
                    print 'zmax not given - defaulting to full zrange!'
                    self.zmin = None
                    self.zmax = None
                else:
                    if type(self.zmax) is not list:
                        try: self.zmax = [self.zmax]
                        except:
                            try:
                                self.zmax = self.zmax.tolist()
                            except:
                                print 'Can\'t convert zmax to list - defaulting to full zrange!'
                                self.zmin = None
                                self.zmax = None
                    if len(self.zmin) != len(self.zmax):
                        print 'Length of zmin and zmax don\'t match - defaulting to full zrange!'
                        self.zmin = None
                        self.zmax = None
        #import pdb; pdb.set_trace()
        if self.npoly is None:
            self.npoly = [4]*len(self.templates)
        else:
            if type(self.npoly) is not list:
                try:
                    self.npoly = [self.npoly]
                except:
                    try:
                        self.npoly = self.npoly.tolist()
                    except:
                        print 'npoly not a list and unable to convert to list - defaulting to npoly=4 for all templates!'
                        self.npoly = [4]*len(self.templates)
            else:
                self.npoly = map(int, self.npoly)
        if self.npixstep is None:
            self.npixstep = [1]*len(self.templates)
        else:
            if type(self.npixstep) is not list:
                try:
                    self.npixstep = [self.npixstep]
                except:
                    try:
                        self.npixstep = self.npixstep.tolist()
                    except:
                        print 'npixstep not a list and unable to convert to list - defaulting to npixstep=1 for all templates!'
                        self.npixstep = [1]*len(self.templates)
            else:
                self.npixstep = map(int, self.npixstep)

    def reduce_plate_mjd(self, plate, mjd, fiberid=None):
        # Check types and try to convert to proper types if necessary
        if fiberid is None: fiberid = [i for i in xrange(1000)]
        else:
            if type(fiberid) is not list:
                try:
                    fiberid = [fiberid]
                    fiberid = map(int, fiberid)
                except:
                    try:
                        fiberid = fiberid.tolist()
                        fiberid = map(int, fiberid)
                    except:
                        print 'fiberid not set properly - running full plate!'
                        fiberid = [i for i in xrange(1000)]
            else: fiberid = map(int, fiberid)

        # Spec
        specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)

        # Zfinder, Zfitter
        zfindobjs = []
        zfitobjs = []
        if (self.zmin is not None) & (self.zmax is not None):
            for i in xrange(len(self.templates)):
                zfindobjs.append( zfinder.Zfinder(fname=self.templates[i], npoly=self.npoly[i], zmin=self.zmin[i], zmax=self.zmax[i]) )
                zfindobjs[i].zchi2( specs.flux, specs.loglambda, specs.ivar, npixstep=self.npixstep[i] )
                zfitobjs.append( zfitter.Zfitter(zfindobjs[i].zchi2arr, zfindobjs[i].zbase) )
                zfitobjs[i].z_refine()
        else:
            for i in xrange(len(self.templates)):
                zfindobjs.append( zfinder.Zfinder(fname=self.templates[i], npoly=self.npoly[i], npixstep=self.npixstep[i]) )
                zfindobjs[i].zchi2( specs.flux, specs.loglambda, specs.ivar, npixstep=self.npixstep[i] )
                zfitobjs.append( zfitter.Zfitter(zfindobjs[i].zchi2arr, zfindobjs[i].zbase) )
                zfitobjs[i].z_refine()

        # Flags
        flags = []
        for i in xrange(len(zfindobjs)):
            flags.append( misc.comb_flags(specs, zfindobjs[i], zfitobjs[i]) )

        # Zpicker
        if len(self.templates) == 1: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0])
        elif len(self.templates) == 2: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0], zfindobjs[1], zfitobjs[1], flags[1])
        elif len(self.templates) == 3: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0], zfindobjs[1], zfitobjs[1], flags[1],
                                                          zfindobjs[2], zfitobjs[2], flags[2])
        elif len(self.templates) == 4: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0], zfindobjs[1], zfitobjs[1], flags[1],
                                                          zfindobjs[2], zfitobjs[2], flags[2], zfindobjs[3], zfitobjs[3], flags[3])
        elif len(self.templates) == 5: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0], zfindobjs[1], zfitobjs[1], flags[1],
                                                          zfindobjs[2], zfitobjs[2], flags[2], zfindobjs[3], zfitobjs[3], flags[3],
                                                          zfindobjs[4], zfitobjs[4], flags[4])

        # Write output
        if self.dest is None:
            output = io.Write_Redmonster(zpick, clobber=self.clobber)
        else:
            if type(self.dest) is str:
                output = io.Write_Redmonster(zpick, dest=self.dest, clobber=self.clobber)
            else:
                try:
                    self.dest = str(self.dest)
                    output = io.Write_Redmonster(zpick, dest=self.dest, clobber=self.clobber)
                except:
                    print 'Could not convert dest to string - writing to default directory and NOT clobbering old files!'
                    output = io.Write_Redmonster(zpick, clobber=False)




        












