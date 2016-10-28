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

from os.path import exists
from time import gmtime, strftime
import sys
try:
    from configparser import SafeConfigParser
except ImportError:
    from ConfigParser import SafeConfigParser

import numpy as n
from astropy.io import fits

from redmonster.datamgr import spec, io, io2
from redmonster.physics import zfinder, zfitter, zpicker, zpicker2
from redmonster.physics import misc


class ZFind:

    def __init__(self, num_z=5, inifile=None, dest=None, nproc=1, clobber=True):
        self.num_z = num_z
        self.inifile = inifile
        self.dest = dest
        self.clobber = clobber
        if self.inifile: self.set_templates_from_inifile()
        self.nproc = nproc

    def set_templates_from_inifile(self):
        self.labels = []
        self.templates = []
        self.zmin = []
        self.zmax = []
        self.npoly = []
        self.npixstep = []
        self.group = []

        if exists(self.inifile):
            self.option = SafeConfigParser()
            self.option.optionxform = str
            r = self.option.read(self.inifile)
            if len(r) == 1:
                for section in self.option.sections():
                    self.labels.append(section)
                    if self.option.has_option(section,'template'):
                        self.templates.append(self.option.get(section,
                                                              'template'))
                    if self.option.has_option(section,'zmin'):
                        self.zmin.append(self.option.getfloat(section,'zmin'))
                    if self.option.has_option(section,'zmax'):
                        self.zmax.append(self.option.getfloat(section,'zmax'))
                    if self.option.has_option(section,'npoly'):
                        self.npoly.append(self.option.getint(section,'npoly'))
                    if self.option.has_option(section,'npixstep'):
                        self.npixstep.append(self.option.getint(section,
                                                                'npixstep'))
                    if self.option.has_option(section,'group'):
                        self.group.append(self.option.get(section,'group'))

            else: print("Cannot parse ini file %r" % self.inifile)
            if not self.labels: self.labels = None
            if not self.templates: self.templates = None
            if not self.zmin: self.zmin = None
            if not self.zmax: self.zmax = None
            if not self.npoly: self.npoly = None
            if not self.npixstep: self.npixstep = None
            if not self.group: self.group = None

            self.set_templates()
        else: print("WARNING: %r does not exist" % self.inifile)


    def set_templates(self, templates=None, zmin=None, zmax=None, npoly=None,
                      npixstep=None, group=None):
        if templates: self.templates = templates
        if zmin: self.zmin = zmin
        if zmax: self.zmax = zmax
        if npoly: self.npoly = npoly
        if npixstep: self.npixstep = npixstep
        if group: self.group = group

        if type(self.templates) is str:
            try: self.templates = [self.templates]
            except Exception as e:
                print('Templates not a list and unable to convert to list! \
                Exception: %r' % e)
                sys.exit(1)
        if type(self.templates) is list:
                self.templates = list(map(str, self.templates))
        if self.zmin is not None:
            if type(self.zmin) is not list:
                try:
                    self.zmin = [self.zmin]
                except:
                    try:
                        self.zmin = self.zmin.tolist()
                    except Exception as e:
                        print('Can\'t convert zmin to list - defaulting to \
                                full zrange! Exception: %r' % e)
                        self.zmin = None
                        self.zmax = None
            if type(self.zmin) is list:
                if len(self.zmin) != len(self.templates):
                    print('Length of zmin doesn\'t match length of templates - \
                            defaulting to full zrange!')
                    self.zmin = None
                    self.zmax = None
                if self.zmax is None:
                    print('zmax not given - defaulting to full zrange!')
                    self.zmin = None
                    self.zmax = None
                else:
                    if type(self.zmax) is not list:
                        try: self.zmax = [self.zmax]
                        except:
                            try:
                                self.zmax = self.zmax.tolist()
                            except Exception as e:
                                print('Can\'t convert zmax to list - \
                                        defaulting to full zrange! \
                                        Exception: %r' % e)
                                self.zmin = None
                                self.zmax = None
                    if len(self.zmin) != len(self.zmax):
                        print('Length of zmin and zmax don\'t match - \
                                defaulting to full zrange!')
                        self.zmin = None
                        self.zmax = None
        #import pdb; pdb.set_trace()
        if self.npoly is None:
            self.npoly = [4]*len(self.templates)
        if self.group is None:
            self.group = [0]*len(self.templates)
        else:
            if type(self.npoly) is not list:
                try:
                    self.npoly = [self.npoly]
                except:
                    try:
                        self.npoly = self.npoly.tolist()
                    except Exception as e:
                        print('npoly not a list and unable to convert to \
                                list - defaulting to npoly=4 for all \
                                templates! Exception: %r' % e)
                        self.npoly = [4]*len(self.templates)
            else:
                self.npoly = list(map(int, self.npoly))
        if self.npixstep is None:
            self.npixstep = [1]*len(self.templates)
        else:
            if type(self.npixstep) is not list:
                try:
                    self.npixstep = [self.npixstep]
                except:
                    try:
                        self.npixstep = self.npixstep.tolist()
                    except Exception as e:
                        print('npixstep not a list and unable to convert to \
                                list - defaulting to npixstep=1 for all \
                                templates! Exception: %r' % e)
                        self.npixstep = [1]*len(self.templates)
            else:
                self.npixstep = list(map(int, self.npixstep))


    def reduce_plate_mjd(self, plate=None, mjd=None, fiberid=None, data_range=None,
                         chi2file=False, platepath=None):
        self.chi2file = chi2file
        # Check types and try to convert to proper types if necessary
        if fiberid is not None:
            if type(fiberid) is not list:
                try:
                    fiberid = [fiberid]
                    fiberid = list(map(int, fiberid))
                except ValueError:
                    try:
                        fiberid = fiberid.tolist()
                        fiberid = list(map(int, fiberid))
                    except ValueError:
                        print('fiberid not set properly - running full plate!')
            else: fiberid = list(map(int, fiberid))

        # Spec
        specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid, platepath=platepath)
        fiberid = specs.fiberid

        # ZFinder, ZFitter
        zfindobjs = []
        zfitobjs = []
        if (self.zmin is not None) & (self.zmax is not None):
            for i in range(len(self.templates)):
                zfindobjs.append( zfinder.ZFinder(fname=self.templates[i],
                                                  group=self.group[i],
                                                  npoly=self.npoly[i],
                                                  zmin=self.zmin[i],
                                                  zmax=self.zmax[i],
                                                  nproc=self.nproc) )
                zfindobjs[i].zchi2( specs.flux, specs.loglambda, specs.ivar,
                                   npixstep=self.npixstep[i], plate=plate,
                                   mjd=mjd, fiberid=fiberid[0],
                                   chi2file=self.chi2file )
                zfitobjs.append( zfitter.ZFitter(zfindobjs[i].zchi2arr,
                                                 zfindobjs[i].zbase) )
                zfitobjs[i].z_refine2()
        else:
            for i in range(len(self.templates)):
                zfindobjs.append( zfinder.ZFinder(fname=self.templates[i],
                                                  group=self.group[i],
                                                  npoly=self.npoly[i],
                                                  npixstep=self.npixstep[i],
                                                  nproc=self.nproc) )
                zfindobjs[i].zchi2( specs.flux, specs.loglambda, specs.ivar,
                                   npixstep=self.npixstep[i], plate=plate,
                                   mjd=mjd, fiberid=fiberid[0],
                                   chi2file=self.chi2file )
                zfitobjs.append( zfitter.ZFitter(zfindobjs[i].zchi2arr,
                                                 zfindobjs[i].zbase) )
                zfitobjs[i].z_refine2()

        # Flags
        flags = []
        for i in range(len(zfindobjs)):
            flags.append( misc.comb_flags(specs, zfindobjs[i], zfitobjs[i]) )

        # ZPicker
        zpick = zpicker2.ZPicker(specs, zfindobjs, zfitobjs, flags)



        output = None

        # Write output
        if self.dest is None:
            output = io2.WriteRedmonster(zpick, clobber=True)
        else:
            if type(self.dest) is str:
                output = io2.WriteRedmonster(zpick, dest=self.dest,
                                              clobber=True)
            else:
                try:
                    self.dest = str(self.dest)
                    output = io2.WriteRedmonster(zpick, dest=self.dest,
                                                  clobber=True)
                except Exception as e:
                    print('Could not convert dest to string - writing to \
                            default directory and NOT clobbering old files! \
                            Exception: %r' % e)
                    output = io2.WriteRedmonster(zpick, clobber=True)

        if output:
            if len(zpick.fiberid) == 1: output.write_fiber()
            else: output.write_plate()
