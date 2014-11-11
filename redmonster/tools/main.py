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
#                $BOSS_SPECTRO_REDUX/$RUN2D/pppp/$RUN1D/ , where pppp is the 4 digit plate id.
# clobber (Boolean): Default behavior is to overwrite older output files for same plate/mjd.  Setting to false will cause new
#                    version to be written.
#
# Tim Hutchinson, University of Utah, October 2014
# t.hutchinson@utah.edu

import numpy as n
from astropy.io import fits
from redmonster.datamgr import spec, io
from redmonster.physics import zfinder, zfitter, zpicker
from redmonster.sandbox import yanny as y
from redmonster.physics import misc
from time import gmtime, strftime
import matplotlib.pyplot as p
p.interactive(True)
import sys

def main(plate, mjd, templates, fiberid=None, zmin=None, zmax=None, npoly=None, npixstep=None, dest=None, clobber=True):
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
    if type(templates) is str:
        try: templates = [templates]
        except:
            print 'Templates not a list and unable to convert to list!'
            sys.exit(1)
    if type(templates) is list: templates = map(str, templates)
    if zmin is not None:
        if type(zmin) is not list:
            try:
                zmin = [zmin]
            except:
                try:
                    zmin = zmin.tolist()
                except:
                    print 'Can\'t convert zmin to list - defaulting to full zrange!'
                    zmin = None
                    zmax = None
        if type(zmin) is list:
            if len(zmin) != len(templates):
                print 'Length of zmin doesn\'t match length of templates - defaulting to full zrange!'
                zmin = None
                zmax = None
            if zmax is None:
                print 'zmax not given - defaulting to full zrange!'
                zmin = None
                zmax = None
            else:
                if type(zmax) is not list:
                    try: zmax = [zmax]
                    except:
                        try:
                            zmax = zmax.tolist()
                        except:
                            print 'Can\'t convert zmax to list - defaulting to full zrange!'
                            zmin = None
                            zmax = None
                if len(zmin) != len(zmax):
                    print 'Length of zmin and zmax don\'t match - defaulting to full zrange!'
                    zmin = None
                    zmax = None
#import pdb; pdb.set_trace()
    if npoly is None:
        npoly = [4]*len(templates)
    else:
        if type(npoly) is not list:
            try:
                npoly = [npoly]
            except:
                try:
                    npoly = npoly.tolist()
                except:
                    print 'npoly not a list and unable to convert to list - defaulting to npoly=4 for all templates!'
                    npoly = [4]*len(templates)
        else:
            npoly = map(int, npoly)
    if npixstep is None:
        npixstep = [1]*len(templates)
    else:
        if type(npixstep) is not list:
            try:
                npixstep = [npixstep]
            except:
                try:
                    npixstep = npixstep.tolist()
                except:
                    print 'npixstep not a list and unable to convert to list - defaulting to npixstep=1 for all templates!'
                    npixstep = [1]*len(templates)
        else:
            npixstep = map(int, npixstep)

    # Spec
    specs = spec.Spec(plate=plate, mjd=mjd, fiberid=fiberid)

    # Zfinder, Zfitter
    zfindobjs = []
    zfitobjs = []
    if (zmin is not None) & (zmax is not None):
        for i in xrange(len(templates)):
            zfindobjs.append( zfinder.Zfinder(fname=templates[i], npoly=npoly[i], zmin=zmin[i], zmax=zmax[i]) )
            zfindobjs[i].zchi2( specs.flux, specs.loglambda, specs.ivar, npixstep=npixstep[i] )
            zfitobjs.append( zfitter.Zfitter(zfindobjs[i].zchi2arr, zfindobjs[i].zbase) )
            zfitobjs[i].z_refine()
    else:
        for i in xrange(len(templates)):
            zfindobjs.append( zfinder.Zfinder(fname=templates[i], npoly=npoly[i], npixstep=npixstep[i]) )
            zfindobjs[i].zchi2( specs.flux, specs.loglambda, specs.ivar, npixstep=npixstep[i] )
            zfitobjs.append( zfitter.Zfitter(zfindobjs[i].zchi2arr, zfindobjs[i].zbase) )
            zfitobjs[i].z_refine()

    # Flags
    flags = []
    for i in xrange(len(zfindobjs)):
        flags.append( misc.comb_flags(specs, zfindobjs[i], zfitobjs[i]) )

    # Zpicker
    if len(templates) == 1: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0])
    elif len(templates) == 2: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0], zfindobjs[1], zfitobjs[1], flags[1])
    elif len(templates) == 3: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0], zfindobjs[1], zfitobjs[1], flags[1],
                                                      zfindobjs[2], zfitobjs[2], flags[2])
    elif len(templates) == 4: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0], zfindobjs[1], zfitobjs[1], flags[1],
                                                      zfindobjs[2], zfitobjs[2], flags[2], zfindobjs[3], zfitobjs[3], flags[3])
    elif len(templates) == 5: zpick = zpicker.Zpicker(specs, zfindobjs[0], zfitobjs[0], flags[0], zfindobjs[1], zfitobjs[1], flags[1],
                                                      zfindobjs[2], zfitobjs[2], flags[2], zfindobjs[3], zfitobjs[3], flags[3],
                                                      zfindobjs[4], zfitobjs[4], flags[4])

    # Write output
    if dest is None:
        if clobber is False:
            output = io.Write_Redmonster(zpick, clobber=False)
        else:
            output = io.Write_Redmonster(zpick)
    else:
        if type(dest) is str:
            if clobber is False:
                output = io.Write_Redmonster(zpick, dest=dest, clobber=False)
            else:
                output = io.Write_Redmonster(zpick, dest=dest)
        else:
            try:
                dest = str(dest)
                if clobber is False:
                    output = io.Write_Redmonster(zpick, dest=dest, clobber=False)
                else:
                    output = io.Write_Redmonster(zpick, dest=dest)
            except:
                print 'Could not convert dest to string - writing to default directory and NOT clobbering old files!'
                output = io.Write_Redmonster(zpick, clobber=False)

















