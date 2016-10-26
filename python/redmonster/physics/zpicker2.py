# Use chi2 surfaces from zfinder.py to choose and classify objects on each fiber.  Each zfind
# argument is the ENTIRE object created by ZFinder for a single template, and each zfit argument
# is the ENTIRE object created by ZFitter for a single template. Zfiti must correspond with
# zfindi.  Each 'flag' argument is the combined flag for a template, combined with
# redmonster.misc.comb_flags function.  Output is in zpick.type, a list of len()=nfibers where
# each entry is a string of objclassN.type for the best reduced chi2 amongst all input classes.
#
# N.B. this version of zpicker only works with zfitter objects created with z_refine2()
#
# Tim Hutchinson, University of Utah, July 2014
# Significantly rewritten by TH, July 2015
# t.hutchinson@utah.edu

from os import environ
from os.path import join

import numpy as n
from scipy import linalg
from scipy.optimize import nnls
from astropy.io import fits

from redmonster.physics.misc import poly_array
from redmonster.datamgr.io import read_ndArch

class ZPicker:

    def __init__(self, specobj, zfindobjs, zfitobjs, flags, num_z=5):
        self.num_z = num_z # Number of redshifts to retain
        self.npixflux = specobj.npix if hasattr(specobj,'npix') else \
                specobj.flux.shape[-1]
        self.flux = specobj.flux
        self.ivar = specobj.ivar
        if hasattr(specobj,'plate'): self.plate = specobj.plate
        if hasattr(specobj,'mjd'): self.mjd = specobj.mjd
        if hasattr(specobj,'fiberid'): self.fiberid = specobj.fiberid
        if hasattr(specobj,'hdr'): self.hdr = specobj.hdr
        else: self.hdr = fits.Header()
        self.threshold = zfitobjs[0].threshold
        self.fname = []
        self.group = []
        self.type = []
        self.subtype = []
        self.minvector = []
        #self.z = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        #self.z_err = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        #self.zwarning = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        self.z = []
        self.z_err = []
        self.zwarning = []
        self.npoly = []
        self.dof = specobj.dof.copy()
        self.npixstep = []
        self.models = n.zeros( (specobj.flux.shape[0],num_z,self.npixflux) )
        self.fs = []
        self.nclass = len(zfindobjs)
        #self.minrchi2 = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        self.minrchi2 = []
        self.rchi2diff = []
        self.chi2_null = []
        self.sn2_data = []
        if hasattr(specobj,'boss_target1'):
            self.boss_target1 = specobj.boss_target1
        if hasattr(specobj,'eboss_target0'):
            self.eboss_target0 = specobj.eboss_target0
        if hasattr(specobj, 'eboss_target1'):
            self.eboss_target1 = specobj.eboss_target1
        self.classify_obj(zfindobjs, zfitobjs, flags)

    def classify_obj(self, zfindobjs, zfitobjs, flags, rchi2threshold=0.005):
        # Build dictionary of template position in object lists vs
        # position in ordered lists (i.e. if the chi2.argmin() = 14, and
        # tempdict[14] = 3, then the min chi2 is from the template
        # corresponding to zfindobjs[3])
        # While looping over templates, also build out an array of rchi2
        # values for each
        self.rchi2threshold = rchi2threshold
        tempdict = {}
        # Builds list of 0:num_z repeated ntemps times (i.e.,
        # [0,1,2,3,0,1,2,3,0,1,2,3] for num_z=4 and 3 templates)
        poslist = list(range(self.num_z))*len(zfindobjs)
        rchi2s = []
        for itemp in range(len(zfindobjs)):
            for i in range(self.num_z):
                tempdict[ i+(itemp*self.num_z) ] = itemp
            '''
            rchi2s.append( zfindojbs[itemp].zchi2arr.copy() ) # Add copy
            of chi2 array to rchi2s list
            for ifiber in xrange(zfindobjs[0].zchi2arr.shape[0]):
                rchi2s[itemp][ifiber] /= ( specobj.dof[ifiber] -
                zfindobjs[itemp].npoly ) # Convert from chi2 to rchi2 by
                dividing by (number of pixels - number of poly terms)
            '''
        for ifiber in range(zfindobjs[0].zchi2arr.shape[0]):
            fibermins = []
            fiberminvecs = []
            # Create tuples for all the values we're going to carry forward.
            # One tuple per fiber, each with num_z redshifts,
            # classifications, etc.
            ztuple = ()
            zerrtuple = ()
            fnametuple = ()
            grouptuple = ()
            typetuple = ()
            subtypetuple = ()
            minchi2tuple = ()
            vectortuple = ()
            npolytuple = ()
            npixsteptuple = ()
            fstuple = ()
            # Catch spectra that are all 0's and return null result
            if n.all(self.flux[ifiber] == 0.0) or n.all(self.ivar[ifiber] == 0.0):
                ztuple = (-1,)*self.num_z
                zerrtuple = (-1,)*self.num_z
                fnametuple = ('noSpectrum',)*self.num_z
                grouptuple = (-1,)*self.num_z
                typetuple = ('noSpectrum',)*self.num_z
                subtypetuple = ('noSpectrum',)*self.num_z
                minchi2tuple = (-1,)*self.num_z
                vectortuple = (-1,)*self.num_z
                npolytuple = (-1,)*self.num_z
                npixsteptuple = (-1,)*self.num_z
                fstuple = (-1,)*self.num_z
                self.zwarning.append(int(flags[0][0]))
                self.flag_small_dchi2(ifiber)
            else:
                # Build temporary array of num_z lowest minima for each template
                for itemp in range(self.nclass):
                    for imin in range(self.num_z):
                        try:
                            # Add num_z best chi2s found in zfitter divided by
                            # (number of pixels - number of poly terms) to
                            # convert to rchi2
                            fibermins.append(
                                    zfitobjs[itemp].chi2vals[ifiber][imin] /
                                             (self.dof[ifiber] -
                                              zfindobjs[itemp].npoly) )
                            # Add num_z vectors for location of each chi2 in
                            # the above step
                            fiberminvecs.append(
                                    zfitobjs[itemp].minvectors[ifiber][imin])
                        except IndexError as e:
                            print(repr(e))
                            fiberminvecs.append( (-1,) )
                            if len(zfitobjs[itemp].chi2vals[ifiber]) > 0:
                                fibermins.append( \
                                       n.max(zfitobjs[itemp].chi2vals[ifiber]) / \
                                       (self.dof[ifiber] - zfindobjs[itemp].npoly))
                            else:
                                fibermins.append(100000.)

                            # this is ok, might be no solution
                            print("WARNING no z fitted for fiber #%d, class #%d, zid #%d"%(ifiber,itemp,imin))

                # Build tuples of num_z best redshifts and classifications
                # for this fiber
                #for iz in xrange(self.num_z):
                iz = 0
                while iz < self.num_z:
                    # Location of this best redshfit in fibermins array - to
                    # be fed into tempdict to find template
                    zpos = n.asarray(fibermins).argmin()
                    # Location in lists of template objects of this redshift
                    # classification
                    tempnum = tempdict[zpos]
                    # Location in zfitobj[tempnum] of this z
                    znum = poslist[zpos]
                    # Check for repeat z
                    if fiberminvecs[zpos] != (-1,):
                        ztuple += (zfitobjs[tempnum].z[ifiber][znum],)
                        zerrtuple += (zfitobjs[tempnum].z_err[ifiber][znum],)
                        fnametuple += (zfindobjs[tempnum].fname,)
                        grouptuple += (zfindobjs[tempnum].group,)
                        typetuple += (zfindobjs[tempnum].type,)
                        vectortuple += (fiberminvecs[zpos],)
                        d = {} # Dictionary for subtype
                        for j in range( len(vectortuple[-1][:-1]) ):
                            d[zfindobjs[tempnum].infodict['par_names'][j]] = \
                                    zfindobjs[tempnum].baselines[j] \
                                            [vectortuple[-1][j]]
                        subtypetuple += (d,)
                        minchi2tuple += (fibermins[zpos],)
                        npolytuple += (zfindobjs[tempnum].npoly,)
                        npixsteptuple += (zfindobjs[tempnum].npixstep,)
                        if iz == 0: # Only the first flag is kept
                            self.zwarning.append( flags[tempnum][ifiber] )
                            self.chi2_null.append(
                                    zfindobjs[tempnum].chi2_null[ifiber])
                            self.sn2_data.append(
                                    zfindobjs[tempnum].sn2_data[ifiber])
                        fibermins[zpos] = 1e9
                        self.models[ifiber,iz], f = self.create_model(
                                fnametuple[iz], npolytuple[iz],
                                npixsteptuple[iz], vectortuple[iz],
                                zfindobjs[tempnum], self.flux[ifiber],
                                self.ivar[ifiber])
                        fstuple += (f,)
                        iz += 1
                    else:
                        fibermins[zpos] = 1e9

            self.z.append(ztuple)
            self.z_err.append(zerrtuple)
            self.fname.append(fnametuple)
            self.group.append(grouptuple)
            self.type.append(typetuple)
            self.subtype.append(subtypetuple)
            self.minrchi2.append(minchi2tuple)
            self.minvector.append(vectortuple)
            self.npoly.append(npolytuple)
            self.npixstep.append(npixsteptuple)
            self.rchi2diff.append( self.minrchi2[ifiber][1] - \
                                  self.minrchi2[ifiber][0])
            self.fs.append( fstuple )
            c_kms = 299792.458
            if self.rchi2diff[ifiber] < self.rchi2threshold:
                #if (n.abs(self.z[ifiber][0] - self.z[ifiber][1]) /
                    #n.sqrt(self.z_err[ifiber][0]**2 +
                            #self.z_err[ifiber][1]**2)) > 1:
                if self.z[ifiber][0] != -1:
                    if (c_kms*n.abs(self.z[ifiber][0] - self.z[ifiber][1])) / \
                            (1 + self.z[ifiber][0]) > 1000:
                        if not any(a in self.group[ifiber][0] for
                                   a in self.group[ifiber][1]):
                            self.flag_small_dchi2(ifiber)
                else:
                    self.flag_small_dchi2(ifiber)
            if n.isnan(self.rchi2diff[ifiber]):
                self.flag_small_dchi2(ifiber)
            self.flag_null_fit(ifiber, flags)
        self.zwarning = list(map(int, self.zwarning))

    def flag_small_dchi2(self, ifiber):
        """Set the small delta chi**2 zwarning flag."""
        flag_val = int('0b100',2) # From BOSS zwarning flag definitions
        self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val

    def flag_null_fit(self, ifiber, flags):
        """Set flag if any template classes had a null fit."""
        null_fit_flag = int('0b100000000',2)
        for template in flags:
            if int(template[ifiber]) & null_fit_flag > 0:
                self.zwarning[ifiber] = int(self.zwarning[ifiber]) | \
                        null_fit_flag

    def create_model(self, fname, npoly, npixstep, minvector, zfindobj,
                     flux, ivar):
        """Return the best fit model for a given template at a
            given redshift.
        """
        try:
            pixoffset = zfindobj.pixoffset
            temps = read_ndArch( join( environ['REDMONSTER_TEMPLATES_DIR'],
                                      fname ) )[0]
            pmat = n.zeros( (self.npixflux, npoly+1) )
            this_temp = temps[minvector[:-1]]
            pmat[:,0] = this_temp[(minvector[-1]*npixstep)+pixoffset:\
                                  (minvector[-1]*npixstep)+pixoffset + \
                                  self.npixflux]
            polyarr = poly_array(npoly, self.npixflux)
            pmat[:,1:] = n.transpose(polyarr)
            ninv = n.diag(ivar)
            f = linalg.solve( n.dot(n.dot(n.transpose(pmat),ninv),pmat),
                               n.dot( n.dot(n.transpose(pmat),ninv),flux) ); \
                    f = n.array(f)
            if f[0] < 0:
                try:
                    f = nnls( n.dot(n.dot(n.transpose(pmat),ninv),pmat),
                             n.dot( n.dot(n.transpose(pmat),ninv),flux) )[0]; \
                            f = n.array(f)
                    return n.dot(pmat,f), tuple(f)
                except Exception as e:
                    print("Exception: %r" % e)
                    return n.zeros(self.npixflux), (0,)
            else:
                return n.dot(pmat,f), tuple(f)
        except Exception as e:
            print("Exception: %r" % e)
            return n.zeros(self.npixflux), (0,)
