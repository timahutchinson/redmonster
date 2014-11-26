# Use chi2 surfaces from zfinder.py to choose and classify objects on each fiber.  Each zfind
# argument is the ENTIRE object created by Zfinder for a single template, and each zfit argument
# is the ENTIRE object created by Zfitter for a single template. Zfiti must correspond with
# zfindi.  Each 'flag' argument is the combined flag for a template, combined with
# redmonster.misc.comb_flags function.  Output is in zpick.type, a list of len()=nfibers where
# each entry is a string of objclassN.type for the best reduced chi2 amongst all input classes.
#
# Tim Hutchinson, University of Utah, July 2014
# Significantly rewritten by TH, November 2014
# t.hutchinson@utah.edu

import numpy as n

class Zpicker:

    def __init__(self, specobj, zfindobjs, zfitobjs, flags):
        self.npixflux = specobj.npix
        if specobj.plate: self.plate = specobj.plate
        if specobj.mjd: self.mjd = specobj.mjd
        if specobj.fiberid: self.fiberid = specobj.fiberid
        if specobj.hdr: self.hdr = specobj.hdr
        self.fname = []
        self.type = []
        self.subtype = []
        self.minvector = []
        self.z = n.zeros( (zfindobjs[0].zchi2arr.shape[0],5) )
        self.z_err = n.zeros( (zfindobjs[0].zchi2arr.shape[0],5) )
        self.zwarning = n.zeros( (zfindobjs[0].zchi2arr.shape[0],5) )
        self.npoly = []
        self.dof = specobj.dof.copy()
        self.npixstep = []
        self.models = n.zeros( (specobj.flux.shape) )
        self.nclass = len(zfindobjs)
        self.minrchi2 = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.nclass*5) )
        self.classify_obj(zfindobjs, zfitobjs, flags)
    
    def classify_obj(self, zfindobjs, zfitobjs, flags):
        flag_val = int('0b100',2) # From BOSS zwarning flag definitions for small dchi2
        objtypes = []
        for ifiber in xrange(zfindobjs[0].zchi2arr.shape[0]):
            for itemp in xrange(self.nclass):
                for i in xrange(5):
                    import pdb; pdb.set_trace()
                    self.minrchi2[ifiber,(itemp*5)+i] = n.min(zfindobjs[itemp].zchi2arr[ifiber]) / (self.dof[ifiber] - zfindobjs[itemp].npoly)
                    loc = n.unravel_index(zfindobjs[itemp].zchi2arr[ifiber].argmin(), zfindobjs[itemp].zchi2arr[ifiber].shape)
                    zfindobjs[itemp].zchi2arr[ifiber][loc[:-1]][(loc[-1]-zfitobjs[itemp].width):(loc[-1]+zfitobjs[itemp].width)] = 10000000.
                for i in xrange(5):
                    objtypes.append(zfindobjs[itemp].type)

            
            minpos = []
            for ipos in xrange(5):
                minpos.append(self.minrchi2[ifiber].argmin())
                self.minrchi2[ifiber,minpos[ipos]] = 10000000.
                while (int(flags[minpos[ipos]/5][ifiber]) & 8) == 8:
                    minpos[ipos] = self.minrchi2[ifiber].argmin()
                    self.minrchi2[ifiber,minpos[ipos]] = 10000000.
                self.type.append(objtypes[ipos])
            
#            minpos1 = self.minrchi2[ifiber].argmin() # Location of best chi2 of array of best (individual template) chi2s
#            self.minrchi2[ifiber,minpos1] = 10000000.
#            while (flags[minpos1/5][ifiber] & 8) == 8:
#                minpos1 = self.minrchi2[ifiber].argmin()
#                minpos1 = self.minrchi2[ifiber,minpos1] = 10000000.
#
#            minpos2 = self.minrchi2[ifiber].argmin()
#            self.minrchi2[ifiber,minpos2] = 10000000.
#            while (flags[minpos2/5][ifiber] & 8) == 8:
#                minpos2 = self.minrchi2[ifiber].argmin()
#                minpos2 = self.minrchi2[ifiber,minpos2] = 10000000.
#
#            minpos3 = self.minrchi2[ifiber].argmin()
#            self.minrchi2[ifiber,minpos3] = 10000000.
#            while (flags[minpos3/5][ifiber] & 8) == 8:
#                minpos3 = self.minrchi2[ifiber].argmin()
#                minpos3 = self.minrchi2[ifiber,minpos3] = 10000000.
#
#            minpos4 = self.minrchi2[ifiber].argmin()
#            self.minrchi2[ifiber,minpos4] = 10000000.
#            while (flags[minpos4/5][ifiber] & 8) == 8:
#                minpos4 = self.minrchi2[ifiber].argmin()
#                minpos4 = self.minrchi2[ifiber,minpos4] = 10000000.
#
#            minpos5 = self.minrchi2[ifiber].argmin()
#            self.minrchi2[ifiber,minpos5] = 10000000.
#            while (flags[minpos5/5][ifiber] & 8) == 8:
#                minpos5 = self.minrchi2[ifiber].argmin()
#                minpos5 = self.minrchi2[ifiber,minpos5] = 10000000.













            if minpos == 0: # Means overall chi2 minimum came from template 1
                self.type.append(zfind1.type)
                self.minvector.append(zfit1.minvector[ifiber][1:])
                self.z[ifiber] = zfit1.z[ifiber]
                self.z_err[ifiber] = zfit1.z_err[ifiber]
                self.fname.append(zfind1.fname)
                minloc = n.unravel_index(zfind1.zchi2arr[ifiber].argmin(), zfind1.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind1.infodict['par_names'][i]] = zfind1.baselines[i][minloc[i]]
                self.subtype.append(d)
                self.zwarning = n.append(self.zwarning, flags1[ifiber])
                self.dof[ifiber] -= zfind1.npoly
                self.npoly.append(zfind1.npoly)
                self.npixstep.append(zfind1.npixstep)
                self.models[ifiber] = zfind1.models[ifiber]
                argsort = self.minrchi2[ifiber].argsort()
                if len(argsort) > 1:
                    if argsort[1] == 1:
                        if ( n.min(zfind2.zchi2arr[ifiber]) - n.min(zfind1.zchi2arr[ifiber]) ) < zfit1.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val #THIS THRESHOLD PROBABLY ISN'T RIGHT AND NEEDS TO BE CHANGED
                    elif argsort[1] == 2:
                        if ( n.min(zfind3.zchi2arr[ifiber]) - n.min(zfind1.zchi2arr[ifiber]) ) < zfit1.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                    elif argsort[1] == 3:
                        if ( n.min(zfind4.zchi2arr[ifiber]) - n.min(zfind1.zchi2arr[ifiber]) ) < zfit1.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                    elif argsort[1] == 4:
                        if ( n.min(zfind5.zchi2arr[ifiber]) - n.min(zfind1.zchi2arr[ifiber]) ) < zfit1.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
        
            elif minpos == 1: # Means overall chi2 minimum came from template 2
                self.type.append(zfind2.type)
                self.minvector.append(zfit2.minvector[ifiber][1:])
                self.z[ifiber] = zfit2.z[ifiber]
                self.z_err[ifiber] = zfit2.z_err[ifiber]
                self.fname.append(zfind2.fname)
                minloc = n.unravel_index(zfind2.zchi2arr[ifiber].argmin(), zfind2.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind2.infodict['par_names'][i]] = zfind2.baselines[i][minloc[i]]
                self.subtype.append(d)
                self.zwarning = n.append(self.zwarning, flags2[ifiber])
                self.dof[ifiber] = self.dof[ifiber] - zfind2.npoly
                self.npoly.append(zfind2.npoly)
                self.npixstep.append(zfind2.npixstep)
                self.models[ifiber] = zfind2.models[ifiber]
                argsort = self.minrchi2[ifiber].argsort()
                if argsort[1] == 0:
                    if ( n.min(zfind1.zchi2arr[ifiber]) - n.min(zfind2.zchi2arr[ifiber]) ) < zfit2.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 2:
                    if ( n.min(zfind3.zchi2arr[ifiber]) - n.min(zfind2.zchi2arr[ifiber]) ) < zfit2.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 3:
                    if ( n.min(zfind4.zchi2arr[ifiber]) - n.min(zfind2.zchi2arr[ifiber]) ) < zfit2.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 4:
                    if ( n.min(zfind5.zchi2arr[ifiber]) - n.min(zfind2.zchi2arr[ifiber]) ) < zfit2.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val

            elif minpos == 2: # Means overall chi2 minimum came from template 3
                self.type.append(zfind3.type)
                self.minvector.append(zfit3.minvector[ifiber][1:])
                self.z[ifiber] = zfit3.z[ifiber]
                self.z_err[ifiber] = zfit3.z_err[ifiber]
                self.fname.append(zfind3.fname)
                minloc = n.unravel_index(zfind3.zchi2arr[ifiber].argmin(), zfind3.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind3.infodict['par_names'][i]] = zfind3.baselines[i][minloc[i]]
                self.subtype.append(d)
                self.zwarning = n.append(self.zwarning, flags3[ifiber])
                self.dof[ifiber] = self.dof[ifiber] - zfind3.npoly
                self.npoly.append(zfind3.npoly)
                self.npixstep.append(zfind3.npixstep)
                self.models[ifiber] = zfind3.models[ifiber]
                argsort = self.minrchi2[ifiber].argsort()
                if argsort[1] == 0:
                    if ( n.min(zfind1.zchi2arr[ifiber]) - n.min(zfind3.zchi2arr[ifiber]) ) < zfit3.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 1:
                    if ( n.min(zfind2.zchi2arr[ifiber]) - n.min(zfind3.zchi2arr[ifiber]) ) < zfit3.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 3:
                    if ( n.min(zfind4.zchi2arr[ifiber]) - n.min(zfind3.zchi2arr[ifiber]) ) < zfit3.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 4:
                    if ( n.min(zfind5.zchi2arr[ifiber]) - n.min(zfind3.zchi2arr[ifiber]) ) < zfit3.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val

            elif minpos == 3: # Means overall chi2 minimum came from template 4
                self.type.append(zfind4.type)
                self.minvector.append(zfit4.minvector[ifiber][1:])
                self.z[ifiber] = zfit4.z[ifiber]
                self.z_err[ifiber] = zfit4.z_err[ifiber]
                self.fname.append(zfind4.fname)
                minloc = n.unravel_index(zfind4.zchi2arr[ifiber].argmin(), zfind4.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind4.infodict['par_names'][i]] = zfind4.baselines[i][minloc[i]]
                self.subtype.append(d)
                self.zwarning = n.append(self.zwarning, flags4[ifiber])
                self.dof[ifiber] = self.dof[ifiber] - zfind4.npoly
                self.npoly.append(zfind4.npoly)
                self.npixstep.append(zfind4.npixstep)
                self.models[ifiber] = zfind4.models[ifiber]
                argsort = self.minrchi2[ifiber].argsort()
                if argsort[1] == 0:
                    if ( n.min(zfind1.zchi2arr[ifiber]) - n.min(zfind4.zchi2arr[ifiber]) ) < zfit4.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 1:
                    if ( n.min(zfind2.zchi2arr[ifiber]) - n.min(zfind4.zchi2arr[ifiber]) ) < zfit4.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 2:
                    if ( n.min(zfind3.zchi2arr[ifiber]) - n.min(zfind4.zchi2arr[ifiber]) ) < zfit4.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 4:
                    if ( n.min(zfind5.zchi2arr[ifiber]) - n.min(zfind4.zchi2arr[ifiber]) ) < zfit4.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val


            elif minpos == 4: # Means overall chi2 minimum came from template 5
                self.type.append(zfind5.type)
                self.minvector.append(zfit5.minvector[ifiber][1:])
                self.z[ifiber] = zfit5.z[ifiber]
                self.z_err[ifiber] = zfit5.z_err[ifiber]
                self.fname.append(zfind5.fname)
                minloc = n.unravel_index(zfind5.zchi2arr[ifiber].argmin(), zfind5.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind5.infodict['par_names'][i]] = zfind5.baselines[i][minloc[i]]
                self.subtype.append(d)
                self.zwarning = n.append(self.zwarning, flags5[ifiber])
                self.dof[ifiber] = self.dof[ifiber] - zfind5.npoly
                self.npoly.append(zfind5.npoly)
                self.npixstep.append(zfind5.npixstep)
                self.models[ifiber] = zfind5.models[ifiber]
                argsort = self.minrchi2[ifiber].argsort()
                if argsort[1] == 0:
                    if ( n.min(zfind1.zchi2arr[ifiber]) - n.min(zfind5.zchi2arr[ifiber]) ) < zfit5.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 1:
                    if ( n.min(zfind2.zchi2arr[ifiber]) - n.min(zfind5.zchi2arr[ifiber]) ) < zfit5.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 2:
                    if ( n.min(zfind3.zchi2arr[ifiber]) - n.min(zfind5.zchi2arr[ifiber]) ) < zfit5.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 3:
                    if ( n.min(zfind4.zchi2arr[ifiber]) - n.min(zfind5.zchi2arr[ifiber]) ) < zfit5.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val


















