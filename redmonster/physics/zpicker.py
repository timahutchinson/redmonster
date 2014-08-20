# Use chi2 surfaces from zfinder.py to choose and classify objects on each fiber.  Each zfind
# argument is the ENTIRE object created by Zfinder for a single template, and each zfit argument
# is the ENTIRE object created by Zfitter for a single template. Zfiti must correspond with
# zfindi.  Output is in zpick.type, a list of len()=nfibers where each entry is a string of
# objclassN.type for the best reduced chi2 amongst all input classes.
#
# Tim Hutchinson, University of Utah, July 2014
# t.hutchinson@utah.edu

import numpy as n

class Zpicker:

    def __init__(self, specobj, zfind1, zfit1, zfind2=None, zfit2=None, zfind3=None, zfit3=None, zfind4=None, zfit4=None, zfind5=None, zfit5=None):
        self.npixflux = specobj.npix
        self.plate = specobj.plate
        self.mjd = specobj.mjd
        self.fiberid = specobj.fiberid
        self.hdr = specobj.hdr
        self.type = []
        self.subtype = []
        self.minvector = []
        self.z = n.zeros( (zfind1.zchi2arr.shape[0],2) )
        self.z_err = n.zeros( (zfind1.zchi2arr.shape[0],2) )
        if not zfind2: self.nclass = 1
        elif zfind2 and not zfind3: nclass = 2
        elif zfind3 and not zfind4: nclass = 3
        elif zfind4 and not zfind5: nclass = 4
        else: nclass = 5
        self.minrchi2 = n.zeros( (zfind1.zchi2arr.shape[0],nclass) )
        self.classify_obj(zfind1, zfit1, zfind2, zfit2, zfind3, zfit3, zfind4, zfit4, zfind5, zfit5)
    
    def classify_obj(self, zfind1, zfit1, zfind2, zfit2, zfind3, zfit3, zfind4, zfit4, zfind5, zfit5):
        for ifiber in xrange(zfind1.zchi2arr.shape[0]):
            self.minrchi2[ifiber,0] = n.min(zfind1.zchi2arr[ifiber]) / (self.npixflux - zfind1.npoly)
            if zfind2: self.minrchi2[ifiber,1] = n.min(zfind2.zchi2arr[ifiber]) / (self.npixflux - zfind2.npoly)
            if zfind3: self.minrchi2[ifiber,2] = n.min(zfind3.zchi2arr[ifiber]) / (self.npixflux - zfind3.npoly)
            if zfind4: self.minrchi2[ifiber,3] = n.min(zfind4.zchi2arr[ifiber]) / (self.npixflux - zfind4.npoly)
            if zfind5: self.minrchi2[ifiber,4] = n.min(zfind5.zchi2arr[ifiber]) / (self.npixflux - zfind5.npoly)
            minpos = self.minrchi2[ifiber].argmin()
            if minpos == 0:
                self.type.append(zfind1.type)
                self.minvector.append(zfit1.minvector[ifiber])
                self.z[ifiber] = zfit1.z[ifiber]
                self.z_err[ifiber] = zfit1.z_err[ifiber]
                minloc = n.unravel_index(zfind1.zchi2arr[ifiber].argmin(), zfind1.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind1.infodict['par_names'][i]] = zfind1.baselines[i][minloc[i]]
                self.subtype.append(d)
            elif minpos == 1:
                self.type.append(zfind2.type)
                self.minvector.append(zfit2.minvector[ifiber])
                self.z[ifiber] = zfit2.z[ifiber]
                self.z_err[ifiber] = zfit2.z_err[ifiber]
                minloc = n.unravel_index(zfind2.zchi2arr[ifiber].argmin(), zfind2.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind2.infodict['par_names'][i]] = zfind2.baselines[i][minloc[i]]
                self.subtype.append(d)
            elif minpos == 2:
                self.type.append(zfind3.type)
                self.minvector.append(zfit3.minvector[ifiber])
                self.z[ifiber] = zfit3.z[ifiber]
                self.z_err[ifiber] = zfit3.z_err[ifiber]
                minloc = n.unravel_index(zfind3.zchi2arr[ifiber].argmin(), zfind3.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind3.infodict['par_names'][i]] = zfind3.baselines[i][minloc[i]]
                self.subtype.append(d)
            elif minpos == 3:
                self.type.append(zfind4.type)
                self.minvector.append(zfit4.minvector[ifiber])
                self.z[ifiber] = zfit4.z[ifiber]
                self.z_err[ifiber] = zfit4.z_err[ifiber]
                minloc = n.unravel_index(zfind4.zchi2arr[ifiber].argmin(), zfind4.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind4.infodict['par_names'][i]] = zfind4.baselines[i][minloc[i]]
                self.subtype.append(d)
            elif minpos == 4:
                self.type.append(zfind5.type)
                self.minvector.append(zfit5.minvector[ifiber])
                self.z[ifiber] = zfit5.z[ifiber]
                self.z_err[ifiber] = zfit5.z_err[ifiber]
                minloc = n.unravel_index(zfind5.zchi2arr[ifiber].argmin(), zfind5.zchi2arr[ifiber].shape)[:-1]
                d = {}
                for i in xrange(len(minloc)):
                    d[zfind5.infodict['par_names'][i]] = zfind5.baselines[i][minloc[i]]
                self.subtype.append(d)
