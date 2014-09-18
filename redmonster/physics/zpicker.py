# Use chi2 surfaces from zfinder.py to choose and classify objects on each fiber.  Each zfind
# argument is the ENTIRE object created by Zfinder for a single template, and each zfit argument
# is the ENTIRE object created by Zfitter for a single template. Zfiti must correspond with
# zfindi.  Each 'flag' argument is the combined flag for a template, combined with
# redmonster.misc.comb_flags function.  Output is in zpick.type, a list of len()=nfibers where
# each entry is a string of objclassN.type for the best reduced chi2 amongst all input classes.
#
# Tim Hutchinson, University of Utah, July 2014
# t.hutchinson@utah.edu

import numpy as n

class Zpicker:

    def __init__(self, specobj, zfind1, zfit1, flags1, zfind2=None, zfit2=None, flags2=None, zfind3=None, zfit3=None, flags3=None, zfind4=None, zfit4=None, flags4=None, zfind5=None, zfit5=None, flags5=None):
        self.npixflux = specobj.npix
        if specobj.plate: self.plate = specobj.plate
        if specobj.mjd: self.mjd = specobj.mjd
        if specobj.fiberid: self.fiberid = specobj.fiberid
        if specobj.hdr: self.hdr = specobj.hdr
        self.fname = []
        self.type = []
        self.subtype = []
        self.minvector = []
        self.z = n.zeros( (zfind1.zchi2arr.shape[0],2) )
        self.z_err = n.zeros( (zfind1.zchi2arr.shape[0],2) )
        self.zwarning = []
        self.npoly = []
        self.dof = specobj.dof.copy()
        if not zfind2: self.nclass = 1
        elif zfind2 and not zfind3: nclass = 2
        elif zfind3 and not zfind4: nclass = 3
        elif zfind4 and not zfind5: nclass = 4
        else: nclass = 5
        self.nclass = nclass
        self.minrchi2 = n.zeros( (zfind1.zchi2arr.shape[0],nclass) )
        self.classify_obj(zfind1, zfit1, flags1, zfind2, zfit2, flags2, zfind3, zfit3, flags3, zfind4, zfit4, flags4, zfind5, zfit5, flags5)
    
    def classify_obj(self, zfind1, zfit1, flags1, zfind2, zfit2, flags2, zfind3, zfit3, flags3, zfind4, zfit4, flags4, zfind5, zfit5, flags5):
        flag_val = int('0b100',2) # From BOSS zwarning flag definitions
        for ifiber in xrange(zfind1.zchi2arr.shape[0]):
            #self.minrchi2[ifiber,0] = n.min(zfind1.zchi2arr[ifiber]) / (self.npixflux - zfind1.npoly) # Calculate reduced chi**2 values to compare templates of differing lengths
            #if zfind2: self.minrchi2[ifiber,1] = n.min(zfind2.zchi2arr[ifiber]) / (self.npixflux - zfind2.npoly)
            #if zfind3: self.minrchi2[ifiber,2] = n.min(zfind3.zchi2arr[ifiber]) / (self.npixflux - zfind3.npoly)
            #if zfind4: self.minrchi2[ifiber,3] = n.min(zfind4.zchi2arr[ifiber]) / (self.npixflux - zfind4.npoly)
            #if zfind5: self.minrchi2[ifiber,4] = n.min(zfind5.zchi2arr[ifiber]) / (self.npixflux - zfind5.npoly)
            # Try using specobj.dof instead, because it accounts for masked pixels
            self.minrchi2[ifiber,0] = n.min(zfind1.zchi2arr[ifiber]) / (self.dof[ifiber] - zfind1.npoly) # Calculate reduced chi**2 values to compare templates of differing lengths
            if zfind2: self.minrchi2[ifiber,1] = n.min(zfind2.zchi2arr[ifiber]) / (self.dof[ifiber] - zfind2.npoly)
            if zfind3: self.minrchi2[ifiber,2] = n.min(zfind3.zchi2arr[ifiber]) / (self.dof[ifiber] - zfind3.npoly)
            if zfind4: self.minrchi2[ifiber,3] = n.min(zfind4.zchi2arr[ifiber]) / (self.dof[ifiber] - zfind4.npoly)
            if zfind5: self.minrchi2[ifiber,4] = n.min(zfind5.zchi2arr[ifiber]) / (self.dof[ifiber] - zfind5.npoly)

            minpos = self.minrchi2[ifiber].argmin() # Location of best chi2 of array of best (individual template) chi2s
            
            if minpos == 0: # Means overall chi2 minimum came from template 1
                self.type.append(zfind1.type)
                self.minvector.append(zfit1.minvector[ifiber])
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
                self.minvector.append(zfit2.minvector[ifiber])
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
                self.minvector.append(zfit3.minvector[ifiber])
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
                self.minvector.append(zfit4.minvector[ifiber])
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
                self.minvector.append(zfit5.minvector[ifiber])
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
                argsort = self.minrchi2[ifiber].argsort()
                if argsort[1] == 0:
                    if ( n.min(zfind1.zchi2arr[ifiber]) - n.min(zfind5.zchi2arr[ifiber]) ) < zfit5.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 1:
                    if ( n.min(zfind2.zchi2arr[ifiber]) - n.min(zfind5.zchi2arr[ifiber]) ) < zfit5.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 2:
                    if ( n.min(zfind3.zchi2arr[ifiber]) - n.min(zfind5.zchi2arr[ifiber]) ) < zfit5.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val
                elif argsort[1] == 3:
                    if ( n.min(zfind4.zchi2arr[ifiber]) - n.min(zfind5.zchi2arr[ifiber]) ) < zfit5.threshold: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val












