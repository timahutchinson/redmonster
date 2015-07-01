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

    def __init__(self, specobj, zfindobjs, zfitobjs, flags, num_z=5):
        self.num_z = num_z # Number of redshifts to retain
        self.npixflux = specobj.npix
        if specobj.plate: self.plate = specobj.plate
        if specobj.mjd: self.mjd = specobj.mjd
        if specobj.fiberid: self.fiberid = specobj.fiberid
        if specobj.hdr: self.hdr = specobj.hdr
        self.fname = []
        self.type = []
        self.subtype = []
        self.minvector = []
        self.z = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        self.z_err = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        #self.zwarning = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        #self.z = []
        #self.z_err = []
        self.zwarning = []
        self.npoly = []
        self.dof = specobj.dof.copy()
        self.npixstep = []
        self.models = n.zeros( (specobj.flux.shape) )
        self.nclass = len(zfindobjs)
        self.minrchi2 = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.nclass*self.num_z) )
        self.chi2diff = []
        try:
            self.boss_target1 = specobj.boss_target1
        except:
            try:
                self.eboss_target1 = specobj.eboss_target1
            except:
                pass

        self.classify_obj(zfindobjs, zfitobjs, flags)


    def classify_obj(zfindobjs, zfitobjs, flags):
        for ifiber in xrange(zfindobjbs[0].zchi2arr.shape[0]):
            # Convert to rchi2 for comparing templates of differing lengths
            for itemp in xrange(len(zfindobjs)): # Loop over template classes
                for iz in xrange(self.num_z): # Build out num_z options for each template class, in case the overall best num_z are from the same
                    minloc = n.unravel_index(zfindobjs[itemp].zchi2arr[ifiber].argmin(),zfindobjs[itemp].zchi2arr[ifiber].shape)
                    self.minrchi2[ifiber, iz+(itemp*self.num_z)] = n.min(zfindobjs[itemp].zchi2arr[ifiber]) / (self.dof[ifiber]-zfindobjs[itemp].npoly)
    



















    



'''

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

'''
















