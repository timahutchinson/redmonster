# Use chi2 surfaces from zfinder.py to choose and classify objects on each fiber.  Each zfind
# argument is the ENTIRE object created by Zfinder for a single template, and each zfit argument
# is the ENTIRE object created by Zfitter for a single template. Zfiti must correspond with
# zfindi.  Each 'flag' argument is the combined flag for a template, combined with
# redmonster.misc.comb_flags function.  Output is in zpick.type, a list of len()=nfibers where
# each entry is a string of objclassN.type for the best reduced chi2 amongst all input classes.
#
# N.B. this version of zpicker only works with zfitter objects created with z_refine2()
#
# Tim Hutchinson, University of Utah, July 2014
# Significantly rewritten by TH, July 2015
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
        #self.z = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        #self.z_err = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        #self.zwarning = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        self.z = []
        self.z_err = []
        self.zwarning = []
        self.npoly = []
        self.dof = specobj.dof.copy()
        self.npixstep = []
        self.models = n.zeros( (specobj.flux.shape) )
        self.nclass = len(zfindobjs)
        #self.minrchi2 = n.zeros( (zfindobjs[0].zchi2arr.shape[0],self.num_z) )
        self.minrchi2 = []
        self.rchi2diff = []
        self.zwarning = []
        try:
            self.boss_target1 = specobj.boss_target1
        except:
            try:
                self.eboss_target1 = specobj.eboss_target1
            except:
                pass

        self.classify_obj(zfindobjs, zfitobjs, flags)


    def classify_obj(self, zfindobjs, zfitobjs, flags, rchi2threshold=0.01):
        # Build dictionary of template position in object lists vs position in ordered lists (i.e. if the chi2.argmin() = 14, and tempdict[14] = 3, then the min chi2 is from the template corresponding to zfindobjs[3])
        # While looping over templates, also build out an array of rchi2 values for each
        tempdict = {}
        poslist = range(self.num_z)*len(zfindobjs) # Builds list of 0:num_z repeated ntemps times (i.e., [0,1,2,3,0,1,2,3,0,1,2,3] for num_z=4 and 3 templates)
        rchi2s = []
        for itemp in xrange(len(zfindobjs)):
            for i in xrange(self.num_z):
                tempdict[ i+(itemp*self.num_z) ] = itemp
            '''
            rchi2s.append( zfindojbs[itemp].zchi2arr.copy() ) # Add copy of chi2 array to rchi2s list
            for ifiber in xrange(zfindobjs[0].zchi2arr.shape[0]):
                rchi2s[itemp][ifiber] /= ( specobj.dof[ifiber] - zfindobjs[itemp].npoly ) # Convert from chi2 to rchi2 by dividing by (number of pixels - number of poly terms)
            '''
        for ifiber in xrange(zfindobjs[0].zchi2arr.shape[0]):
            fibermins = []
            fiberminvecs = []
            # Create tuples for all the values we're going to carry forward.  One tuple per fiber, each with num_z redshifts, classifications, etc.
            ztuple = ()
            zerrtuple = ()
            fnametuple = ()
            typetuple = ()
            minchi2tuple = ()
            vectortuple = ()
            npolytuple = ()
            npixsteptuple = ()
            for itemp in xrange(self.nclass): # Build temporary array of num_z lowest minima for each template
                for imin in xrange(self.num_z):
                    fibermins.append( zfitobjs[itemp].chi2vals[ifiber][imin] / (self.dof[ifiber] - zfindobjs[itemp].npoly) ) # Add num_z best chi2s found in zfitter divided by (number of pixels - number of poly terms) to convert to rchi2
                    fiberminvecs.append( zfitobjs[itemp].minvectors[ifiber][imin]) # Add num_z vectors for location of each chi2 in the above step
            for iz in xrange(self.num_z): # Build tuples of num_z best redshifts and classifications for this fiber
                zpos = n.asarray(fibermins).argmin() # Location of this best redshfit in fibermins array - to be fed into tempdict to find template
                tempnum = tempdict[zpos] # Location in lists of template objects of this redshift classification
                znum = poslist[zpos] # Location in zfitobj[tempnum] of this z
                '''
                self.z[ifiber][iz] = zfitobjs[tempnum].z[ifiber][znum]
                self.z_err[ifiber][iz] = zfitobjs[tempnum].z_err[ifiber][znum]
                '''
                ztuple += (zfitobjs[tempnum].z[ifiber][znum],)
                zerrtuple += (zfitobjs[tempnum].z_err[ifiber][znum],)
                fnametuple += (zfindobjs[tempnum].fname,)
                typetuple += (zfindobjs[tempnum].type,)
                #subtype
                minchi2tuple += (fibermins[zpos],)
                vectortuple += (fiberminvecs[zpos],)
                npolytuple += (zfindobjs[tempnum].npoly,)
                npixsteptuple += (zfindobjs[tempnum].npixstep,)
                if iz == 0: # Only the first flag and model are kept
                    self.models[ifiber] = zfindobjs[tempnum].models[ifiber]
                    self.zwarning.append( flags[tempnum][ifiber] )
                fibermins[zpos] = 1e9
            self.z.append(ztuple)
            self.z_err.append(zerrtuple)
            self.fname.append(fnametuple)
            self.type.append(typetuple)
            #subtype
            self.minrchi2.append(minchi2tuple)
            self.minvector.append(vectortuple)
            self.npoly.append(npolytuple)
            self.npixstep.append(npixsteptuple)
            self.rchi2diff.append( self.minrchi2[ifiber][1] - self.minrchi2[ifiber][0])
            if self.rchi2diff < rchi2threshold:
                self.flag_small_dchi2(ifiber)
        self.zwarning = map(int, self.zwarning)


        def flag_small_dchi2(self, ifiber):
            flag_val = int('0b100',2) # From BOSS zwarning flag definitions
            self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val



                

'''
    for imin in xrange(self.num_z):
    minpos = n.unravel_index(rchi2s[itemp][ifiber].argmin(), rchi2s[itemp][ifiber].shape) # Location of minimum in zfindobjs[itemp][ifiber]
    fiberminvecs.append( minpos )
    fibermins.append( rchi2s[itemp][ifiber][minpos] ) # Add this particular minimum to the temp array of minima
    rchi2s[itemp][ifiber][...,minpos[-1]-zfitobjs[itemp].width:minpos[-1]+zfitobjs[itemp].width+1] = 1e9 # Block out a window of 'width' pixels (set in zfitter) in redshift around this minimum
    '''

























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
















