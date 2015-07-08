# Subgrid refinement and error estimation of redshift value found by redmonster.physics.zfinder.py .
# Interpolates both between redshift pixel lags and between model parameters.
#
# Tim Hutchinson, University of Utah @ IAC, May 2014
# Significant updates to z_refine() by TH, June 2015
# t.hutchinson@utah.edu

import numpy as n
from redmonster.physics.misc import quadfit
import matplotlib as m
from matplotlib import pyplot as p
m.interactive(True)
from redmonster.physics import grid_spline as gs

class Zfitter:

    def __init__(self, zchi2, zbase):
        self.zchi2 = zchi2
        self.zbase = zbase
        self.z = n.zeros((zchi2.shape[0],5))
        self.z_err = n.zeros((zchi2.shape[0],5))
        self.minvector = []
        self.zwarning = n.zeros(zchi2.shape[0])
    

    def z_refine(self, threshold=46.6, width=15): # Default threshold of 46.6 is delta (chi_r)**2 = .01, and width is 1000 km/s
        self.threshold = threshold
        self.width = width
        for ifiber in xrange(self.zchi2.shape[0]):
            self.minvector.append( (ifiber,) + n.unravel_index(self.zchi2[ifiber].argmin(),self.zchi2[ifiber].shape))
            bestzvec = n.zeros( self.zchi2.shape[-1])
            for iz in xrange(self.zchi2.shape[-1]):
                bestzvec[iz] = n.min( self.zchi2[ifiber,...,iz] )
            posinvec = n.where( bestzvec == n.min(bestzvec) )[0][0]
            if (posinvec == 0) or (posinvec == bestzvec.shape[0]-1): # Flag and skip interpolation fit if best chi2 is at edge of z-range
                self.flag_z_fitlimit(ifiber)
                self.z[ifiber] = -1.
                self.z_err[ifiber] = -1.
            else:
                # Find best redshift using global chi2 minimum
                xp = n.linspace(self.zbase[posinvec-1], self.zbase[posinvec+1], 1000)
                f = quadfit(self.zbase[posinvec-1:posinvec+2], bestzvec[posinvec-1:posinvec+2])
                fit = quad_for_fit(xp, f[0], f[1], f[2])
                #p.plot(xp, fit, color='red', hold=False)
                #p.plot(self.zbase[posinvec-1:posinvec+2], bestzvec[posinvec-1:posinvec+2], 'ko', hold=True)
                self.z[ifiber,0] = xp[n.where(fit == n.min(fit))[0][0]]
                self.z_err[ifiber,0] = self.estimate_z_err(xp, fit)
                # Find second- and third-, and fourth-best redshifts, just to have them
                zspline = gs.GridSpline(bestzvec)
                zminlocs = n.round(zspline.get_min())
                zminvals = zspline.get_val(zminlocs)
                if len(zminvals) == 1:
                    self.z[ifiber,1] = -1.
                    self.z_err[ifiber,1] = -1.
                    self.z[ifiber,2] = -1
                    self.z_err[ifiber,2] = -1
                    self.z[ifiber,3] = -1
                    self.z_err[ifiber,3] = -1
                    self.z[ifiber,4] = -1
                    self.z_err[ifiber,4] = -1
                elif len(zminvals) == 2:
                    self.z[ifiber,2] = -1
                    self.z_err[ifiber,2] = -1
                    self.z[ifiber,3] = -1
                    self.z_err[ifiber,3] = -1
                    self.z[ifiber,4] = -1
                    self.z_err[ifiber,4] = -1
                elif len(zminvals) == 3:
                    self.z[ifiber,3] = -1
                    self.z_err[ifiber,3] = -1
                    self.z[ifiber,4] = -1
                    self.z_err[ifiber,4] = -1
                elif len(zminvals) == 4:
                    self.z[ifiber,4] = -1
                    self.z_err[ifiber,4] = -1
                imin = 0
                while (self.z[ifiber,1] == 0) | (self.z[ifiber,2] == 0) | (self.z[ifiber,3] == 0) | (self.z[ifiber,4] == 0):
                    imin += 1
                    nextpos = zminlocs[n.where(zminvals == n.sort(zminvals)[imin])[0][0]]
                    if (abs(posinvec - nextpos) > width) & (nextpos != 0) & (nextpos != bestzvec.shape[0]-1) :
                        xp = n.linspace(self.zbase[nextpos-1], self.zbase[nextpos+1], 1000)
                        f = quadfit(self.zbase[nextpos-1:nextpos+2], bestzvec[nextpos-1:nextpos+2])
                        fit = quad_for_fit(xp, f[0], f[1], f[2])
                        self.z[ifiber,imin] = xp[n.where(fit == n.min(fit))[0][0]]
                        self.z_err[ifiber,imin] = self.estimate_z_err(xp, fit)
                    else:
                        self.z[ifiber,imin] = -1.
                        self.z_err[ifiber,imin] = -1.
                
                self.flag_small_dchi2(ifiber, bestzvec, threshold=threshold, width=width) # Flag fibers with small delta chi2 in redshift


    def estimate_z_err(self, xp, fit):
        fitminloc = n.where(fit == n.min(fit))[0]
        if len( fitminloc ) is not 1: fitminloc = fitminloc[0] # Rarely, fit will have two equal global minima - use the first
        z_err = abs(xp[fitminloc]-xp[abs(n.min(fit)+1-fit).argmin()]) # abs() of difference between z_(chi2_min) and z_(chi2_min_+1)
        return z_err
    

    def flag_small_dchi2(self, ifiber, zvector, threshold, width): # zvector: vector of minimum chi2 in parameter-space at each redshift
        flag_val = int('0b100',2) # From BOSS zwarning flag definitions
        do_flag = False
        globminloc = n.where(zvector == n.min(zvector))[0][0]
        globmin = zvector[globminloc]
        zspline = gs.GridSpline(zvector)
        zminlocs = n.round(zspline.get_min())
        zminvals = zspline.get_val(zminlocs)
        small_dchi2 = n.where(zminvals < (globmin+threshold))[0]
        if len(small_dchi2) > 0:
            for i in small_dchi2:
                #if abs(zminvals[i] - globminloc) < threshold: do_flag = False
                if abs(zminlocs[i] - globminloc) > width: do_flag = True
        if do_flag: self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val


    def flag_z_fitlimit(self, ifiber):
        flag_val = int('0b100000',2) # From BOSS zwarning flag definitions
        self.zwarning[ifiber] = int(self.zwarning[ifiber]) ^ flag_val



#-------------------------------------------------------------------------------------------------------------------


    def z_refine2(self, threshold=46.6, width=15, num_z=5): # Default threshold of 46.6 is delta (chi_r)**2 = .01, and width is 1000 km/s
        self.num_z = num_z
        self.threshold = threshold
        self.width = width
        self.z = n.zeros((self.zchi2.shape[0],num_z))
        self.z_err = n.zeros((self.zchi2.shape[0],num_z))
        self.minvectors = []
        self.chi2vals = []
        for ifiber in xrange(self.zchi2.shape[0]):
            allminvectors = []
            bestminvectors = []
            bestchi2vals = []
            bestzvec = n.zeros( self.zchi2.shape[-1] )
            for iz in xrange(self.zchi2.shape[-1]): # Make vector of minimum chi2 at each trial redshift
                bestzvec[iz] = n.min( self.zchi2[ifiber,...,iz] )
                allminvectors.append( n.unravel_index( self.zchi2[ifiber,...,iz].argmin(), self.zchi2[ifiber,...,iz].shape ) ) # Location in zchi2[ifiber,...,iz] of min
            zspline = gs.GridSpline(bestzvec) # Spline up bestzvec
            zminlocs = map(int,n.round(zspline.get_min())) # Minimun locations
            zminvals = zspline.get_val(zminlocs) # Minimum values
            z_ind = 0
            stop = False
            while (z_ind < num_z) & (stop == False):
                thisminloc = zminvals.argmin() # Location in zminvals vector of minimum
                posinvec = zminlocs[thisminloc] # Location in bestzvec of minimum
                if (posinvec == 0) or (posinvec == bestzvec.shape[0]-1): # Flag and skip interpretation if best fit chi2 is at edge of z-range
                    if z_ind == 0: # If it's the first redshift, set all num_z redshifts and errors to -1
                        self.flag_z_fitlimit(ifiber)
                        self.z[ifiber] = -1.
                        self.z_err[ifiber] = -1.
                        self.minvectors.append( [(-1,)]*num_z )
                        self.chi2vals.append( [1e9]*num_z )
                        stop = True
                    elif posinvec == 0: # If it's not the first, and if it's on the left edge, set the 15 pixels to the right = 1e9 and move on
                        k = 0
                        while True:
                            k += 1
                            try:
                                if n.abs(zminlocs[thisminloc+k]-zminlocs[thisminloc]) < self.width:
                                    zminvals[thisminloc+k] = 1e9
                                else:
                                    break
                            except:
                                break
                        z_ind += 1
                    elif posinvec == bestzvec.shape[0]-1: # If it's not the first, and if it's on the right edge, set the 15 pixels to the left = 1e9 and move on
                        k = 0
                        while True:
                            k += 1
                            try:
                                if n.abs(zminlocs[thisminloc-k]-zminlocs[thisminloc]) < self.width:
                                    zminvals[thisminloc-k] = 1e9
                                else:
                                    break
                            except:
                                break
                        z_ind += 1
                else: # Fit with quadratic and find minimum and error based on that
                    xp = n.linspace(self.zbase[posinvec-1], self.zbase[posinvec+1], 1000)
                    f = quadfit(self.zbase[posinvec-1:posinvec+2], bestzvec[posinvec-1:posinvec+2])
                    fit = quad_for_fit(xp, f[0], f[1], f[2])
                    if not xp[n.where(fit == n.min(fit))[0][0]] in self.z[ifiber]:
                        self.z[ifiber,z_ind] = xp[n.where(fit == n.min(fit))[0][0]]
                        self.z_err[ifiber,z_ind] = self.estimate_z_err(xp, fit)
                        bestchi2vals.append( bestzvec[posinvec] )
                        bestminvectors.append( allminvectors[posinvec] + (posinvec,) )
                    else:
                        self.z[ifiber,z_ind] = -1.
                        self.z_err[ifiber,z_ind] = -1.
                        bestchi2vals.append(1e9)
                        bestminvectors.append((-1,))
                    zminvals[thisminloc] = 1e9
                    #self.minvectors.append( n.unravel_index(self.zchi2[ifiber,...,posinvec].argmin(), self.zchi2[ifiber,...,posinvec].shape) + (posinvec,) )
                    k = 0
                    while True: # Set width points to the right and left of minimum to 1e9
                        k += 1
                        try:
                            if n.abs(zminlocs[thisminloc+k]-zminlocs[thisminloc]) < self.width:
                                zminvals[thisminloc+k] = 1e9
                            else:
                                break
                        except:
                            break
                    k = 0
                    while True:
                        k += 1
                        try:
                            if n.abs(zminlocs[thisminloc-k]-zminlocs[thisminloc]) < self.width:
                                zminvals[thisminloc-k] = 1e9
                            else:
                                break
                        except:
                            break
                    z_ind += 1
            self.minvectors.append( bestminvectors )
            self.chi2vals.append( bestchi2vals )
        self.flag_small_dchi2_2()


    def flag_small_dchi2_2(self):
        flag_val = int('0b100',2) # From BOSS zwarning flag definitions
        if self.num_z > 1:
            for ifiber in xrange(self.zchi2.shape[0]):
                if n.abs(self.chi2vals[ifiber][1] - self.chi2vals[ifiber][0]) < self.threshold:
                    self.zwarning[ifiber] = int(self.zwarning[ifiber]) | flag_val

            



            
'''
            
            
            
            
            
            self.minvector.append( (ifiber,) + n.unravel_index(self.zchi2[ifiber].argmin(),self.zchi2[ifiber].shape))
            bestzvec = n.zeros( self.zchi2.shape[-1])
            for iz in xrange(self.zchi2.shape[-1]):
                bestzvec[iz] = n.min( self.zchi2[ifiber,...,iz] )
            posinvec = n.where( bestzvec == n.min(bestzvec) )[0][0]
            if (posinvec == 0) or (posinvec == bestzvec.shape[0]-1): # Flag and skip interpolation fit if best chi2 is at edge of z-range
                self.flag_z_fitlimit(ifiber)
                self.z[ifiber] = -1.
                self.z_err[ifiber] = -1.
            else:
                # Find best redshift using global chi2 minimum
                xp = n.linspace(self.zbase[posinvec-1], self.zbase[posinvec+1], 1000)
                f = quadfit(self.zbase[posinvec-1:posinvec+2], bestzvec[posinvec-1:posinvec+2])
                fit = quad_for_fit(xp, f[0], f[1], f[2])
                #p.plot(xp, fit, color='red', hold=False)
                #p.plot(self.zbase[posinvec-1:posinvec+2], bestzvec[posinvec-1:posinvec+2], 'ko', hold=True)
                self.z[ifiber,0] = xp[n.where(fit == n.min(fit))[0][0]]
                self.z_err[ifiber,0] = self.estimate_z_err(xp, fit)
                # Find second- and third-, and fourth-best redshifts, just to have them
                zspline = gs.GridSpline(bestzvec)
                zminlocs = n.round(zspline.get_min())
                zminvals = zspline.get_val(zminlocs)
                if len(zminvals) == 1:
                    self.z[ifiber,1] = -1.
                    self.z_err[ifiber,1] = -1.
                    self.z[ifiber,2] = -1
                    self.z_err[ifiber,2] = -1
                    self.z[ifiber,3] = -1
                    self.z_err[ifiber,3] = -1
                    self.z[ifiber,4] = -1
                    self.z_err[ifiber,4] = -1
                elif len(zminvals) == 2:
                    self.z[ifiber,2] = -1
                    self.z_err[ifiber,2] = -1
                    self.z[ifiber,3] = -1
                    self.z_err[ifiber,3] = -1
                    self.z[ifiber,4] = -1
                    self.z_err[ifiber,4] = -1
                elif len(zminvals) == 3:
                    self.z[ifiber,3] = -1
                    self.z_err[ifiber,3] = -1
                    self.z[ifiber,4] = -1
                    self.z_err[ifiber,4] = -1
                elif len(zminvals) == 4:
                    self.z[ifiber,4] = -1
                    self.z_err[ifiber,4] = -1
                imin = 0
                while (self.z[ifiber,1] == 0) | (self.z[ifiber,2] == 0) | (self.z[ifiber,3] == 0) | (self.z[ifiber,4] == 0):
                    imin += 1
                    nextpos = zminlocs[n.where(zminvals == n.sort(zminvals)[imin])[0][0]]
                    if (abs(posinvec - nextpos) > width) & (nextpos != 0) & (nextpos != bestzvec.shape[0]-1) :
                        xp = n.linspace(self.zbase[nextpos-1], self.zbase[nextpos+1], 1000)
                        f = quadfit(self.zbase[nextpos-1:nextpos+2], bestzvec[nextpos-1:nextpos+2])
                        fit = quad_for_fit(xp, f[0], f[1], f[2])
                        self.z[ifiber,imin] = xp[n.where(fit == n.min(fit))[0][0]]
                        self.z_err[ifiber,imin] = self.estimate_z_err(xp, fit)
                    else:
                        self.z[ifiber,imin] = -1.
                        self.z_err[ifiber,imin] = -1.
                
                self.flag_small_dchi2(ifiber, bestzvec, threshold=threshold, width=width) # Flag fibers with small delta chi2 in redshift
'''
#-------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------

def quad_for_fit(x, a, b, c):
    return a*(x**2) + b*x + c
