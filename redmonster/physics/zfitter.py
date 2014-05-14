# Routine to refine redshifts found in zfinder.py

import numpy as n
#from scipy.optimize import curve_fit
from redmonster.math.misc import quadfit
import matplotlib as m
from matplotlib import pyplot as p
m.interactive(True)

class Zfitter:

    def __init__(self, zchi2, zbase):
        self.zchi2 = zchi2
        self.zbase = zbase
        self.best_z = n.zeros(zchi2.shape[0])
        self.z_err = n.zeros(zchi2.shape[0])
    #import pdb; pdb.set_trace()

    def z_refine(self):
        for ifiber in xrange(self.zchi2.shape[0]):
            zminpos = n.where(self.zchi2[ifiber] == n.min(self.zchi2[ifiber]))
            vecpos = ()
            for i in xrange(len(zminpos)-1):
                vecpos += (zminpos[i][0],)
            bestzvec = self.zchi2[(ifiber,)+vecpos]
            posinvec = zminpos[-1][0]
            #popt, pcov = curve_fit(quad_for_fit, self.zbase[posinvec-1:posinvec+2]*1000., bestzvec[posinvec-1:posinvec+2]) # Factor of 1000 because numerical curve_fit sometimes doesn't like small numbers
            xp = n.linspace(self.zbase[posinvec-1], self.zbase[posinvec+1], 1000)
            #fit = quad_for_fit(xp, popt[0]*(1000.**2), popt[1]*(1000.), popt[2]) # See above comment
            f = quadfit(self.zbase[posinvec-1:posinvec+2], bestzvec[posinvec-1:posinvec+2])
            fit = quad_for_fit(xp, f[0], f[1], f[2])
            #p.plot(xp, fit, color='red')
            #p.plot(self.zbase[posinvec-1:posinvec+2], bestzvec[posinvec-1:posinvec+2], 'ko', hold=True)
            self.best_z[ifiber] = xp[n.where(fit == n.min(fit))[0][0]]
            self.z_err[ifiber] = self.estimate_z_err(xp, fit)
            print self.best_z[ifiber]

    def estimate_z_err(self, xp, fit):
        fitminloc = n.where(fit == n.min(fit)) # Index of lowest chi2
        z_err = abs(xp[fitminloc]-xp[abs(n.min(fit)+1-fit).argmin()]) # abs() of difference between z_(chi2_min) and z_(chi2_min_+1)
        return z_err

# -----------------------------------------------------------------------------

def quad_for_fit(x, a, b, c):
    return a*(x**2) + b*x + c