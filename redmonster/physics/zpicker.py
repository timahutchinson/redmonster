# Use chi2 surfaces from zfinder.py to choose and classify objects on each fiber.  Each objclass
# argument is the ENTIRE object created by Zfinder for a single template.  Output is in
# zpick.type, a list of len()=nfibers where each entry is a string of objclassN.type for the
# best reduced chi2 amongst all input classes.
#
# Tim Hutchinson, July 2014
# t.hutchinson@utah.edu

import numpy as n

class Zpicker:

    def __init__(self, npixflux, objclass1=None, objclass2=None, objclass3=None, objclass4=None, objclass5=None):
        self.type = []
        self.subtype = []
        self.npixflux = npixflux
        self.objclass1 = objclass1
        if not objclass2: self.nclass = 1
        elif objclass2 and not objclass3: nclass = 2
        elif objclass3 and not objclass4: nclass = 3
        elif objclass4 and not objclass5: nclass = 4
        else: nclass = 5
        self.minrchi2 = n.zeros( (objclass1.zchi2arr.shape[0],nclass) )
        self.classify_obj(objclass1, objclass2, objclass3, objclass4, objclass5)
    
    def classify_obj(self, objclass1, objclass2, objclass3, objclass4, objclass5):
        for ifiber in xrange(objclass1.zchi2arr.shape[0]):
            self.minrchi2[ifiber,0] = n.min(objclass1.zchi2arr[ifiber]) / (self.npixflux - objclass1.npoly)
            if objclass2: self.minrchi2[ifiber,1] = n.min(objclass2.zchi2arr[ifiber]) / (self.npixflux - objclass2.npoly)
            if objclass3: self.minrchi2[ifiber,2] = n.min(objclass3.zchi2arr[ifiber]) / (self.npixflux - objclass3.npoly)
            if objclass4: self.minrchi2[ifiber,3] = n.min(objclass4.zchi2arr[ifiber]) / (self.npixflux - objclass4.npoly)
            if objclass5: self.minrchi2[ifiber,4] = n.min(objclass5.zchi2arr[ifiber]) / (self.npixflux - objclass5.npoly)
            minpos = self.minrchi2[ifiber].argmin()
            if minpos == 0: self.type.append(objclass1.type)
            elif minpos == 1: self.type.append(objclass2.type)
            elif minpos == 2: self.type.append(objclass3.type)
            elif minpos == 3: self.type.append(objclass4.type)
            elif minpos == 4: self.type.append(objclass5.type)
