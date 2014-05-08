# Routine to refine redshifts found in zfinder.py
import numpy as n

class Zfitter:

    def __init__(self, zchi2):
        self.zchi2 = zchi2

    def z_refine(self):
        zminpos = n.where(self.zchi2 == n.min(self.zchi2))
        vecpos = ()
        for i in xrange(len(zminpos)-1):
            vecpos += (pos[i][0],)
        self.bestzvec = self.zchi2[vecpos]


