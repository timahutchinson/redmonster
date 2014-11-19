# Write redmonster output files
#
# Tim Hutchinson, , University of Utah, April 2014
# t.hutchinson@utah.edu
#
# July 2014 - now probably defunct, replaced by redmonster.datamgr.io.py

class Output:

    def __init__(self, spec=None):
        if spec: self.set_filenames()

    def set_filenames(self):
        self.zallfile = 'spZall-%s-%s.fits' % (spec.plate, spec.mjd)
        self.zbestfile = 'spZbest-%s-%s.fits' % (spec.plate, spec.mjd)
        self.zlinefile = 'spZline-%s-%s.fits' % (spec.plate, spec.mjd)
        self.logfile = 'spDiag1d-%s-%s.log' % (spec.plate, spec.mjd)
        self.plotfile = 'spDiag1d-%s-%s.ps' % (spec.plate, spec.mjd)