# Import statements

class Output:

    def __init__(self, spec=None):
        if spec: self.set_filenames()

    def set_filenames(self):
        self.zallfile = 'spZall-%s-%s.fits' % (spec.plate, spec.mjd)
        self.zbestfile = 'spZbest-%s-%s.fits' % (spec.plate, spec.mjd)
        self.zlinefile = 'spZline-%s-%s.fits' % (spec.plate, spec.mjd)
        self.logfile = 'spDiag1d-%s-%s.log' % (spec.plate, spec.mjd)
        self.plotfile = 'spDiag1d-%s-%s.ps' % (spec.plate, spec.mjd)