class spec:

    def __init__(self,plate=None,mjd=None):
        self.plate = plate
        self.mjd = mjd
        self.flux = None if self.mjd and self.plate else self.set_flux()
        
    def set_flux(self):
        #readspec
        self.flux = 0