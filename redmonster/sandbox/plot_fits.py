# GUI used for quickly plotting BOSS spectra.  Also allows overplotting of best-fit template as
# determined by redmonster pipeline.  Sort of a redmonster version of plotspec.pro, though currently
# with less bells and whistles.
#
# Tim Hutchinson, University of Utah, April 2014
# Signifcantly updated by TH, October 2014
#
# thutchinson@utah.edu

from Tkinter import *
import numpy as n
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from astropy.io import fits
from os import environ
from os.path import join, exists
from redmonster.math.misc import poly_array
from astropy.convolution import convolve, Box1DKernel
'''
from __future__ import print_function
import glob
import os
os.chdir("/mydir")
for file in glob.glob("*.txt"):
    print(file)
'''

class Plot_Fit:

    def __init__(self, master):
        self.master = master
        self.plate = None
        self.mjd = None
        L1 = Label(master, text='Plate')
        L1.grid(sticky=E)
        L2 = Label(master, text='MJD')
        L2.grid(sticky=E)
        L3 = Label(master, text='Fiber')
        L3.grid(stick=E)
        self.e1 = Entry(master)
        self.e1.grid(row=0, column=1)
        self.e2 = Entry(master)
        self.e2.grid(row=1, column=1)
        self.e3 = Entry(master)
        self.e3.grid(row=2, column=1)
        self.var = IntVar()
        c = Checkbutton(master, text='Overplot best-fit model', variable=self.var)
        c.grid(row=3, column=1)
        plot = Button(master, text='Plot', command=self.do_plot)
        plot.grid(row=4, column=1)
        qbutton = Button(master, text='QUIT', fg='red', command=master.destroy)
        qbutton.grid(row=5, column=1)
        nextfiber = Button(master, text='>', command=self.next_fiber)
        nextfiber.grid(row=1, column=4)
        prevfiber = Button(master, text='<', command=self.prev_fiber)
        prevfiber.grid(row=1, column=3)

    def do_plot(self):
        if self.plate != int(self.e1.get()) or self.mjd != int(self.e2.get()):
            self.plate = int(self.e1.get())
            self.mjd = int(self.e2.get())
            self.fiber = int(self.e3.get())
            self.platepath = join(environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'], '%s' % self.plate, 'spPlate-%s-%s.fits' % (self.plate, self.mjd))
            hdu = fits.open(self.platepath)
            self.specs = hdu[0].data
            self.wave = 10**(hdu[0].header['COEFF0'] + n.arange(hdu[0].header['NAXIS1'])*hdu[0].header['COEFF1'])
            self.models = fits.open(join(environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'], '%s' % self.plate, environ['RUN1D'], 'redmonster-%s-%s.fits' % (self.plate, self.mjd)))[2].data
            self.fiberid = fits.open(join(environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'], '%s' % self.plate, environ['RUN1D'], 'redmonster-%s-%s.fits' % (self.plate, self.mjd)))[1].data.FIBERID
        else:
            self.fiber = int(self.e3.get())
        f = Figure(figsize=(10,6), dpi=100)
        a = f.add_subplot(111)
        if self.var.get() == 0:
            a.plot(self.wave, self.specs[self.fiber], color='red')
        elif self.var.get() == 1:
            a.plot(self.wave, convolve(self.specs[self.fiber], Box1DKernel(5)), color='red')
            
            # Overplot model
            loc = n.where(self.fiberid == self.fiber)[0]
            if len(loc) is not 0:
                a.plot(self.wave, self.models[loc[0]], color='black')
            else:
                print 'Fiber %s is not in redmonster-%s-%s.fits' % (self.fiber, self.plate, self.mjd)

        a.set_xlabel('Wavelength (Angstroms)')
        a.set_ylabel('Flux in some units')
        canvas = FigureCanvasTkAgg(f, master=self.master)
        canvas.show()
        canvas.get_tk_widget().grid(row=5, column=5)

    def next_fiber(self):
        self.fiber += 1
        self.e3.delete(0, END)
        self.e3.insert(0, str(self.fiber))
        self.do_plot()

    def prev_fiber(self):
        self.fiber -= 1
        self.e3.delete(0, END)
        self.e3.insert(0, str(self.fiber))
        self.do_plot()

root = Tk()
app = Plot_Fit(root)
root.mainloop()
