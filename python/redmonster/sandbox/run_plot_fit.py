# Command line execution for plot_fit.py
#
# Tim Hutchinson, University of Utah, April 2014
# t.hutchinson@utah.edu

from Tkinter import *
import numpy as n
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from astropy.io import fits
from os import environ
from os.path import join
from redmonster.physics.misc import poly_array
from plot_fits import PlotFit

zchi2_ssp = fits.open('/Users/Tim/zchi2_ssp.fits')[0].data
zchi2_star = fits.open('/Users/Tim/zchi2_star.fits')[0].data
zchi2_cap = fits.open('/Users/Tim/zchi2_cap.fits')[0].data

root = Tk()
#app = PlotFit(root, zchi2_star=zchi2_star, zchi2_cap=zchi2_cap, fiber_offset=39)
app = PlotFit(root, zchi2_ssp=zchi2_ssp, fiber_offset=39)
#NavigationToolbar2TkAgg(app, root)
#app.mainloop()
