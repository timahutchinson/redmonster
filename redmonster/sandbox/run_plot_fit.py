from Tkinter import *
import numpy as n
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from astropy.io import fits
from os import environ
from os.path import join
from redmonster.math.misc import poly_array
from plot_fits import Plot_Fit

zchi2_ssp = fits.open('/Users/Tim/zchi2_ssp.fits')[0].data
zchi2_star = fits.open('/Users/Tim/zchi2_star.fits')[0].data
zchi2_cap = fits.open('/Users/Tim/zchi2_cap.fits')[0].data

root = Tk()
app = Plot_fit(root, zchi2_star=zchi2_star, zchi2_cap=zchi2_cap, fiber_offset=39)
#NavigationToolbar2TkAgg(app, root)
#app.mainloop()