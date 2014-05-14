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
from plot_fits import Plot_fit

zchi2_ssp = fits.open('/Users/Tim/zchi2_ssp.fits')[0].data
zchi2_star = fits.open('/Users/Tim/zchi2_star.fits')[0].data

root = Tk()
app = Plot_fit(root, zchi2_ssp, zchi2_star)
#app.mainloop()