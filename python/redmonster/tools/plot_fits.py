# GUI used for quickly plotting BOSS spectra.  Also allows overplotting of best-fit template as
# determined by redmonster pipeline.  Sort of a redmonster version of plotspec.pro, though currently
# with less bells and whistles.
#
# Tim Hutchinson, University of Utah, April 2014
# Signifcantly updated by TH, October 2014
#
# thutchinson@utah.edu

from os import environ
from os.path import join, exists
try:
    from tkinter import *
except ImportError:
    from Tkinter import *
import numpy as n
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
import seaborn as sns
sns.set_style('whitegrid')

from redmonster.physics.misc import poly_array


class PlotFit(Frame):
    def __init__ (self):
        self.root = Tk()
        self.ablinelist = [
                           3890.2, 3933.7, 3968.5, 4102.9, 4307, 4341.7, 4862.7,
                           5175, 5889, 5895
                           ]
        self.ablinenames = [
                            r'H$\epsilon$','Ca K', 'Ca H', r'H$\delta$', 'Ca G',
                            r'H$\gamma$', r'H$\beta$', 'Mg I', 'Na I', 'Na I'
                            ]
        self.emlinelist = [2500]
        self.emlinenames = ['test emline']
        self.plate = None
        self.mjd = None
        #
        plate = StringVar()
        plate.set('7848')
        mjd = StringVar()
        mjd.set('56959')
        #
        L1 = Label(self.root, text='Plate')
        L1.grid(sticky=E)
        L2 = Label(self.root, text='MJD')
        L2.grid(sticky=E)
        L3 = Label(self.root, text='Fiber')
        L3.grid(stick=E)
        L5 = Label(self.root, text='z num')
        L5.grid(stick=E)
        self.e1 = Entry(self.root, textvariable=plate)
        self.e1.bind()
        self.e1.grid(row=0, column=1)
        self.e2 = Entry(self.root, textvariable=mjd)
        self.e2.grid(row=1, column=1)
        fiber = StringVar()
        fiber.set('0')
        self.e3 = Entry(self.root, textvariable=fiber)
        self.e3.grid(row=2, column=1)
        znum = StringVar()
        znum.set('1')
        self.e5 = Entry(self.root, textvariable=znum)
        self.e5.grid(row=3, column=1)
        nextz = Button(self.root, text='+', command=self.next_z)
        nextz.grid(row=3, column=4)
        prevz = Button(self.root, text='-', command=self.prev_z)
        prevz.grid(row=3, column=3)
        self.var = BooleanVar()
        self.var.set(1)
        self.restframe = BooleanVar()
        self.restframe.set(0)
        self.ablines = BooleanVar()
        self.ablines.set(0)
        self.emlines = BooleanVar()
        self.emlines.set(0)
        c = Checkbutton(self.root, text='Overplot best-fit model',
                        variable=self.var)
        c.grid(row=4, column=1)
        restframe = Checkbutton(self.root, text='Rest-frame wavelength',
                                variable=self.restframe)
        restframe.grid(row=5,column=1)
        ablines = Checkbutton(self.root, text='Show absorption lines ',
                              variable=self.ablines)
        ablines.grid(row=6, column=1)
        emlines = Checkbutton(self.root, text='Show emission lines    ',
                              variable=self.emlines)
        emlines.grid(row=7, column=1)
        #
        smooth = StringVar()
        smooth.set('5')
        L4 = Label(self.root, text='Smooth')
        L4.grid(sticky=E)
        self.e4 = Entry(self.root, textvariable=smooth)
        self.e4.grid(row=8, column=1)
        plot = Button(self.root, text='Plot', command=self.do_plot)
        plot.grid(row=9, column=1)
        qbutton = Button(self.root, text='QUIT', fg='red',
                         command=self.root.destroy)
        qbutton.grid(row=10, column=1)
        nextfiber = Button(self.root, text='>', command=self.next_fiber)
        nextfiber.grid(row=2, column=4)
        prevfiber = Button(self.root, text='<', command=self.prev_fiber)
        prevfiber.grid(row=2, column=3)
        Frame.__init__(self,self.root)
        self.root.mainloop()

    def do_plot(self):
        if self.plate != int(self.e1.get()) or self.mjd != int(self.e2.get()):
            self.plate = int(self.e1.get())
            self.mjd = int(self.e2.get())
            self.fiber = int(self.e3.get())
            self.znum = int(self.e5.get())
            self.platepath = join(environ['BOSS_SPECTRO_REDUX'],
                                  environ['RUN2D'], '%s' % self.plate,
                                  'spPlate-%s-%s.fits' % (self.plate, self.mjd))
            hdu = fits.open(self.platepath)
            self.specs = hdu[0].data
            self.wave = 10**(hdu[0].header['COEFF0'] +
                             n.arange(hdu[0].header['NAXIS1']) *
                             hdu[0].header['COEFF1'])
            self.models = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                         environ['RUN2D'], '%s' % self.plate,
                                         environ['RUN1D'],
                                         'redmonster-%s-%s.fits' %
                                         (self.plate, self.mjd)))[2].data
            self.fiberid = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                          environ['RUN2D'], '%s' % self.plate,
                                          environ['RUN1D'],
                                          'redmonster-%s-%s.fits' %
                                          (self.plate,
                                           self.mjd)))[1].data.FIBERID
            self.type1 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                        environ['RUN2D'], '%s' % self.plate,
                                        environ['RUN1D'],
                                        'redmonster-%s-%s.fits' %
                                        (self.plate, self.mjd)))[1].data.CLASS1
            self.type2 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                        environ['RUN2D'], '%s' % self.plate,
                                        environ['RUN1D'],
                                        'redmonster-%s-%s.fits' %
                                        (self.plate, self.mjd)))[1].data.CLASS2
            self.type3 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                        environ['RUN2D'], '%s' % self.plate,
                                        environ['RUN1D'],
                                        'redmonster-%s-%s.fits' %
                                        (self.plate, self.mjd)))[1].data.CLASS3
            self.type4 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                        environ['RUN2D'], '%s' % self.plate,
                                        environ['RUN1D'],
                                        'redmonster-%s-%s.fits' %
                                        (self.plate, self.mjd)))[1].data.CLASS4
            self.type5 = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                        environ['RUN2D'], '%s' % self.plate,
                                        environ['RUN1D'],
                                        'redmonster-%s-%s.fits' %
                                        (self.plate, self.mjd)))[1].data.CLASS5
            self.z = n.zeros((self.fiberid.shape[0],5))
            self.z[:,0] = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                         environ['RUN2D'], '%s' % self.plate,
                                         environ['RUN1D'],
                                         'redmonster-%s-%s.fits' %
                                         (self.plate, self.mjd)))[1].data.Z1
            self.z[:,1] = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                         environ['RUN2D'], '%s' % self.plate,
                                         environ['RUN1D'],
                                         'redmonster-%s-%s.fits' %
                                         (self.plate, self.mjd)))[1].data.Z2
            self.z[:,2] = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                         environ['RUN2D'], '%s' % self.plate,
                                         environ['RUN1D'],
                                         'redmonster-%s-%s.fits' %
                                         (self.plate, self.mjd)))[1].data.Z3
            self.z[:,3] = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                         environ['RUN2D'], '%s' % self.plate,
                                         environ['RUN1D'],
                                         'redmonster-%s-%s.fits' %
                                         (self.plate, self.mjd)))[1].data.Z4
            self.z[:,4] = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                         environ['RUN2D'], '%s' % self.plate,
                                         environ['RUN1D'],
                                         'redmonster-%s-%s.fits' %
                                         (self.plate, self.mjd)))[1].data.Z5
            self.zwarning = fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],
                                           environ['RUN2D'], '%s' % self.plate,
                                           environ['RUN1D'],
                                           'redmonster-%s-%s.fits' %
                                           (self.plate,
                                            self.mjd)))[1].data.ZWARNING
        else:
            self.fiber = int(self.e3.get())
        f = Figure(figsize=(10,6), dpi=100)
        a = f.add_subplot(111)
        loc = n.where(self.fiberid == self.fiber)[0]
        if self.znum == 1:
            z = self.z[loc[0],0]
            thistype = self.type1[loc[0]]
        elif self.znum == 2:
            z = self.z[loc[0],1]
            thistype = self.type2[loc[0]]
        elif self.znum == 3:
            z = self.z[loc[0],2]
            thistype = self.type3[loc[0]]
        elif self.znum == 4:
            z = self.z[loc[0],3]
            thistype = self.type4[loc[0]]
        elif self.znum == 5:
            z = self.z[loc[0],4]
            thistype = self.type5[loc[0]]
        if self.var.get() == 0:
            if self.restframe.get() == 0:
                a.plot(self.wave, self.specs[self.fiber], color='black')
            elif self.restframe.get() == 1:
                a.plot(self.wave/(1+self.z[loc][0]), self.specs[self.fiber],
                       color='black')
        elif self.var.get() == 1:
            smooth = self.e4.get()
            if smooth is '':
                if self.restframe.get() == 0:
                    a.plot(self.wave, self.specs[self.fiber], color='black')
                elif self.restframe.get() == 1:
                    a.plot(self.wave/(1+z), self.specs[self.fiber],
                           color='black')
            else:
                if self.restframe.get() == 0:
                    a.plot(self.wave, convolve(self.specs[self.fiber],
                                               Box1DKernel(int(smooth))),
                           color='black')
                elif self.restframe.get() == 1:
                    a.plot(self.wave/(1+z), convolve(self.specs[self.fiber],
                                                     Box1DKernel(int(smooth))),
                           color='black')
            # Overplot model
            if len(loc) is not 0:
                if self.restframe.get() == 0:
                    #a.plot(self.wave, self.models[loc[0]], color='black')
                    # This for when multiple models are in redmonster file
                    a.plot(self.wave, self.models[loc[0],self.znum-1],
                           color='cyan')
                    if self.ablines.get() == 1:
                        for i, line in enumerate(self.ablinelist):
                            if ((line*(1+z) > self.wave[0]) &
                                    (line*(1+z) < self.wave[-1])):
                                a.axvline(line*(1+z), color='blue',
                                          linestyle='--',
                                          label=self.ablinenames[i])
                    if self.emlines.get() == 1:
                        for i, line in enumerate(self.emlinelist):
                            if (line*(1+z) > self.wave[0]) & (line*(1+z) < \
                                                              self.wave[-1]):
                                a.axvline(line*(1+z), color='red',
                                          linestyle='--',
                                          label=self.emlinenames[i])
                    if self.ablines.get() == 1 or self.emlines.get() == 1:
                        a.legend(prop={'size':10})
                elif self.restframe.get() == 1:
                    a.plot(self.wave/(1+z), self.models[loc[0],self.znum-1],
                           color='cyan')
                    if self.ablines.get() == 1:
                        for i, line in enumerate(self.ablinelist):
                            if (line > self.wave[0]) & (line < self.wave[-1]):
                                a.axvline(line, color='blue', linestyle='--',
                                          label=self.ablinenames[i])
                    if self.emlines.get() == 1:
                        for i, line in enumerate(self.emlinelist):
                            if (line > self.wave[0]) & (line < self.wave[-1]):
                                a.axvline(line, color='red', linestyle='--',
                                          label=self.emlinenames[i])
                    if self.ablines.get() == 1 or self.emlines.get() == 1:
                        a.legend(prop={'size':10})
                a.set_title('Plate %s Fiber %s: z=%s class=%s zwarning=%s' %
                            (self.plate, self.fiber, z, thistype,
                             self.zwarning[loc[0]]))
            else:
                print('Fiber %s is not in redmonster-%s-%s.fits' % \
                        (self.fiber, self.plate, self.mjd))
                a.set_title('Plate %s Fiber %s' % (self.plate, self.fiber))

        if self.restframe.get() == 1:
            lower_data, upper_data = self.set_limits()
            a.axis([self.wave[0]/(1+z)-100,self.wave[-1]/(1+z)+100,
                    lower_data,upper_data])
        elif self.restframe.get() == 0:
            lower_data, upper_data = self.set_limits()
            a.axis([self.wave[0]-100,self.wave[-1]+100,lower_data,upper_data])
        a.set_xlabel('Wavelength ($\AA$)')
        a.set_ylabel('Flux ($10^{-17} erg\ cm^2 s^{-1} \AA^{-1}$)')
        canvas = FigureCanvasTkAgg(f, master=self.root)
        canvas.get_tk_widget().grid(row=0, column=5, rowspan=20)
        toolbar_frame = Frame(self.root)
        toolbar_frame.grid(row=20,column=5)
        toolbar = NavigationToolbar2TkAgg( canvas, toolbar_frame )
        canvas.show()

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

    def next_z(self):
        if (self.znum >= 1) & (self.znum < 5):
            self.znum += 1
            self.e5.delete(0, END)
            self.e5.insert(0, str(self.znum))
            self.do_plot()
        else:
            if self.znum < 1:
                self.znum = 1
                self.e5.delete(0, END)
                self.e5.insert(0, str(self.znum))
                self.do_plot()
            elif self.znum >= 5:
                self.znum = 5
                self.e5.delete(0, END)
                self.e5.insert(0, str(self.znum))
                self.do_plot()
            else:
                self.znum = 1
                self.e5.delete(0, END)
                self.e5.insert(0, str(self.znum))
                self.do_plot()

    def prev_z(self):
        if (self.znum > 1) & (self.znum <= 5):
            self.znum -= 1
            self.e5.delete(0, END)
            self.e5.insert(0, str(self.znum))
            self.do_plot()
        else:
            if self.znum <= 1:
                self.znum = 1
                self.e5.delete(0, END)
                self.e5.insert(0, str(self.znum))
                self.do_plot()
            elif self.znum > 5:
                self.znum = 5
                self.e5.delete(0, END)
                self.e5.insert(0, str(self.znum))
                self.do_plot()
            else:
                self.znum = 1
                self.e5.delete(0, END)
                self.e5.insert(0, str(self.znum))
                self.do_plot()

    def set_limits(self, percentile=.95):
        sorted_flux = n.sort( self.specs[self.fiber] )
        bottom_ind = int(n.floor((1-percentile)/2. * sorted_flux.shape[0]))
        top_ind = n.ceil(sorted_flux.shape[0] - bottom_ind)
        return sorted_flux[bottom_ind], sorted_flux[top_ind]



app = PlotFit()
