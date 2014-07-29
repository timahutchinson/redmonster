# GUI used for quickly plotting BOSS spectra.  Also allows overplotting of best-fit template as
# determined by redmonster pipeline.  Sort of a redmonster version of plotspec.pro, though currently
# with less bells and whistles.
#
# Tim Hutchinson, April 2014
# thutchinson@utah.edu

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
from scipy.signal import boxcar

class Plot_Fit:

    def __init__(self, master, zchi2_ssp=None, zchi2_star=None, zchi2_cap=None, fiber_offset=0):
        self.master = master
        self.plate = None
        self.zchi2_ssp = zchi2_ssp
        self.zchi2_star = zchi2_star
        self.zchi2_cap = zchi2_cap
        self.fos = fiber_offset
        L1 = Label(master, text='Plate')
        L1.grid(sticky=E)
        L2 = Label(master, text='Fiber')
        L2.grid(sticky=E)
        self.e1 = Entry(master)
        self.e1.grid(row=0, column=1)
        self.e2 = Entry(master)
        self.e2.grid(row=1, column=1)
        self.var = IntVar()
        c = Checkbutton(master, text='Overplot best fit SSP and star', variable=self.var)
        c.grid(row=2, column=1)
        plot = Button(master, text='Plot', command=self.do_plot)
        plot.grid(row=3, column=1)
        qbutton = Button(master, text='QUIT', fg='red', command=master.destroy)
        qbutton.grid(row=4, column=1)
        nextfiber = Button(master, text='>', command=self.next_fiber)
        nextfiber.grid(row=1, column=4)
        prevfiber = Button(master, text='<', command=self.prev_fiber)
        prevfiber.grid(row=1, column=3)
        #self.startemps2 = fits.open(join(environ['REDMONSTER_DIR'],'templates','ndArch-spEigenStar-55734.fits'))[0].data
        #self.startemps1 = fits.open(join(environ['REDMONSTER_DIR'],'templates','ndArch-all-CAP-grids.fits'))[0].data
        self.ssptemps = fits.open(join(environ['REDMONSTER_DIR'],'templates','ndArch-ssp_em_galaxy-v000.fits'))[0].data

    def do_plot(self):
        if self.plate != int(self.e1.get()):
            self.plate = int(self.e1.get())
            self.fiber = int(self.e2.get())
            self.platepath = '/Users/Tim/Documents/Workstuff/BOSS/redux/run2d/4389/spPlate-4389-55539.fits'
            hdu = fits.open(self.platepath)
            self.specs = hdu[0].data
            self.ivar = hdu[1].data
            self.wave = 10**(hdu[0].header['COEFF0'] + n.arange(4663)*hdu[0].header['COEFF1'])
        else:
            self.fiber = int(self.e2.get())
        f = Figure(figsize=(10,6), dpi=100)
        a = f.add_subplot(111)
        if self.var.get() == 0:
            a.plot(self.wave, self.specs[self.fiber], color='red')
        elif self.var.get() == 1:
            a.plot(self.wave, n.convolve(self.specs[self.fiber], boxcar(5)/5., mode='same'), color='red')
            
            # Overplot best-fit SSP
            ssploc = []
            for i in xrange(len(self.zchi2_ssp.shape)-1):
                ssploc.append(n.where(self.zchi2_ssp[self.fiber-self.fos] == n.min(self.zchi2_ssp[self.fiber-self.fos]))[i][0])
            pmat = n.zeros((4663,5))
            pmat[:,0] = self.ssptemps[tuple(ssploc[:-1])][ssploc[-1]+667:ssploc[-1]+4663+667]
            pmat[:,1:] = n.transpose(poly_array(4, 4663))
            ninv = n.diag(self.ivar[self.fiber])
            coeffvec = n.linalg.solve(n.dot(n.transpose(pmat),n.dot(ninv,pmat)),n.dot(n.transpose(pmat),n.dot(ninv,self.specs[self.fiber])))
            model1 = n.dot(pmat,coeffvec)
            a.plot(self.wave, model1, color='black')
            
            # Overplot best fit spEigenStar
            #ssploc = []
            #for i in xrange(len(self.zchi2_star.shape)-1):
            #    ssploc.append(n.where(self.zchi2_star[self.fiber-self.fos] == n.min(self.zchi2_star[self.fiber-self.fos]))[i][0])
            #pmat = n.zeros((4663,5))
            #pmat[:,0] = self.startemps2[tuple(ssploc[:-1])][ssploc[-1]+653:ssploc[-1]+4663+653] #653 for spEigenStar, 864 for CAP grids, 667 for Charlie's SSPs
            #pmat[:,1:] = n.transpose(poly_array(4, 4663))
            #ninv = n.diag(self.ivar[self.fiber])
            #coeffvec = n.linalg.solve(n.dot(n.transpose(pmat),n.dot(ninv,pmat)),n.dot(n.transpose(pmat),n.dot(ninv,self.specs[self.fiber])))
            #model2 = n.dot(pmat,coeffvec)
            #a.plot(self.wave, model2, color='blue')
            #a.set_ylim([n.min(model2)-5, n.max(model2)+5])
            
            # Overplot best-fit CAP star
            #ssploc = []
            #for i in xrange(len(self.zchi2_cap.shape)-1):
            #    ssploc.append(n.where(self.zchi2_cap[self.fiber-self.fos] == n.min(self.zchi2_cap[self.fiber-self.fos]))[i][0])
            #pmat = n.zeros((4663,5))
            #pmat[:,0] = self.startemps1[tuple(ssploc[:-1])][ssploc[-1]+864:ssploc[-1]+4663+864] #653 for spEigenStar, 864 for CAP grids, 667 for Charlie's SSPs
            #pmat[:,1:] = n.transpose(poly_array(4, 4663))
            #ninv = n.diag(self.ivar[self.fiber])
            #coeffvec = n.linalg.solve(n.dot(n.transpose(pmat),n.dot(ninv,pmat)),n.dot(n.transpose(pmat),n.dot(ninv,self.specs[self.fiber])))
            #model3 = n.dot(pmat,coeffvec)
            #a.plot(self.wave, model3, color='black')

        a.set_xlabel('Wavelength (Angstroms)')
        a.set_ylabel('Flux in some units')
        canvas = FigureCanvasTkAgg(f, master=self.master)
        canvas.show()
        canvas.get_tk_widget().grid(row=5, column=5)

    def next_fiber(self):
        self.fiber += 1
        self.e2.delete(0, END)
        self.e2.insert(0, str(self.fiber))
        self.do_plot()

    def prev_fiber(self):
        self.fiber -= 1
        self.e2.delete(0, END)
        self.e2.insert(0, str(self.fiber))
        self.do_plot()

#root = Tk()
#app = Plot_fit(root)
#root.mainloop()
