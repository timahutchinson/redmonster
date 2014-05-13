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

class Plot_fit:

    def __init__(self, master, zchi2_ssp, zchi2_star):
        self.master = master
        self.plate = None
        self.zchi2_ssp = zchi2_ssp
        self.zchi2_star = zchi2_star
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
        self.startemps = fits.open(join(environ['REDMONSTER_DIR'],'templates','ndArch-spEigenStar-55734.fits'))[0].data
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
                ssploc.append(n.where(self.zchi2_ssp[self.fiber] == n.min(self.zchi2_ssp[self.fiber]))[i][0])
            pmat = n.zeros((4663,4))
            pmat[:,0] = self.ssptemps[tuple(ssploc[:-1])][ssploc[-1]+667:ssploc[-1]+4663+667]
            pmat[:,1:] = n.transpose(poly_array(3, 4663))
            ninv = n.diag(self.ivar[self.fiber])
            coeffvec = n.linalg.solve(n.dot(n.transpose(pmat),n.dot(ninv,pmat)),n.dot(n.transpose(pmat),n.dot(ninv,self.specs[self.fiber])))
            model = n.dot(pmat,coeffvec)
            a.plot(self.wave, model, color='black')
            # Overplot best-fit star
            ssploc = []
            for i in xrange(len(self.zchi2_star.shape)-1):
                ssploc.append(n.where(self.zchi2_star[self.fiber] == n.min(self.zchi2_star[self.fiber]))[i][0])
            pmat = n.zeros((4663,5))
            pmat[:,0] = self.startemps[tuple(ssploc[:-1])][ssploc[-1]+653:ssploc[-1]+4663+653]
            pmat[:,1:] = n.transpose(poly_array(4, 4663))
            ninv = n.diag(self.ivar[self.fiber])
            coeffvec = n.linalg.solve(n.dot(n.transpose(pmat),n.dot(ninv,pmat)),n.dot(n.transpose(pmat),n.dot(ninv,self.specs[self.fiber])))
            model2 = n.dot(pmat,coeffvec)
            a.plot(self.wave, model2, color='blue')
        a.set_xlabel('Wavelength (Angstroms)')
        a.set_ylabel('Flux in some units')
        a.set_ylim([n.min(model)-5, n.max(model)+5])
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
