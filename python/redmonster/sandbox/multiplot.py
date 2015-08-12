import numpy as n
from astropy.io import fits
from matplotlib import pyplot as p
from astropy.convolution import convolve, Box1DKernel
from os.path import join
from os import environ

# Open spPlate file
path = '/Users/timhutchinson/work/boss/run2d/7848/spPlate-7848-56959.fits'
spp = fits.open(path)
sigma = 1. / n.sqrt( spp[1].data )
wave = 10**( spp[0].header['COEFF0'] + n.arange(spp[0].header['NAXIS1'])*spp[0].header['COEFF1'] )

# Open redmonster file
path = '/Users/timhutchinson/work/boss/run2d/7848/run1d/redmonster-7848-56959.fits'
redm = fits.open(path)

f = p.figure()
# z~0.6 passive galaxy (fiber 100)
ax1 = f.add_subplot(321)
p.plot(wave, sigma[100], color='red')
p.plot(wave, convolve(spp[0].data[100]+.5, Box1DKernel(5)), color='black')
p.plot(wave, redm[2].data[100]+.5, color='cyan')
p.axis([wave[0],10000,0,2])
p.xlabel('Observed Wavelength ($\AA$)')
p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
p.text(wave[0]+350, 2-.17, '(a) 7848-56959-101', fontsize=12)
ax1.set_yticks([0.0,0.5,1.0,1.5,2.0])

# z~1.0 passive galaxy (fiber 106)
ax2 = f.add_subplot(322)
p.plot(wave, sigma[106], color='red')
p.plot(wave, convolve(spp[0].data[106]+1, Box1DKernel(5)), color='black')
p.plot(wave, redm[2].data[106]+1, color='cyan')
p.axis([wave[0],10000,0,2.5])
p.xlabel('Observed Wavelength ($\AA$)')
p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
p.text(wave[0]+350, 2.5-.22, '(b) 7848-56959-107', fontsize=12)
ax2.set_yticks([0.0,0.5,1.0,1.5,2.0,2.5])

# EL galaxy (fiber 96)
ax3 = f.add_subplot(323)
p.plot(wave, sigma[96], color='red')
p.plot(wave, convolve(spp[0].data[96]+.5, Box1DKernel(5)), color='black')
p.plot(wave, redm[2].data[96]+.5, color='cyan')
p.axis([wave[0],10000,-0.1,5.5])
p.xlabel('Observed Wavelength ($\AA$)')
p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
p.text(wave[0]+350, 5.5-.5, '(c) 7848-56959-097', fontsize=12)
ax3.set_yticks([0.5,1.5,2.5,3.5,4.5,5.5])

# Star (fiber 64)
ax4 = f.add_subplot(324)
p.plot(wave, sigma[64], color='red')
p.plot(wave, convolve(spp[0].data[64]+.5, Box1DKernel(5)), color='black')
p.plot(wave, redm[2].data[64]+.5, color='cyan')
p.axis([wave[0],10000,0,2.5])
p.xlabel('Observed Wavelength ($\AA$)')
p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
p.text(wave[0]+350, 2.5-.22, '(d) 7848-56959-065', fontsize=12)

# z~1.0 QSO (fiber 63)
ax5 = f.add_subplot(325)
p.plot(wave, sigma[63], color='red')
p.plot(wave, convolve(spp[0].data[63]+.5, Box1DKernel(5)), color='black')
p.plot(wave, redm[2].data[63]+.5, color='cyan')
p.axis([wave[0],10000,-0.1,6.5])
p.xlabel('Observed Wavelength ($\AA$)')
p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
p.text(wave[0]+350, 6.5-.58, '(e) 7848-56959-064', fontsize=12)
ax5.set_yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5])

# z~3.0 QSO (fiber 75)
ax6 = f.add_subplot(326)
p.plot(wave, sigma[75], color='red')
p.plot(wave, convolve(spp[0].data[75]+.5, Box1DKernel(5)), color='black')
p.plot(wave, redm[2].data[75]+.5, color='cyan')
p.axis([wave[0],10000,-0.1,5.5])
p.xlabel('Observed Wavelength ($\AA$)')
p.ylabel('$f_\lambda$ $10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$')
p.text(wave[0]+350, 5.5-.48, '(f) 7848-56959-076', fontsize=12)
ax6.set_yticks([0.5,1.5,2.5,3.5,4.5,5.5])








p.subplots_adjust(wspace = .3, hspace = .3)
#p.savefig('/Users/timhutchinson/Desktop/6plot.pdf')


# Addtional code for manually scrolling through to find good example plots
'''
from os import environ; from os.path import join; from astropy.io import fits; import numpy as n;import matplotlib;from matplotlib import pyplot as p;p.interactive(True); from astropy.convolution import convolve, Box1DKernel
plate = 7851
mjd = 56932
fiber = 487
p.plot(convolve(fits.open(join(environ['BOSS_SPECTRO_REDUX'],environ['RUN2D'],'%s'%plate,'spPlate-%s-%s.fits'%(plate,mjd)))[0].data[fiber],Box1DKernel(5)), 'k'); p.plot(fits.open(join(environ['REDMONSTER_SPECTRO_REDUX'],environ['RUN2D'],'%s'%plate,environ['RUN1D'],'redmonster-%s-%s-%03d.fits'%(plate,mjd,fiber)))[2].data[0,0],color='cyan')
'''



