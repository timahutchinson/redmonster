import numpy as n
#import matplotlib as m
#m.interactive(True)
#from matplotlib import pyplot as p
from astropy.io import fits
#from astropy import wcs
#import os
#import pixelsplines as pxs

# Script to make test file for "read_ndArch" file-reader routine:

# Set the plan:
#naxis1 = 3000 # log-wavelength
naxis1 = 8 # log-wavelength
naxis2 = 5    # irregular
naxis3 = 10   # PS-labeled
#naxis4 = 20   # regular
naxis4 = 6   # regular
naxis5 = 2    # N-labeled
naxis6 = 3    # Default index

# Generate the simulation data and put it into an HDU object:
image = n.zeros((naxis6, naxis5, naxis4, naxis3, naxis2, naxis1), dtype='float32')
image[...] = n.random.uniform(size=(naxis3, naxis2, naxis1))
hdu = fits.PrimaryHDU(image)

# 0. Indicate the units:
hdu.header.set('BUNIT', value='10^(-17) erg/cm^2/s/Angstrom', comment='Data unit')

# 1. Set the wavelength basis information:
hdu.header.set('CNAME1', value='loglam', comment='Axis 1 name')
hdu.header.set('CUNIT1', value='log10(Angstroms)', comment='Axis 1 unit')
hdu.header.set('CRPIX1', value=1., comment='Axis 1 reference pixel')
hdu.header.set('CRVAL1', value=3.5, comment='Axis 1 reference value')
hdu.header.set('CDELT1', value=0.0001, comment='Axis 1 increment')

# 2. Set the irregular grid parameters:
hdu.header.set('CNAME2', value='Metallicity', comment='Axis 2 name')
hdu.header.set('CUNIT2', value='Z', comment='Axis 2 unit')
this_axis = str(2)
metvals = [0.0006, 0.00012, 0.0049, 0.0190, 0.03]
for i in xrange(len(metvals)):
    hdu.header.set('PV' + this_axis + '_' + str(i+1), metvals[i],
                   comment='Axis 2 value at pixel ' + str(i+1))

# 3. Set the PS-labeled grid parameters:
hdu.header.set('CNAME3', value='Spectral Type', comment='Axis 3 name')
hdu.header.set('CUNIT3', value='OBAFGKMLT', comment='Axis 3 unit')
this_axis = str(3)
labels = ['O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'T', 'Y']
for i in xrange(len(labels)):
    hdu.header.set('PS' + this_axis + '_' + str(i+1), labels[i],
                   comment='Axis 3 label at pixel ' + str(i+1))

# 4. Set the regular parameter grid:
hdu.header.set('CNAME4', value='Velocity Dispersion', comment='Axis 4 name')
hdu.header.set('CUNIT4', value='km/s', comment='Axis 4 unit')
hdu.header.set('CRPIX4', value=1., comment='Axis 4 reference pixel')
hdu.header.set('CRVAL4', value=100, comment='Axis 4 reference value')
hdu.header.set('CDELT4', value=25.0, comment='Axis 4 increment')

# 5. Set the name-labeled axis info:
hdu.header.set('CNAME5', value='Familiar Name', comment='Axis 5 name')
this_axis = str(5)
names = ['Foo', 'Bar']
for i in xrange(len(names)):
    hdu.header.set('N' + this_axis + '_' + str(i+1), names[i],
                   comment='Axis 5 name at pixel ' + str(i+1))

# 6. For the final axis, we default to an index, and include nothing.

# Write it out:
ofile = 'ndArch-TEST-v00.fits'
hdu.writeto(ofile, clobber=True)

