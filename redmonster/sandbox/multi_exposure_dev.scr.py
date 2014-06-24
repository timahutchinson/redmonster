#
# multi_exposure_dev.scr.py
#
# script for developping methods that handle multiple
# exposures of a single fiber.
#
import numpy as n
import matplotlib as m
m.interactive(True)
from matplotlib import pyplot as p
from astropy.io import fits
import os
from redmonster.math import misc
#from scipy import signal as sig

#import pixelsplines as pxs

# Set the following:
# export BOSS_SPECTRO_REDUX='/data/BOSS/redux/dr10mini'
# export RUN2D=v5_5_12
# export RUN1D=v5_5_12

plate = 3686
mjd = 55268
fiberid = 265

# Code to get the individual exposures for this plate:
spf = os.getenv('BOSS_SPECTRO_REDUX') + '/' + os.getenv('RUN2D') + '/' + \
      str(plate).strip() + '/spPlate-' + str(plate).strip() + '-' + \
      str(mjd).strip() + '.fits'
hdr = fits.getheader(spf)

exp_keys = [this_key for this_key in hdr.keys() if this_key[:5] == 'EXPID']

exp_ids = [hdr[this_key][:11] for this_key in exp_keys]

# Build the individual file exposure names:
path_to_spectra = os.getenv('BOSS_SPECTRO_REDUX') + '/' + \
                  os.getenv('RUN2D') + '/' + str(plate).strip() + '/'
spCFrame_list = [path_to_spectra + 'spCFrame-' + this_id + '.fits' for this_id in exp_ids]

# Get the data:
data_list = [fits.getdata(this_file) for this_file in spCFrame_list]
invvar_list = [fits.getdata(this_file, 1) for this_file in spCFrame_list]
loglam_list = [fits.getdata(this_file, 3) for this_file in spCFrame_list]
sigma_list = [fits.getdata(this_file, 4) for this_file in spCFrame_list]
plug_list = [fits.getdata(this_file, 5) for this_file in spCFrame_list]

# Find the indices of the fiberid of interest:
i_exp = []
j_row = []
for i in xrange(len(plug_list)):
    wh_fib = n.where(plug_list[i].FIBERID == fiberid)[0]
    for j in xrange(len(wh_fib)):
        i_exp.append(i)
        j_row.append(wh_fib[j])

n_spec = len(i_exp)

# Pull out the bits we actually want:
data_list_one = [data_list[i_exp[k]][j_row[k]] for k in xrange(n_spec)]
invvar_list_one = [invvar_list[i_exp[k]][j_row[k]] for k in xrange(n_spec)]
loglam_list_one = [loglam_list[i_exp[k]][j_row[k]] for k in xrange(n_spec)]
sigma_list_one = [sigma_list[i_exp[k]][j_row[k]] for k in xrange(n_spec)]

# Derived vectors:
wave_list_one = [10.**loglam_this for loglam_this in loglam_list_one]
logbound_list_one = [misc.cen2bound(loglam_this) for loglam_this in loglam_list_one]
wavebound_list_one = [10.**logbound_this for logbound_this in logbound_list_one]
sigwave_list_one = [10.**(-4) * n.log(10.) * wave_list_one[i] * sigma_list_one[i]
                    for i in xrange(n_spec)]


# That basically does it for the data.
# Next thing is to get the models into a compatible form!

# Although perhaps we should encapsulate this spCFrame
# handling code into a slicker interface...


# Plot this and see if it makes sense:
hold_val = n_spec * [True]
hold_val[0] = False

for i in xrange(n_spec):
    p.plot(wave_list_one[i], data_list_one[i], hold=hold_val[i], color='k')



for i in xrange(n_spec):
    p.plot(wave_list_one[i], sigwave_list_one[i], hold=hold_val[i], color='k')





wh_fib = [n.where(this_plug.FIBERID == fiberid)[0] for this_plug in plug_list]

[
