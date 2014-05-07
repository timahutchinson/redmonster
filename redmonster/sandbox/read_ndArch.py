# read_ndArch.py
#
# Code for reading ndArch files
#
# bolton@utah@iac 2014mayo
#

import numpy as n
from astropy.io import fits

def read_ndArch(fname):
    """
    Read in an ndArch archetype file, parsing parameter baselines
    """
    # fname = 'ndArch-TEST-v00.fits'
    # Parse class and version from the filename:
    fn_ROOT = fname[:fname.rfind('.fits')]
    fn_CLASS = fn_ROOT.split('-')[1]
    fn_VERSION = fn_ROOT.split('-')[-1]
    # Get the data and header:
    data = fits.getdata(fname).copy()
    header = fits.getheader(fname)
    # Identify how many parameters:
    npars = len(data.shape) - 1
    emptyparlist = ['']
    # Initialize output info dictionary:
    infodict = {'class': fn_CLASS,
                'version': fn_VERSION,
                'coeff0': header['CRVAL1'],
                'coeff1': header['CDELT1'],
                'nwave': header['NAXIS1'],
                'fluxunit': '',
                'par_names': ['']*npars,
                'par_units': ['']*npars,
                'par_axistype': ['index']*npars}
    if ('BUNIT' in header): infodict['fluxunit'] = header['BUNIT']
    # Initialize list of baselines with index defaults:
    baselines = [n.arange(this_size)+1 for this_size in data.shape[:-1]]
    # Loop over parameters and construct baselines:
    for ipar in xrange(npars):
        # Translate Python axis index integer to FITS axis index string:
        ax = str(npars + 1 - ipar)
        # Populate name & units for this axis, if available:
        if ('CNAME'+ax in header): infodict['par_names'][ipar] = header['CNAME'+ax]
        if ('CUNIT'+ax in header): infodict['par_units'][ipar] = header['CUNIT'+ax]
        # The axis condition tests -- maybe inefficient to always compute
        # all of these, but makes for nicer code:
        is_regular = ('CRPIX'+ax in header) and \
                     ('CRVAL'+ax in header) and \
                     ('CDELT'+ax in header)
        pv_base = ['PV'+ax+'_'+str(j+1) for j in xrange(data.shape[ipar])]
        pv_test = n.asarray([this_pv in header for this_pv in pv_base])
        is_irregular = pv_test.prod() > 0
        ps_base = ['PS'+ax+'_'+str(j+1) for j in xrange(data.shape[ipar])]
        ps_test = n.asarray([this_ps in header for this_ps in ps_base])
        is_labeled = ps_test.prod() > 0
        n_base = ['N'+ax+'_'+str(j+1) for j in xrange(data.shape[ipar])]
        n_test = n.asarray([this_n in header for this_n in n_base])
        is_named = n_test.prod() > 0
        if is_regular:
            baselines[ipar] = (n.arange(data.shape[ipar]) + 1 - header['CRPIX'+ax]) \
                            * header['CDELT'+ax] + header['CRVAL'+ax]
            infodict['par_axistype'][ipar] = 'regular'
        if is_irregular:
            baselines[ipar] = n.asarray([header[this_pv] for this_pv in pv_base])
            infodict['par_axistype'][ipar] = 'irregular'
        if is_labeled:
            baselines[ipar] = n.asarray([header[this_ps] for this_ps in ps_base])
            infodict['par_axistype'][ipar] = 'labeled'
        if is_named:
            baselines[ipar] = n.asarray([header[this_n] for this_n in n_base])
            infodict['par_axistype'][ipar] = 'named'
    return data, baselines, infodict
