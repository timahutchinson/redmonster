#
# multifit.py
#
# Module of routines for enabling fitting of multiple
# spectra at one time
#
# bolton@utah@iac 2014junio
#

import numpy as n
from redmonster.math import misc

def multi_projector(wavebound_list, sigma_list, coeff0, coeff1):
    """
    Function to take a list of spectro wavelength baselines and
    associated instrumental dispersion parameters and return
    a list of projection matrices from a common uniform
    log-lambda grid into the frame of the individual exposures,
    including broadening by instrumental dispersion.

    Arguments:
      wavebound_list: List of 1D vectors each containing the
        pixel boundaries of the individual spectra.  Each
        vector is taken to be of length npix_j + 1, where
        npix_j is the number of pixels in the j'th spectrum
        Must be monotonically increasing or something will crash.
      sigma_list: List of instrumental Gaussian sigma dispersion
        parameter vectors.  Each vector to be of length npix_j
      coeff0: log10-angstroms (vacuum, rest frame) of the zeroth
        pixel of the constant log10-lambda grid from which the
        projection matrices should broadcast to the individual
        spectra.
      coeff1: delta log10-angstroms per pixel of the constant
        log10-lambda grid from which the projection matrices
        should broadcast to the individual spectra.

    Returns:
      (matrix_list, idx_list, nsamp_list)
    where
      matrix_list is a list of (scipy.sparse) projection
        matrices that broadcast from the model baseline to
        the individual spectrum baselines, and
      idx_list is a list of the indices of the zeroth pixel of
        the model-grid matrix dimensions within the baseline
        specified by the input coeff0/coeff1 values, and
      nsamp_list is a list of the number of model-space
        sample-pixels encompassed by each matrix.

    To make that more explicit,
      matrix_list[j] will have formal dimensions
        npix_j X nsamp_j
      and the dimesion of length nsamp_j aligns with
        loglam[idx_list[j]:idx_list[j]+nsamp_list[j]]
      where
        loglam = coeff0 + coeff1 + n.arange(...)
      is the log10-lambda baseline of the model grid.

      This all ensures that we don't waste time on any
      dimensional coverage that we don't need, and that
      we can slide things in redshift by pixel-shifting
      the projection matrices within the models.

    Written: bolton@utah@iac 2014junio
    """
    # Number of spectra:
    nspec = len(wavebound_list)
    # Number of pixels in each spectrum:
    npix_list = [len(this_sigma) for this_sigma in sigma_list]
    # Wavelengths of a 6-sigma buffer at the high and low ends:
    wavelim_lo = [wavebound_list[k][0] - 10. * sigma_list[k][0] for k in xrange(nspec)]
    wavelim_hi = [wavebound_list[k][-1] + 10. * sigma_list[k][-1] for k in xrange(nspec)]
    # Translate these into indices within the nominal full model baseline
    idx_list = [int(round((n.log10(this_wave) - coeff0) / coeff1)) for this_wave in wavelim_lo]
    idx_hi = [int(round((n.log10(this_wave) - coeff0) / coeff1)) for this_wave in wavelim_hi]
    nsamp_list = [idx_hi[k] - idx_list[k] + 1 for k in xrange(nspec)]
    # Compute the nominal wavelength arrays for the spectra:
    wave_list = [0.5 * (this_bound[1:] + this_bound[:-1]) for this_bound in wavebound_list]
    # Compute the various model-space wavelength baselines that we need:
    modloglam_list = [coeff0 + coeff1 * (n.arange(nsamp_list[k]) + idx_list[k]) for k in xrange(nspec)]
    modlogbound_list = [misc.cen2bound(this_loglam) for this_loglam in modloglam_list]
    modwave_list = [10.**this_loglam for this_loglam in modloglam_list]
    modwavebound_list = [10.**this_logbound for this_logbound in modlogbound_list]
    # Interpolate the spectrum-frame sigmas onto the model-frame grids:
    modsigma_list = [n.interp(modwave_list[k], wave_list[k], sigma_list[k]) for k in xrange(nspec)]
    # Compute the projection matrices:
    matrix_list = [misc.gaussproj(modwavebound_list[k], modsigma_list[k],
                                  wavebound_list[k]) for k in xrange(nspec)]
    # Return results:
    return matrix_list, idx_list, nsamp_list
