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
        parameter vectors.  Each vector to be of length npix_j.
        Should be in same units as wavebound_list.
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

class MultiProjector:
    """
    Class to take a list of spectro wavelength baselines and
    associated instrumental dispersion parameters and return
    an object that implements projection from uniform log10-lambda
    baseline model grid.

    This is basically an object-wrapper to the function
      multi_projector
    See the documentation for that function for more information
    on the internals of this class.
    
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

    Written: bolton@utah@iac 2014junio
    """
    def __init__(self, wavebound_list, sigma_list, coeff0, coeff1):
        """
        Constructor for the MultiProjector object.

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
        """
        self.nspec = len(wavebound_list)
        self.npix_list = [len(this_sigma) for this_sigma in sigma_list]
        self.coeff0 = coeff0
        self.coeff1 = coeff1
        self.matrix_list, self.idx_list, self.nsamp_list = \
            multi_projector(wavebound_list, sigma_list, coeff0, coeff1)
    def project_model_grid(self, model_grid, pixlag=0, coeff0=None):
        """
        Function to project a grid of constant-log10-lambda models onto
        the individual frames of multiple spectra, with pixel-redshift
        
        Arguments:
          model_grid: ndarray of models gridded in constant log10-Angstroms.
            Wavelength dimension must be final dimension.  Any number of
            leading dimensions is allowed.
          pixlag: pixel-redshift to apply when placing projection matrices
            within model loglam grids.  By convention, pixlag > 0 is redshift
            and pixlag < 0 is blieshift.
          coeff0: Override value if the model grid has a zero-pixel
            log10-Angstrom value other than that which is in self.coeff0.
            This will be converted to pixels and rounded to an integer value.
        """
        # Do we need to offset for a different coeff0?
        if coeff0 is None:
            ishift = 0
        else:
            # If argument coeff0 is greater than self.coeff0, then
            # we need to index our matrices into lower-numbered
            # indices within the model baseline.  If this corresponds
            # to a positive 'ishift' value as it does here, then it
            # corresponds to subtraction when building the slices
            # in the code immediately below.
            ishift = int(round((coeff0 - self.coeff0) / self.coeff1))
        # Build a list of slices within the model grid:
        slice_list = [slice(self.idx_list[k]-pixlag-ishift,
                            self.idx_list[k]+self.nsamp_list[k]-pixlag-ishift)
                      for k in xrange(self.nspec)]
        # How many pixels in the model grids?
        npix_model = model_grid.shape[-1]
        # Dimensionality of the model-grid space:
        dimshape_model = model_grid.shape[:-1]
        # Total number of models for looping:
        nmodels = model_grid.size // npix_model
        # Make a flattened view of the model grids for looping:
        model_flatgrid = model_grid.reshape((nmodels, npix_model))
        # Initialize output list, in flattened form:
        outgrid_list = [n.zeros((nmodels, this_npix), dtype=float)
                        for this_npix in self.npix_list]
        # Now loop over exposures and models:
        for j_spec in xrange(self.nspec):
            for i_mod in xrange(nmodels):
                outgrid_list[j_spec][i_mod] = self.matrix_list[j_spec] \
                  * model_flatgrid[i_mod,slice_list[j_spec]]
            # Resize the output grid to match the input model-space dimensions:
            outgrid_list[j_spec].resize(dimshape_model + (self.npix_list[j_spec],))
        return outgrid_list
    def single_poly_nonneg(self, npoly):
        """
        Method to generate a single global (model-space) observed-frame
        non-negative polynomial basis and project it through the projection matrices
        into the frames of the individual spectra.
        """
        idx_lo = min(self.idx_list)
        idx_hi = max(n.asarray(self.idx_list) + n.asarray(self.nsamp_list))
        npix_poly = idx_hi - idx_lo
        poly_base = n.arange(npix_poly) / float(npix_poly-1)
        poly_grid = n.zeros((2*int(round(npoly)), npix_poly), dtype=float)
        for ipoly in xrange(int(round(npoly))):
            poly_grid[2*ipoly] = poly_base**ipoly
            poly_grid[2*ipoly+1] = - poly_base**ipoly
        return self.project_model_grid(poly_grid, pixlag=idx_lo)
