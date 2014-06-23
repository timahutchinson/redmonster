import numpy as n
from scipy import special as spc
from scipy import sparse

# Function to find where S/N is unreasonably large or where flux is unphysically negative
def flux_check(flux, ivars):
    for i in range(flux.shape[0]):
        ct = n.where(abs(flux[i]) * n.sqrt(ivars[i]) > 200.)[0].shape[0]
        # CHANGE NEXT LINE SO IT ADDS TO LOG FILE RATHER THAN PRINTS
        if ct > 0: print 'WARNING: Fiber #%s has %s pixels with S/N > 200' % (i+1,ct)
        badpix = n.where(flux[i] * n.sqrt(ivars[i]) < -10.)[0]
        if len(badpix) > 0:
            # ALSO CHANGE TO ADD TO LOG
            print 'WARNING: Fiber #%s has %s pixels with Flux < -10*Noise' % (i+1,len(badpix))
            ivars[i] = mask_pixels(badpix, ivars[i])
    return ivars

# Mask unphysically negative pixels + neighboring two pixels in both directions
def mask_pixels(badpix, ivars):
    for j in badpix:
        ivars[j-2:j+3] = 0
    return ivars

# Function to transform pixel centers to pixel boundaries
def cen2bound(pixelcen):
    pixbound = 0.5 * (pixelcen[1:] + pixelcen[:-1])
    pixbound = n.append( n.append( 2.*pixbound[0]-pixbound[1], pixbound ), 2.*pixbound[-1]-pixbound[-2] )
    return pixbound

# Function to transform from pixel boundaries to pixel centers
def bound2cen(pixbound):
    pixbound = .5 * (pixbound[:-1] + pixbound[1:])
    return pixbound

# Create poly array to add polynomial terms in fitting
def poly_array(npoly, npix):
    arr = n.zeros(shape=(npoly,npix))
    xvec = n.arange(npix) / float(npix)
    for i in range(npoly): arr[i] = xvec**i
    return arr

def two_pad(npix):
    i = 0
    res = 0
    while (res == 0):
        i += 1
        res = (2**i) // int(abs(npix))
    return (2**i)

def multipoly_fit(ind, dep, order=2):
    ndata = n.prod(dep.shape)
    ndim = ind.shape[0]
    A = n.zeros((ndata, ndim*(order+1)))
    b = n.zeros(ndata)
    for i in xrange(ndata):
        for j in xrange(ndim*(order+1)):
            A[i,j] = 0

def quadfit_2d(ind, dep): # Requires (2x3) matrix ind where ind[0] is 3 coordinates along one dimension and ind[1] is three coordinates along other, dep is (3x3) matrix of datapoints where dep[i,j] corresponds to ind[0,i], ind[1,j], returns coeffs for ax^2+by^2+cxy+dx+ey+f=z
    x = ind[0]
    y = ind[1]
    A = n.zeros((9,6))
    b = n.reshape(dep,(9,1))
    for k in range(9):
        for i in xrange(3):
            for j in xrange(3):
                A[k] = n.array([ x[i]**2, y[j]**2, x[i]*y[j], x[i], y[j], 1 ])
    f = n.dot(n.linalg.pinv(A),b)
    return f

def quadfit(ind, dep):
    A = n.zeros((3,3))
    for i in xrange(3):
        A[i] = n.array([ ind[i]**2, ind[i], 1 ])
    f = n.linalg.solve(A,dep)
    return f

def gaussflux(pixbound, cen, sig, h_order=0):
    """
    For monotonically increasing pixel boundaries specified by 'pixbound'
    in some abscissa units, consider a Gaussian with unit integrated
    amplitude that expresses a density per those same abscissa units,
    centered on 'cen' and with sigma parameter 'sig', and return
    the average value of that Gaussian between the boundaries
    (i.e., its pixel-averaged value for possibly non-uniform pixels.)
    
    Can specify a Gauss-Hermite order (a la scipy.special.hermitenorm)
    via the 'h_order' argument, which defaults to zero.
    
    (bolton@utah@iac 2014mayo)
    (Added Hermite orders: bolton@utah@iac 2014junio)
    """
    # Calculate the pixel widths and test for monotonicity:
    pixdiff = pixbound[1:] - pixbound[:-1]
    if (pixdiff.min <= 0):
        print 'pixbound must be monotonically increasing!'
        return 0
    # Make sure scalar arguments are scalars:
    if (n.asarray(cen).size != 1):
        print 'cen argument must be scalar!'
        return 0
    if (n.asarray(sig).size != 1):
        print 'sig argument must be scalar!'
        return 0
    # Compute and return:
    if h_order > 0:
        u = (pixbound - cen) / sig
        int_term = - spc.hermitenorm(h_order-1)(u) * n.exp(-0.5 * u**2) / n.sqrt(2. * n.pi)
    else:
        int_term = 0.5 * spc.erf((pixbound - cen) / (n.sqrt(2.) * sig))
    return (int_term[1:] - int_term[:-1]) / pixdiff

def gaussbasis(pixbound, cen, sig, h_order=0, nsigma=6.0):
    """
    Function to generate a sparse matrix to broadcast integrated
    Gaussian amplitudes with sigma parameter 'sig' and centered on 'cen'
    into pixel-averaged values over the pixel baseline specified by 'pixbound'.
    Uses scipy sparse matrix class to implement.
    Can also accept Gauss-Hermite order parameter in 'h_order'.
    Integrates out to at least +/- 'nsigma' times the sig parameter,
    with nsigma=6.0 by default.
    (bolton@utah@iac 2014junio)
    """
    # How many input Gaussians:
    ngauss = len(cen)
    # How many output pixels?
    npix = len(pixbound) - 1
    # Make sure there are same number of 'cen' and 'sig' values:
    if (len(sig) != ngauss):
        print 'Lengths of cen and sig must match!'
        return 0
    # Check for monotonicity:
    dpix = pixbound[1:] - pixbound[:-1]
    if (dpix.min() <= 0.):
        print 'Pixel boundaries not monotonically increasing!'
        return 0
    # Work out indices of the bins into which +/- nsigma values fall:
    bin_lo = n.digitize(cen - nsigma * sig, pixbound) - 1
    bin_hi = n.digitize(cen + nsigma * sig, pixbound) - 1
    # Limit to valid range:
    bin_lo = n.where((bin_lo >= 0), bin_lo, 0)
    bin_hi = n.where((bin_hi < npix), bin_hi, npix-1)
    # Initialize matrix:
    gbasis = sparse.lil_matrix((ngauss,npix))
    # Loop over Gaussians, compute, and return:
    for i in xrange(ngauss):
        if (bin_hi[i] >= bin_lo[i]):
            gbasis[i,bin_lo[i]:bin_hi[i]+1] = gaussflux(pixbound[bin_lo[i]:bin_hi[i]+2],
                                                        cen[i], sig[i], h_order=h_order).reshape((1,-1))
    return gbasis.tocsr().T

def gaussproj(pixbound_in, sigma_in, pixbound_out, h_order=0, nsigma=6.0):
    """
    Function to generate a matrix that projects from one pixelized
    spectral baseline to another pixelized spectral baseline, also
    implementing wavelength-dependent Gaussian (or Gauss-Hermite)
    convolution along the way.

    Arguments (mandatory):
      pixbound_in: monotonically increasing array of boundaries of
        input pixels (dimension npix_in + 1).
      sigma_in: Gaussian sigma associated with each of the input pixels.
        Must be of dimension npix_in, expressed in same units as pixbound_in.
      pixbound_out: monotonically increasing array of boundaries of
        output/target pixels (dimension npix_out + 1), expressed in same
        units as pixbound_in.

    Arguments (optional):
      h_order: Hermite polynomial to multiply into convolving Gaussian.
        Default is h_order=0.  Uses scipy.special.hermitenorm(), not the
        van der Marel & Franx 1993 Gauss-Hermite convention!
      nsigma: number of sigmas +/- to integrate out for Gaussian broadening
        kernel.  Default is h_order=6.  Increasing will increase accuracy,
        but at the expense of reduced sparseness in the matrix.

    Returns:
      A scipy sparse matrix of dimension (npix_out X npix_in) that,
      when multiplied into an input-space spectrum, broadcasts it
      to an output-space spectrum with the desired resampling and
      broadening.  Input spectra should be in 'flux-per-unit-x',
      where 'x' is the baseline unit of the mandatory arguments.

    Note: this function is especially designed for the case of
    resampling-with-broadening.  If you just want to resample without
    broadening, you probably want to use pixelsplines instead.

    Written: bolton@utah@iac 2014junio
    """
    # How many input and output pixels?
    npix_in = len(pixbound_in) - 1
    npix_out = len(pixbound_out) - 1
    # Check for monotonicity:
    dpix_in = pixbound_in[1:] - pixbound_in[:-1]
    dpix_out = pixbound_out[1:] - pixbound_out[:-1]
    if (dpix_in.min() <= 0.):
        print 'Input pixel boundaries not monotonically increasing!'
        return 0
    if (dpix_out.min() <= 0.):
        print 'Output pixel boundaries not monotonically increasing!'
        return 0
    # Make sure the sigma values are matched to the input pixels:
    if (len(sigma_in) != npix_in):
        print 'Must match length of sigma_in to input pixel space!'
        return 0
    # Generate a diagonal matrix that turns the flux density of the
    # input spectrum into the integrated flux of the input pixels:
    pix_flux_mat = sparse.dia_matrix((dpix_in, 0), shape=(npix_in,npix_in))
    # Generate a sparse matrix that takes unit flux in the input pixels
    # to sigma-distributed flux density in the output pixels:
    center_in = 0.5 * (pixbound_in[1:] + pixbound_in[:-1])
    proj_mat = gaussbasis(pixbound_out, center_in, sigma_in,
                          h_order=h_order, nsigma=nsigma)
    # Compose those two matrices and return the result:
    return proj_mat * pix_flux_mat

