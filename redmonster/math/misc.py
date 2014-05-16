import numpy as n
from scipy import special as spc

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

def gaussflux(pixbound, cen, sig):
    """
    For monotonically increasing pixel boundaries specified by 'pixbound'
    in some abscissa units, consider a Gaussian with unit integrated
    amplitude that expresses a density per those same abscissa units,
    centered on 'cen' and with sigma parameter 'sig', and return
    the average value of that Gaussian between the boundaries
    (i.e., its pixel-averaged value for possibly non-uniform pixels.)
    (bolton@utah@iac 2014mayo)
    """
    # Calculate the pixel widths and test for monotonicity:
    pixdiff = pixbound[1:] - pixbound[:-1]
    if (pixdiff.min <= 0):
        print 'pixbound must be monotonically increasing!'
        return 0
    # Compute scaled argument for error function:
    argscale = (pixbound - cen) / (n.sqrt(2.) * sig)
    # Compute and return the argument:
    return 0.5 * (spc.erf(argscale[1:]) - spc.erf(argscale[:-1])) / pixdiff



















