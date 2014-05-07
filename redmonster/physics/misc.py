import numpy as n

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

# Mask unphysically negative pixels + neighboring two pixels in both directions
def mask_pixels(badpix, ivars):
    for j in badpix:
        ivars[j-2:j+2] = 0
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