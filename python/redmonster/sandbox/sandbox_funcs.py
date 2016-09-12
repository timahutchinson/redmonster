# sandbox_funcs.py
#
# place to develop functions in sandbox

from scipy import special as spc
import numpy as n


def gaussflux(pixbound, cen, sig):
    """
    For monotonically increasing pixel boundaries specified by 'pixbound'
    in some abscissa units, consider a Gaussian with unit integrated
    amplitude that expresses a density per those same abscissa units,
    centered on 'cen' and with sigma parameter 'sig', and return
    the average value of that Gaussian between the boundaries
    (i.e., its pixel-averaged value for possibly non-uniform pixels.)
    """
    # Calculate the pixel widths and test for monotonicity:
    pixdiff = pixbound[1:] - pixbound[:-1]
    if (pixdiff.min <= 0):
        print('pixbound must be monotonically increasing!')
        return 0
    # Compute scaled argument for error function:
    argscale = (pixbound - cen) / (n.sqrt(2.) * sig)
    # Compute and return the argument:
    return 0.5 * (spc.erf(argscale[1:]) - spc.erf(argscale[:-1])) / pixdiff

