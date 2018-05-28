#
# airtovac.py
#
"""

Module to enable inter-conversion between air and vacuum wavelengths.
Uses formula from Ciddor 1996, Applied Optics, 35, 1566.
Adapted from vactoair.pro & airtovac.pro (Lindler/Landsman/Schlegel).
This Python implementation by A. Bolton, U. Utah, 2011 March.

The two functions of interest are:
  air_wave = airtovac.v2a(vac_wave)
  vac_wave = airtovac.a2v(air_wave)
with all wavelengths in Angstroms.

No correction made for wavelengths below 2000 Angstroms.

"""
c0 = 1.0
c1 = 5.792105e-2
c2 = 238.0185e0
c3 = 1.67917e-3
c4 = 57.362e0
wmin = 2000.0


def wave_to_sigma2(wave):
    """
    Internal support routine to convert to wavevnumber^2.
    """
    return (1.e4 / wave)**2


def conv_factor(sigma2):
    """
    Internal support routine to compute conversion factor.
    """
    return c0 + c1 / (c2 - sigma2) + c3 / (c4 - sigma2)


def v2a(vac_wave):
    """
    Convert vacuum wavelengths in Anstroms to air wavelengths.
    No correction for wavelengths blueward of 2000 Angstroms.
    """
    air_wave = vac_wave / conv_factor(wave_to_sigma2(vac_wave))
    # Only apply conversion above 2000Ang:
    wave_test = vac_wave >= wmin
    air_wave = air_wave * wave_test + vac_wave * (1 - wave_test)
    return air_wave


def a2v(air_wave):
    """
    Convert air wavelengths in Anstroms to vacuum wavelengths.
    No correction for wavelengths blueward of 2000 Angstroms.
    """
    vac_wave = air_wave * conv_factor(wave_to_sigma2(air_wave))
    # Iterate once for accuracy:
    vac_wave = air_wave * conv_factor(wave_to_sigma2(vac_wave))
    # Only apply conversion above 2000Ang:
    wave_test = air_wave >= wmin
    vac_wave = vac_wave * wave_test + air_wave * (1 - wave_test)
    return vac_wave
