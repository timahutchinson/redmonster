import numpy as np

__linelist__ = {
                'OII_1' : 3727.09,
                'OII_2' : 3729.88,
                'H_theta' : 3798.98,
                'H_eta' : 3836.47,
                'H_xi' : 3890.16,
                'H_epsilon' : 3971.20,
                'SII_1' : 4072.30,
                'H_delta' : 4102.89,
                'H_gamma' : 4341.68,
                'OIII_1' : 4364.44,
                'H_beta' : 4862.68,
                'OIII_2' : 4932.60,
                'OIII_3' : 4960.30,
                'OIII_4' : 5008.24,
                'HeI' : 5877.65,
                'OI_1' : 6302.05,
                'OI_2' : 6365.54,
                'NI' : 6529.03,
                'NII' : 6549.86,
                'H_alpha' : 6564.61,
                'NII' : 6585.27,
                'SII_2' : 6718.29,
                'SII_3' : 6732.67
                }

def create_mask(lines, wave):
    wave = np.array(wave)
    mask = np.ones(wave.shape[0])
    for line in lines.values():
        loc = np.abs(wave-line).argmin()
        mask[loc-5:loc+5] = 0
    return mask
