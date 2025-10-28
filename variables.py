# variables.py

from modules import np
from constants import pc, M_sun

ttt = 'ge0' #Options - ge0: no mass loss, mdt: with mass loss
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'

n_single = 31 #SED of population (single plotting)

n_array = np.linspace(31,1030, 999, dtype='int') # SED (multiple plotting)

M_gauss = (10**6)*M_sun
d_gauss = 100*pc

save = False
n = 'single'
z = 11.5
R = 7000

cube_length = 0.3
input_scale = 1*10**(-3)
#fits_location =

#### FITS HEADER PARAMETERS ###
SIMPLE =' T / conforms to FITS standard'
BITPIX = '-64'
NAXIS =' 3'
NAXIS1 = '300'
NAXIS2 = '300'
NAXIS3 = '8480'
EXTEND = 'T'
CTYPE1 = 'RA'
CTYPE2 = 'DEC'
CTYPE3 = 'wavelength'
CUNIT1 = 'mas'
CUNIT2 = 'mas'
CUNIT3 = 'angstrom'
CDELT1 = '1'
CDELT2 = '1'
CDELT3 = '0.093'
CRPIX3 = '12562.5'
BUNIT = 'erg/s/cm2/AA/arcsec2'
SPECRES = '1.1625'




#HARMONI spectograph resolution data (microm)

LR_IZJ_min = 0.811
LR_IZJ_max = 1.369
LR_HK_min = 1.450
LR_HK_max = 2.450

MR_IZ_min = 0.830
MR_IZ_max = 1.050
MR_J_min = 1.046
MR_J_max = 1.324
MR_H_min = 1.435
MR_H_max = 1.815
MR_K_min = 1.951
MR_K_max = 2.469
