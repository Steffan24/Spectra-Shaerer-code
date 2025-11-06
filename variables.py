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
output_file = "/home/steff/hsim/HSIM/hsim/output_cubes/27_run/V_27.4_new_reduced.fits"
output_array_init = "/home/steff/hsim/HSIM/hsim/output_cubes"
output_2 = "/home/steff/hsim/zackrisson_pop3_all"
output_array = [f"{output_array_init}/V_24_run3/V_24_reduced.fits",f"{output_array_init}/V25/V_24.9_new_reduced.fits",
                f"{output_2}/code/V_26.2_new_reduced.fits",
                f"{output_array_init}/27_run/V_27.4_new_reduced.fits", f"{output_array_init}/29_run/V_28.7_new_reduced.fits", f"{output_array_init}/30_run/V_29.9_new_reduced.fits" ]

#### FITS HEADER PARAMETERS ###
SIMPLE =' T'
BITPIX = '-32'
NAXIS = 3
NAXIS1 = 300
NAXIS2 = 300
NAXIS3 = 6333
EXTEND =' T'
CTYPE1 = 'RA'
CTYPE2 = 'DEC'
CTYPE3 = 'wavelength'
CUNIT1 = 'mas'
CUNIT2 = 'mas'
CUNIT3 = 'angstrom'
CDELT1 = 1
CDELT2 = 1
CDELT3 = 1.209
CRVAL3 = 17562.5
CRPIX3 = 1
BUNIT = 'erg/s/cm2/AA/arcsec2'
SPECRES = 2.675




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
