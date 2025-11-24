# variables.py

from modules import np
from constants import pc, M_sun

### INPUT CUBE PARAMETERS

# pop 3 model params

ttt = 'ge0' #Options - ge0: no mass loss, mdt: with mass loss
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'
n_single = 31 #SED of population (single plotting)
n_array = np.linspace(31,1030, 999, dtype='int') # SED (multiple plotting)

# input cube
M_gauss = (10**6)*M_sun
d_gauss = 100*pc
save = False
n = 'single'
z = 11.5
R = 7000
cube_length = 0.3
input_scale = 1*10**(-3)

# FITS header params
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

### OUTPUT CUBE PARAMETERS
output_scale = 7*10**(-3)

dir_basic = "/home/steff/hsim/HSIM/hsim/output_cubes/salpeter_10ks" #output directories
output_mass = 8.0
corresponding_V = 23.7

#single file analysis
output_file = f"{dir_basic}/M_{output_mass}_output/V_{corresponding_V}_new_reduced.fits"
output_SNR = f"{dir_basic}/M_{output_mass}_output/V_{corresponding_V}_new_reduced_SNR.fits"
output_flux =  f"{dir_basic}/M_{output_mass}_output/V_{corresponding_V}_reduced_flux_cal.fits"
output_std = f"{dir_basic}/M_{output_mass}_output/V_{corresponding_V}_std.fits"

#multiple file analysis
mass_output_array = [7.5,7.0,6.5,6.0,5.5]
V_array = [24.9, 26.2, 27.4, 28.7, 29.9]
output_counts_array = []
output_flux_array = []
output_std_array = []

for i in range(len(mass_output_array)):
    dir_counts = f"{dir_basic}/M_{mass_output_array[i]}_output/V_{V_array[i]}_new_reduced.fits"
    dir_flux = f"{dir_basic}/M_{mass_output_array[i]}_output/V_{V_array[i]}_new_reduced_flux_cal.fits"
    dir_std = dir = f"{dir_basic}/M_{mass_output_array[i]}_output/V_{V_array[i]}_new_std.fits"
    output_counts_array.append(dir_counts)
    output_flux_array.append(dir_flux)
    output_std_array.append(dir_std)
    
output_array_init = "/home/steff/hsim/HSIM/hsim/output_cubes"
output_2 = "/home/steff/hsim/zackrisson_pop3_all"
output_array = [f"{output_array_init}/V_24_run3/V_24_reduced.fits",f"{output_array_init}/V25/V_24.9_new_reduced.fits",
                f"{output_2}/code/V_26.2_new_reduced.fits",
                f"{output_array_init}/27_run/V_27.4_new_reduced.fits", f"{output_array_init}/29_run/V_28.7_new_reduced.fits", f"{output_array_init}/30_run/V_29.9_new_reduced.fits" ]








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
