from modules import np, plt, os, ascii, fits
import plotting_params
from variables_old import LR_IZJ_min, LR_IZJ_max,LR_HK_min,LR_HK_max,MR_IZ_min,MR_IZ_max,MR_J_min,MR_J_max,MR_H_min,MR_H_max,MR_K_min,MR_K_max
from data_cube_params import SIMPLE, BITPIX, NAXIS, EXTEND, CTYPE1, CTYPE2, CTYPE3, CUNIT1, CUNIT2, CUNIT3, CDELT1, CDELT2, CRPIX3, BUNIT

spectrum_wavelength = np.load('/home/steff/hsim/GNz11/model_spectra/wavelength.npy')
spectrum_flux = np.load('/home/steff/hsim/GNz11/model_spectra/flux.npy')


def import_opacity():
    file_loc = "/home/steff/hsim/zackrisson_pop3_all/chtrans_nir_18504ft_lon30d_pwv0900um_zenith_smoothed2e-3um-1.dat"
    data = ascii.read(file_loc, guess = True, data_start = 2)
    wavelength_opacity = np.array(data['Wavelength'])
    wavelength_opacity_angstroms = wavelength_opacity * (10**4)
    transmittance_opacity = np.array(data['Transmittance'])
    opacity_data = {"wavelength": wavelength_opacity_angstroms,
                    "transmittance": transmittance_opacity}
    return opacity_data

opacity_data = import_opacity()

z_orig = 10.6
z_new = 11.7
spectral_lines_file = ascii.read('/home/steff/hsim/spectra/spectra/GN-z11/lines')
line_name = spectral_lines_file['line']
line_wavelength = spectral_lines_file['l0']
line_flux = spectral_lines_file['flux'] * 10**(-19)
line_EW = spectral_lines_file['EW']

# norm_wavelength = line_wavelength / (1+z_orig)
# new_wavelength = norm_wavelength * (1+z_new)
# plt.figure()
# plt.plot(opacity_data["wavelength"], opacity_data["transmittance"], zorder=3)
# for i in range(len(new_wavelength)):
#     plt.vlines(new_wavelength[i], ymin=-0.05, ymax=1.05, linestyle='--', edgecolors='red', zorder=2)
# plt.ylabel("Transmition")
# plt.xlabel("$Wavelength (\mathring{A})$")
# plt.xlim(14500, 24500)
# plt.ylim(-0.05, 1.2)
# plt.text(15000, 0.05, '$z = 9.7$', bbox=dict(edgecolor='black', fc = 'white'))
# plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/z_9_7.png')
# plt.show()

## REDSHIFT LINES

spectrum_wavelength = (spectrum_wavelength / (1+z_orig)) * (1+z_new)

## plot spectrum + harmoni config

### PLOT SPECTRUM VS HARMONI CONFIGURATIONS

def import_harmoni_res(LR_IZJ_min, LR_IZJ_max,LR_HK_min,LR_HK_max,MR_IZ_min,MR_IZ_max,MR_J_min,MR_J_max,MR_H_min,MR_H_max,MR_K_min,MR_K_max):
    LR_IZJ = np.linspace(LR_IZJ_min, LR_IZJ_max,3)* 10**4
    LR_HK = np.linspace(LR_HK_min, LR_HK_max,3)* 10**4
    MR_IZ = np.linspace(MR_IZ_min, MR_IZ_max,3)* 10**4
    MR_J = np.linspace(MR_J_min, MR_J_max,3)* 10**4
    MR_H = np.linspace(MR_H_min, MR_H_max,3)* 10**4
    MR_K = np.linspace(MR_K_min, MR_K_max,3)* 10**4
    return LR_IZJ, LR_HK, MR_IZ, MR_J, MR_H, MR_K

LR_IZJ, LR_HK, MR_IZ, MR_J, MR_H, MR_K = import_harmoni_res(LR_IZJ_min, LR_IZJ_max,LR_HK_min,LR_HK_max,MR_IZ_min,MR_IZ_max,MR_J_min,MR_J_max,MR_H_min,MR_H_max,MR_K_min,MR_K_max)

yl = -0.2*10**(-19)
ym = -0.4*10**(-19)

plt.figure()
plt.step(spectrum_wavelength, spectrum_flux) # plot spectrum
# add HARMONI configurations
plt.plot(LR_IZJ, np.ones(3)*yl)
plt.text(11000, -0.125*10**(-19), 'IZJ')
plt.plot(LR_HK, np.ones(3)*yl)
plt.text(19500, -0.125*10**(-19), 'HK')
plt.plot(MR_IZ, np.ones(3)*ym)
plt.text(8500, -0.375*10**(-19), 'IZ')
plt.plot(MR_J, np.ones(3)*ym)
plt.text(11750, -0.375*10**(-19), 'J')
plt.plot(MR_H, np.ones(3)*ym)
plt.text(16500, -0.375*10**(-19), 'H')
plt.plot(MR_K, np.ones(3)*ym)
plt.text(22500, -0.375*10**(-19), 'K')
#label emission lines
plt.text(((1.157*10**4)/ (1+z_orig)) * (1+z_new), 0.970*10**(-19), r'$Ly\alpha$')
plt.text(((1.673*10**4)/ (1+z_orig)) * (1+z_new), 1.006*10**(-19), '$NIV$')
plt.text(((2.236*10**4)/ (1+z_orig)) * (1+z_new), 0.482*10**(-19), '$CIII$')
plt.text(((2.010*10**4)/ (1+z_orig)) * (1+z_new), 0.601*10**(-19), '$NIII$')
plt.text(((1.951*10**4)/ (1+z_orig)) * (1+z_new), 1.157*10**(-19), '$OIII$')
plt.text(((1.482*10**4)/ (1+z_orig)) * (1+z_new), 0.799*10**(-19), '$HeII$')
plt.text(((1.466*10**4)/ (1+z_orig)) * (1+z_new), 0.404*10**(-19), '$NIV$')
plt.annotate("", xytext=(((2.053*10**4)/ (1+z_orig)) * (1+z_new), 1.157*10**(-19)), xy=(((1.924*10**4)/ (1+z_orig)) * (1+z_new), 0.326*10**(-19)),
            arrowprops=dict(arrowstyle="->"))
plt.annotate("", xytext=(((1.593*10**4)/ (1+z_orig)) * (1+z_new), 0.799*10**(-19)), xy=(((1.888*10**4)/ (1+z_orig)) * (1+z_new), 0.326*10**(-19)),
            arrowprops=dict(arrowstyle="->"))
plt.annotate("", xytext=(((1.568*10**4)/ (1+z_orig)) * (1+z_new), 0.404*10**(-19)), xy=(((1.79*10**4)/ (1+z_orig)) * (1+z_new), 0.33*10**(-19)),
            arrowprops=dict(arrowstyle="->"))
#redshift box
plt.text(0.7*10**4, 1.25*10**(-19), '$z=11.7$', bbox=dict(edgecolor='black', fc = 'None'))
plt.ylim(-0.5*10**(-19), 1.5*10**(-19))
plt.xlim(6000, 26000)
plt.xlabel('$Wavelength (\mathring{A})$')
plt.ylabel('$Flux\ (erg \cdot cm^{-2} \cdot s^{-1} \cdot \mathring{A}^{-1})$')
plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/spectra_with_harmoni_after_redshift.png')
plt.show()



## CREATE DATE CUBE (1 SPAXEL INJECTION)

mask_wavelength = (spectrum_wavelength > 15000) & (spectrum_wavelength < 24000) #mask to relevant input range
masked_flux = np.array(spectrum_flux[mask_wavelength])
masked_wavelength = np.array(spectrum_wavelength[mask_wavelength])

cube_length = 0.2
input_scale = 1*10**(-3)

#setup
n_spaxels = int(round(cube_length / input_scale))
n_wavelength = len(masked_wavelength)
min_wavelength = min(masked_wavelength)
central_spaxel = n_spaxels // 2
spectral_res = masked_wavelength[26] - masked_wavelength[25]
print(f"Number of wavelength points: {n_wavelength}")
print(f"Number of spaxels: {n_spaxels}")
print(f"Minimum wavelength: {min_wavelength}")
print(f"central spaxel no: {central_spaxel}")
print(f"spectral resolution: {spectral_res}")

print("Creating input cube")
flux_pas = masked_flux / (input_scale**2) # convert to flux per arcsec
cube = np.zeros((n_wavelength, n_spaxels, n_spaxels), dtype=np.float32)
#inject into single spaxel
#cube[:, central_spaxel, central_spaxel] = flux_pas.astype(np.float32)

# hdu = fits.PrimaryHDU(data=cube) # create cube
# #header adjustments
# header = hdu.header
# header.clear()
# header["SIMPLE"] = bool(SIMPLE)
# header["BITPIX"] = bool(BITPIX)
# header["NAXIS"] = NAXIS
# header["NAXIS1"] = cube.shape[2]
# header["NAXIS2"] = cube.shape[1]
# header["NAXIS3"] = cube.shape[0]
# header["EXTEND"] = bool(EXTEND)
# header["CTYPE1"] = CTYPE1
# header["CTYPE2"] = CTYPE2
# header["CTYPE3"] = CTYPE3
# header["CUNIT1"] = CUNIT1
# header["CUNIT2"] = CUNIT2
# header["CUNIT3"] = CUNIT3
# header["CDELT1"] = CDELT1
# header["CDELT2"] = CDELT2
# header["CRVAL3"] = min_wavelength
# header["CDELT3"] = spectral_res
# header["CRPIX3"] = CRPIX3
# header["BUNIT"] = BUNIT
# header["SPECRES"] = 4*spectral_res
# hdul = fits.HDUList([hdu])
# hdul.writeto(f'/home/steff/hsim/GNz11/model_spectra/pt_source_input.fits', overwrite = True)
# print("Input cube saved")

## CREATE RESOLVED SOURCE

from scipy import datasets
from scipy.ndimage import gaussian_filter

#single spaxel injection
cube[:, central_spaxel, central_spaxel] = flux_pas.astype(np.float32)

#new_cube = np.zeros((n_wavelength, n_spaxels, n_spaxels), dtype=np.float32)

half_light_radius = 0.016 #arcsec
resolved_fwhm = half_light_radius*2 #arcsec
sigma = resolved_fwhm / (2*np.sqrt(2*np.log(2))) # arcsec
#convert to pixels using spatial resolution (angstrom/pixel)
spatial_res = 0.001 #arcsec/pixel
sigma_pixel = sigma/spatial_res
print(f"gaussian sigma: {sigma_pixel}")

new_cube = gaussian_filter(cube, sigma=(0, sigma_pixel, sigma_pixel)).astype(np.float32)

hdu = fits.PrimaryHDU(data=new_cube) # create cube
#header adjustments
header = hdu.header
header.clear()
header["SIMPLE"] = bool(SIMPLE)
header["BITPIX"] = bool(BITPIX)
header["NAXIS"] = NAXIS
header["NAXIS1"] = cube.shape[2]
header["NAXIS2"] = cube.shape[1]
header["NAXIS3"] = cube.shape[0]
header["EXTEND"] = bool(EXTEND)
header["CTYPE1"] = CTYPE1
header["CTYPE2"] = CTYPE2
header["CTYPE3"] = CTYPE3
header["CUNIT1"] = CUNIT1
header["CUNIT2"] = CUNIT2
header["CUNIT3"] = CUNIT3
header["CDELT1"] = CDELT1
header["CDELT2"] = CDELT2
header["CRVAL3"] = min_wavelength
header["CDELT3"] = spectral_res
header["CRPIX3"] = CRPIX3
header["BUNIT"] = BUNIT
header["SPECRES"] = 3*spectral_res
hdul = fits.HDUList([hdu])
hdul.writeto(f'/home/steff/hsim/GNz11/model_spectra/gauss_source_input.fits', overwrite = True)
print("Input cube saved")







