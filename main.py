# main.py

from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, mticker
from constants import T_sun, c_m, T_100M, M_sun_kg, G, kb, c, h, pc, AU, d, R_sun, M_sun
from variables import ttt, imf, mup, low, sfh, n_single, n_array, M_gauss, d_gauss, save, n, LR_IZJ_min, LR_IZJ_max,LR_HK_min,LR_HK_max,MR_IZ_min,MR_IZ_max,MR_J_min,MR_J_max,MR_H_min,MR_H_max,MR_K_min,MR_K_max, z
import plotting_params
from functions import import_data, sun_type_star, blackbody, plot, import_lines, gaussian_profile, full_spectra, plot_full_spectra, redshifting, AB_magnitude_conversion, import_harmoni_res, plot_spectra_redshifted, interpolate_SED, import_opacity, import_OH

SED_data = import_data(n, save, ttt, imf, mup, low, sfh, n_single, n_array)
opacity_data = import_opacity()
skyline_data = import_OH()

lambda_sun, B, log_B = sun_type_star()

lambda_blackbody, blackbody, log_blackbody = blackbody(T_100M)

plot(n, SED_data, lambda_sun, B, log_B, lambda_blackbody, blackbody, log_blackbody)

SED_data = interpolate_SED(SED_data, n)

age_log, H_beta, H_lya, H_alpha, H_beta_, HeI_4471, HeII_1640, HeII_4686, HeII_3203, HeII_4541 = import_lines(ttt, imf, mup, low, sfh)

flux_H_beta, flux_H_lya, flux_H_alpha, flux_H_4471, flux_HII_1640, flux_HII_3203, flux_HII_4541, flux_HII_4686 = gaussian_profile(M_gauss, d_gauss, SED_data["wavelengths"], age_log, H_beta, H_lya, H_alpha, H_beta_, HeI_4471, HeII_1640, HeII_4686, HeII_3203, HeII_4541)

total_flux_lines = full_spectra(n, SED_data["wavelengths"], flux_H_beta, flux_H_lya, flux_H_alpha, flux_H_4471, flux_HII_1640, flux_HII_3203, flux_HII_4541, flux_HII_4686, SED_data)

flux_z, wavelength_z = redshifting(n,total_flux_lines, SED_data, z)

plot_full_spectra(n, total_flux_lines, SED_data, wavelength_z)



flux_zab = AB_magnitude_conversion(flux_z, wavelength_z)

LR_IZJ, LR_HK, MR_IZ, MR_J, MR_H, MR_K = import_harmoni_res(LR_IZJ_min, LR_IZJ_max,LR_HK_min,LR_HK_max,MR_IZ_min,MR_IZ_max,MR_J_min,MR_J_max,MR_H_min,MR_H_max,MR_K_min,MR_K_max)

plot_spectra_redshifted(flux_z, wavelength_z, flux_zab, LR_IZJ, LR_HK, MR_IZ, MR_J, MR_H, MR_K, opacity_data, skyline_data)

