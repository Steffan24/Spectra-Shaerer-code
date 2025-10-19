# main.py

from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os
from constants import T_sun, c_m, T_100M, M_sun_kg, G, kb, c, h, pc, AU, d, R_sun, M_sun
from variables import ttt, imf, mup, low, sfh, n_single, n_array, M_gauss, d_gauss, save, n
import plotting_params
from functions import import_data, sun_type_star, blackbody, plot, import_lines, gaussian_profile, full_spectra, plot_full_spectra

SED_data = import_data(n, save, ttt, imf, mup, low, sfh, n_single, n_array)

lambda_sun, B, log_B = sun_type_star()

lambda_blackbody, blackbody, log_blackbody = blackbody(T_100M)

plot(n, SED_data, lambda_sun, B, log_B, lambda_blackbody, blackbody, log_blackbody)

age_log, H_beta, H_lya, H_alpha, H_beta_, HeI_4471, HeII_1640, HeII_4686, HeII_3203, HeII_4541 = import_lines(ttt, imf, mup, low, sfh)

flux_H_beta, flux_H_lya, flux_H_alpha, flux_H_4471, flux_HII_1640, flux_HII_3203, flux_HII_4541, flux_HII_4686 = gaussian_profile(M_gauss, d_gauss, SED_data["wavelengths"], age_log, H_beta, H_lya, H_alpha, H_beta_, HeI_4471, HeII_1640, HeII_4686, HeII_3203, HeII_4541)

total_flux_lines = full_spectra(n, SED_data["wavelengths"], flux_H_beta, flux_H_lya, flux_H_alpha, flux_H_4471, flux_HII_1640, flux_HII_3203, flux_HII_4541, flux_HII_4686, SED_data)

plot_full_spectra(n, total_flux_lines, SED_data)
