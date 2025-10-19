# File 2: import_emission_lines.py
# Import emission lines data and apply a gaussian fit

### MODULES ###
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from astropy.io import ascii
import latex
import os

### PLOTTING PARAMS ###
plt.rcParams["font.size"] = 30
plt.rcParams["legend.fontsize"] = 30
plt.rcParams["figure.frameon"] = False
plt.rcParams["figure.titlesize"] = 16
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.subplot.wspace"] = 0.03
plt.rcParams["legend.frameon"] = False
plt.rcParams["axes.linewidth"] = 1.75

### CONSTANTS ###
M_sun = 1.989*10**30
pc = 3.08*10**16
c = 3*10**8

### VARIABLES (unchanged as of 07/10/25)  ###
ttt = 'ge0' #Options - ge0: no mass loss, mdt: with mass loss
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'
M = 10**6
R = 100*pc 

n = 'single'

n_single = 31 #SED of population (single plotting)

n_array = np.linspace(31, 1030, 999, dytype='int')

save = False


### IMPORT LINES ###

def import_lines():
    #read data
    file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.22"
    if os.path.exists(file_loc):
        data = ascii.read(file_loc,guess = True, data_start = 0)
        print(data)
        age_log = data['col1']
        H_beta = data['col5']
        print(H_beta)
        H_lya = data['col7'] * H_beta
        H_alpha = data['col9'] * H_beta
        H_beta_ = data['col11'] * H_beta
        HeI_4471 = data['col13'] * H_beta
        HeII_1640 = data['col15'] * H_beta
        HeII_4686 = data['col17'] * H_beta
        HeII_3203 = data['col19'] * H_beta
        HeII_4541 = data['col21'] * H_beta
    return age_log, H_beta, H_lya, H_alpha, H_beta_, HeI_4471, HeII_1640, HeII_4686, HeII_3203, HeII_4541

### APPLY A GAUSSIAN PROFILE ###

def gaussian_profile(M, R, wavelength):
    print(f"INPUT PARAMS: M = {M} kg, R = {R} m")
    sigma_gal =np.sqrt(M*G/R) # m/s
    print(f"SIGMA GAL: {sigma_gal} m/s")
    #array in peak [H_beta, H_lya, H_alpha, HEI_4471, HeII1640, HeII_4686, HeII_3203, HeII_4541]]
    lambda_peak_array = [4861, 1215, 6563, 4471, 1640, 4686, 3203, 4541]
    lambda_peak_array = np.array(lambda_peak_array)
    sigma_line = (sigma_gal/c_m) * np.array(lambda_peak_array)
    print(f"SIGMA LINE: {sigma_line} Angstrom")
    H_beta_peak_Intensity = (H_beta[0])/(4*np.pi*(d**2)*M_sun)
    print(f"H_beta_peak_Intensity: {np.log10(H_beta_peak_Intensity) + 30}")
    H_lya_peak_intensity = (H_lya[0])/(4*np.pi*(d**2)*M_sun)
    print(f"H_lya_peak_Intensity: {np.log10(H_lya_peak_intensity) + 30}")
    H_alpha_peak_intensity = (H_alpha[0])/(4*np.pi*(d**2)*M_sun)
    print(f"H_alpha_peak_Intensity: {np.log10(H_alpha_peak_intensity) + 30}")
    HeI_4471_peak_intensity = (HeI_4471[0])/(4*np.pi*(d**2)*M_sun)
    print(f"H_4471_peak_Intensity: {np.log10(HeI_4471_peak_intensity) + 30}")
    HeII_1640_peak_intensity = (HeII_1640[0])/(4*np.pi*(d**2)*M_sun)
    print(f"H_1640_peak_Intensity: {np.log10(HeII_1640_peak_intensity) + 30}")
    HeII_3203_peak_intensity = (HeII_3203[0])/(4*np.pi*(d**2)*M_sun)
    print(f"H_3203_peak_Intensity: {np.log10(HeII_3203_peak_intensity) + 30}")
    HeII_4541_peak_intensity = (HeII_4541[0])/(4*np.pi*(d**2)*M_sun)
    print(f"H_4541_peak_Intensity: {np.log10(HeII_4541_peak_intensity) + 30}")
    HeII_4686_peak_intensity = (HeII_4686[0])/(4*np.pi*(d**2)*M_sun)
    print(f"H_4686_peak_Intensity: {np.log10(HeII_4686_peak_intensity) + 30}")
    flux_H_beta = []
    flux_H_lya = []
    flux_H_alpha = []
    flux_H_4471 = []
    flux_HII_1640 = []
    flux_HII_4686 = []
    flux_HII_3203 = []
    flux_HII_4541 = []
    for i in range(len(wavelength)):
        exponential_1 = np.exp(((wavelength[i] - lambda_peak_array[0])**2)/(-2*(sigma_line[0])**2))
        flux_H_beta_1 = ((H_beta_peak_Intensity/(np.sqrt(2*np.pi)*sigma_line[0])) * exponential_1)
        flux_H_beta.append(flux_H_beta_1)
        exponential_2 = np.exp(((wavelength[i] - lambda_peak_array[1])**2)/(-2*(sigma_line[1])**2))
        flux_H_lya_1 = ((H_lya_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[1]))*exponential_2)
        flux_H_lya.append(flux_H_lya_1)
        exponential_3 = np.exp(((wavelength[i] - lambda_peak_array[2])**2)/(-2*(sigma_line[2])**2))
        flux_H_alpha_1 = ((H_alpha_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[2]))*exponential_3)
        flux_H_alpha.append(flux_H_alpha_1)
        exponential_4 = np.exp(((wavelength[i] - lambda_peak_array[3])**2)/(-2*(sigma_line[3])**2))
        flux_H_4471_1 = ((HeI_4471_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[3]))*exponential_4)
        flux_H_4471.append(flux_H_4471_1)
        exponential_5 = np.exp(((wavelength[i] - lambda_peak_array[4])**2)/(-2*(sigma_line[4])**2))
        flux_H_1640_1 = ((HeII_1640_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[4]))*exponential_5)
        flux_HII_1640.append(flux_H_1640_1)
        exponential_6 = np.exp(((wavelength[i] - lambda_peak_array[5])**2)/(-2*(sigma_line[5])**2))
        flux_H_4686_1 = ((HeII_4686_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[5]))*exponential_6)
        flux_HII_4686.append(flux_H_4686_1)
        exponential_7 = np.exp(((wavelength[i] - lambda_peak_array[6])**2)/(-2*(sigma_line[6])**2))
        flux_H_3203_1 = ((HeII_3203_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[6]))*exponential_7)
        flux_HII_3203.append(flux_H_3203_1)
        exponential_8 = np.exp(((wavelength[i] - lambda_peak_array[7])**2)/(-2*(sigma_line[7])**2))
        flux_H_4541_1 = ((HeII_4541_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[7]))*exponential_8)
        flux_HII_4541.append(flux_H_4541_1)
    flux_H_beta = np.array(flux_H_beta)
    flux_H_lya = np.array(flux_H_lya)
    flux_H_alpha = np.array(flux_H_alpha)
    flux_H_4471 = np.array(flux_H_4471)
    flux_HII_1640 = np.array(flux_HII_1640)
    flux_HII_4686 = np.array(flux_HII_4686)
    flux_HII_3203 = np.array(flux_HII_3203)
    flux_HII_4541 = np.array(flux_HII_4541)
    return flux_H_beta, flux_H_lya, flux_H_alpha, flux_H_4471, flux_HII_1640, flux_HII_3203, flux_HII_4541, flux_HII_4686


### MERGE DATA ###

def full_spectra(n, wavelength):
    if n == 'single':
        total_flux_lines = []
        for i in range(len(wavelength)):
            flux = (SED_data["SED_flux"][i]) + flux_H_beta[i] + flux_H_lya[i] + flux_H_alpha[i] + flux_H_4471[i] + flux_HII_1640[i] + flux_HII_3203[i] + flux_HII_4541[i] + flux_HII_4686[i]
            total_flux_lines.append(flux)
        total_flux_lines = np.array(total_flux_lines)
    if n == 'multi':
        
        
        


### RUN SCRIPT ###

from import_spectrum import import_data

SED_data = import_data(n, save)

wavelength = SED_data["wavelengths"]

import_lines()
flux_H_beta, flux_H_lya, flux_H_alpha, flux_H_4471, flux_HII_1640, flux_HII_3203, flux_HII_4541, flux_HII_4686 = gaussian_profile(M, R, wavelength)
