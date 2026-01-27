import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import ascii
import plotting_params
from modules import interpolate


### CONTINUUM ####

wavelength_angstrom_lyman = np.arange(14132, 60000,1000) #pre Ly-alpha break
wavelength_angstrom_early = np.arange(5000,14132,1000) # post Ly-alpha break

cont_flux_lyman = (2.998*10**21)* (10**(-0.4*(26+48.60)))/(1000*wavelength_angstrom_lyman**2) # from Prof. Swinbank code: estimating continuum post Ly break

cont_flux_early = np.zeros(len(wavelength_angstrom_early)) # assume 0 cont pre-Ly break

wavelength_full_cont = np.concatenate([wavelength_angstrom_early, wavelength_angstrom_lyman]) # create 1 array

median_cont_wavelength = np.median(wavelength_full_cont) # median for resolution calc

input_res = median_cont_wavelength/(5*7000) # x5 oversampling

flux_full_cont = np.concatenate([cont_flux_early, cont_flux_lyman]) # merge flux models

## interpolate to match x5 oversampling
x_new = np.arange(min(wavelength_full_cont), max(wavelength_full_cont), input_res)
interpolated_wavelengths = interpolate.interp1d(wavelength_full_cont, flux_full_cont)
interpolated_fluxes = np.array(interpolated_wavelengths(x_new))

diff = x_new[28] - x_new[27]
print(f'SPECTRAL RESOLUTION: {diff}')

#plot continuum and save
plt.figure()
plt.plot(x_new, interpolated_fluxes)
plt.xlabel('$Wavelength\ (\mathring{A})$')
plt.ylabel('$Flux\ (erg \cdot cm^{-2} \cdot s^{-1} \cdot \mathring{A}^{-1})$')
plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/model_cont.png')
plt.show()

### ADD SPECTRAL LINES (from Bunker et al 2023) ###

spectral_lines_file = ascii.read('/home/steff/hsim/spectra/spectra/GN-z11/lines')
line_name = spectral_lines_file['line']
line_wavelength = spectral_lines_file['l0']
line_flux = spectral_lines_file['flux'] * 10**(-19)
line_EW = spectral_lines_file['EW']

def flux_calc(line_flux, line_wavelength, wavelength_array):
    M = (10**(8.73)) * 1.99*10**30
    G = 6.67*10**(-11)
    r = 64 * (3.086*10**16)
    sigma_gal = np.sqrt(M*G / r)
    line_wavelength = line_wavelength
    sigma_line = (sigma_gal/(3*10**8)) * (line_wavelength)
    print(f'SIGMA GAL: {sigma_gal} m/s')
    print(f'SIGMA LINE: {sigma_line} m/s')
    exp = np.exp(((wavelength_array - line_wavelength)**2)/(-2*(sigma_line**2)))
    flux = (line_flux/(np.sqrt(2*np.pi)*sigma_line)) * exp
    return flux

wavelength_array = x_new
lines_arrays = []
for i in range(len(line_wavelength)):
    flux_i = flux_calc(line_flux[i], line_wavelength[i], wavelength_array)
    print(f'line flux: {flux_i}')
    lines_arrays.append(flux_i)


total_lines = np.sum(lines_arrays, axis=0)
plt.figure()
plt.plot(wavelength_array, total_lines)
plt.show()
print(f'total lines: {total_lines}')

spectra_flux = interpolated_fluxes + total_lines

print(f'total flux: {spectra_flux}')

## import actual spectra ###

files = ascii.read('/home/steff/hsim/spectra/spectra/GN-z11/infiles')
infiles = files['name']
print(infiles)
nf = len(infiles)

wavelength_data = []
flux_data = []

for p in range(0,nf):

    data = np.genfromtxt(f'/home/steff/hsim/spectra/spectra/GN-z11/{infiles[p]}')#,dtype=None)

    # wl in Angstroms
    wl = data[:,0]
    # flux in erg/cm^2/s/Ang
    flux = data[:,1]

    wavelength_data.append(wl)
    flux_data.append(flux)

wl_all = np.concatenate(wavelength_data)
fl_all = np.concatenate(flux_data)

plt.figure()
plt.plot(x_new, spectra_flux)
plt.step(wl_all, fl_all)
plt.xlabel('$Wavelength\ (\mathring{A})$')
plt.ylabel('$Flux\ (erg \cdot cm^{-2} \cdot s^{-1} \cdot \mathring{A}^{-1})$')
plt.xlim(10000, 30000)
plt.ylim(-0.5*10**(-19), 1.5*10**(-19))
plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/model_spectra.png')
plt.show()


np.save('/home/steff/hsim/GNz11/model_spectra/wavelength.npy', x_new)
print('Should have saved wavelength in GNz11 dir')
np.save('/home/steff/hsim/GNz11/model_spectra/flux.npy', spectra_flux)
print('Should have saved flux in GNz11 dir')









