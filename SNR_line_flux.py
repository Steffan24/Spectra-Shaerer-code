from modules import np, plt, curve_fit
from variables import dir_basic, output_mass
from SNR_calculation_chi import calcSNR_chi2, calc_line_flux

output_exp_time = [10]
output_masses = [5.5,6.0,6.5,7.0,7.5,8.0]
flux_array = []
SNR_array = []

for i in range(len(output_masses)):
    spectrum = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_masses[i]}_output/spectrum_counts.npy")
    spectrum_flux = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_masses[i]}_output/spectrum_flux.npy")
    spectrum_std = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_masses[i]}_output/spectrum_std.npy")
    wavelength_angstrom = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_masses[i]}_output/wavelength_angstrom.npy")

    SNR, total_counts = calcSNR_chi2(wavelength_angstrom, spectrum, spectrum_std)
    line_flux = calc_line_flux(wavelength_angstrom, spectrum_flux, spectrum_std)
    flux_array.append(line_flux)
    SNR_array.append(SNR)

plt.figure()
plt.scatter(flux_array, SNR_array)
plt.ylabel("\(SNR\)")
plt.xlabel("\(Line flux (erg \cdot\ s^{-1}\cdot cm^{-2})\)")
plt.show()


