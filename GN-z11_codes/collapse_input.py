from modules import np, plt, os, ascii, fits
import plotting_params

input_file = '/home/steff/hsim/GNz11/model_spectra/pt_source_input.fits'

hdul = fits.open(input_fits)
data = hdul[0].data
flux = np.nansum(data, axis=(1,2))
header = hdul[0].header
crval3 = header.get("CRVAL3", 0)
cdelt3 = header.get("CDELT3", 1) 
crpix3 = header.get("CRPIX3", 1)
n_wave = data.shape[0]
wavelength = crval3 + cdelt3 * (np.arange(n_wave) - crpix3)
wavelength_angstrom = wavelength * (10**4)

plt.figure()
plt.step(wavelength_angstrom, flux)
plt.xlabel('$Wavelength (\mathring{A})$')
plt.ylabel('$Flux (ergs \cdot )$')
plt.show()
