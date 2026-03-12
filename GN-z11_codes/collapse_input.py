from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table, fits, inset_axes, mark_inset,patches, curve_fit, GridSpec, axes3d, CircularAperture,aperture_photometry
import plotting_params

input_file = '/home/steff/hsim/GNz11/model_spectra/logE_crazy_bright.fits'

hdul = fits.open(input_file)
data = hdul[0].data * 10**(-3)
flux = np.nansum(data, axis=(1,2))
header = hdul[0].header
crval3 = header.get("CRVAL3", 0)
cdelt3 = header.get("CDELT3", 1) 
crpix3 = header.get("CRPIX3", 1)
n_wave = data.shape[0]
wavelength = crval3 + cdelt3 * (np.arange(n_wave) - crpix3)
wavelength_angstrom = wavelength #* (10**4)

plt.figure()
plt.step(wavelength_angstrom, flux)
plt.xlabel('$Wavelength (\mathring{A})$')
plt.ylabel('$Flux (ergs \cdot cm^{-2} \cdot s^{-1} \cdot \mathring{A}^{-1})$')
plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/logE_crazy_test/logE_source_check_input_spectra.png')
plt.show()

collapsed_image = np.nansum(data, axis=0)
collapsed_image = np.nan_to_num(collapsed_image)
ny, nx = collapsed_image.shape
pixel_scale_mas = 1
extent = [0, nx * pixel_scale_mas, 0, ny * pixel_scale_mas]
plt.figure()
image = plt.imshow(collapsed_image, origin='lower', cmap='gray', aspect='equal', extent=extent)
plt.xlabel('$x\ (mas)$')
plt.ylabel('$y\ (mas)$')
plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/logE_crazy_test/logE_source_collapsed_image.png')
plt.show()

## collapse around coord

region_size = 2

def extract_central_region(data, region_size, subtract_background=True):
    nlam, ny, nx = data.shape
    cy, cx = 86, 109
    half = region_size // 2

    y1, y2 = cy - half, cy + half
    x1, x2 = cx - half, cx + half

    source_region = data[:, y1:y2, x1:x2]

    positions = [(cx,cy)]
    aperture = CircularAperture(positions, r= region_size // 2)

    counts = np.zeros(nlam)
    #error = np.zeros(nlam)
    for i in range(nlam):
        frame = data[i]     
        phot = aperture_photometry(frame, aperture, method='exact')
        counts[i] = phot['aperture_sum'][0]
        #error[i] = phot['aperture_sum_err'][0]
    return counts#, error

def extract_central_region_2(data, region_size, subtract_background=True):
    nlam, ny, nx = data.shape
    cy, cx = 114, 91
    print(f"cy: {cy}, cx: {cx}")
    half = region_size // 2

    y1, y2 = cy - half, cy + half
    x1, x2 = cx - half, cx + half

    source_region = data[:, y1:y2, x1:x2]

    positions = [(cx,cy)]
    aperture = CircularAperture(positions, r= region_size // 2)

    counts = np.zeros(nlam)
    #error = np.zeros(nlam)
    for i in range(nlam):
        frame = data[i]     
        phot = aperture_photometry(frame, aperture, method='exact')
        counts[i] = phot['aperture_sum'][0]
        #error[i] = phot['aperture_sum_err'][0]
    return counts#, error


He_II = 1640*(1+11.7)
flux_cut = extract_central_region(data, region_size, subtract_background=True)
plt.figure()
plt.plot(wavelength_angstrom, flux_cut)
plt.vlines(x=He_II, ymin=0, ymax=max(flux_cut), colors='red', ls='--')
plt.show()

flux = extract_central_region_2(data, region_size, subtract_background=True)
plt.figure()
plt.plot(wavelength_angstrom, flux)
plt.vlines(x=He_II, ymin=0, ymax=max(flux), colors='red', ls='--')
plt.show()
flux_pop = flux_cut - flux



plt.figure()
plt.plot(wavelength_angstrom, flux_pop)
plt.vlines(x=He_II, ymin=0, ymax=max(flux_pop), colors='red', ls='--')
plt.ylabel("Flux")
plt.xlabel("wavelength")
plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/logE_crazy_test/popII_minus_galaxy.png')
plt.show()


### PLOT EMISSION LINES

# z_orig = 10.6
# z_new = 11.7

# spectral_lines_file = ascii.read('/home/steff/hsim/spectra/spectra/GN-z11/lines')
# line_name = spectral_lines_file['line']
# line_wavelength = spectral_lines_file['l0']

# line_wavelength_shift = (line_wavelength/(1+z_orig)) * (1+z_new)

# def gaussian(x,A, B, C):
#     std = C
#     gaussian = (A/(np.sqrt(2*np.pi)*std))*np.exp(-0.5*((x-B)**2)/std**2)
#     return gaussian


