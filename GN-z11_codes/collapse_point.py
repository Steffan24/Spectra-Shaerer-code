from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table, fits, inset_axes, mark_inset,patches, curve_fit, GridSpec, axes3d, CircularAperture,aperture_photometry

import plotting_params

#SETUP

output_dir = '/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/pt_source'

counts_file = f'{output_dir}/pt_source_input_reduced.fits'
flux_file = f'{output_dir}/pt_source_input_reduced_flux_cal.fits'
std_file = f'{output_dir}/pt_source_input_std.fits'

#IMPORT DATA
print("Importing files:")
print("counts")
#counts
hdul = fits.open(counts_file)
data = hdul[0].data
counts_data = np.nansum(data, axis=(1,2))
header = hdul[0].header
crval3 = header.get("CRVAL3", 0)
cdelt3 = header.get("CDELT3", 1) 
crpix3 = header.get("CRPIX3", 1)
n_wave = data.shape[0]
wavelength = crval3 + cdelt3*(np.arange(n_wave) - crpix3)
counts_wavelength = wavelength*(10**4)
print("flux")
#flux
hdul = fits.open(flux_file)
data = hdul[0].data
flux_data = np.nansum(data, axis=(1,2))
header = hdul[0].header
crval3 = header.get("CRVAL3", 0)
cdelt3 = header.get("CDELT3", 1) 
crpix3 = header.get("CRPIX3", 1)
n_wave = data.shape[0]
wavelength = crval3 + cdelt3*(np.arange(n_wave) - crpix3)
flux_wavelength = wavelength*(10**4)
print("standard deviations")
#std
hdul = fits.open(std_file)
data = hdul[0].data
std_data = np.nansum(data, axis=(1,2))
header = hdul[0].header
crval3 = header.get("CRVAL3", 0)
cdelt3 = header.get("CDELT3", 1) 
crpix3 = header.get("CRPIX3", 1)
n_wave = data.shape[0]
wavelength = crval3 + cdelt3*(np.arange(n_wave) - crpix3)
std_wavelength = wavelength*(10**4)

## APERTURE

extraction_region = 6

def extract_central_region(data, region_size, data_std, subtract_background=True):
    nlam, ny, nx = data.shape
    cy, cx = ny // 2, nx // 2
    half = region_size // 2

    y1, y2 = cy - half, cy + half
    x1, x2 = cx - half, cx + half

    source_region = data[:, y1:y2, x1:x2]

    positions = [(cx,cy)]
    aperture = CircularAperture(positions, r= region_size // 2)

    counts = np.zeros(nlam)
    error = np.zeros(nlam)
    for i in range(nlam):
        frame = data[i]     
        phot = aperture_photometry(frame, aperture, error = data_std[i], method='exact')
        counts[i] = phot['aperture_sum'][0]
        error[i] = phot['aperture_sum_err'][0]
    return counts, error

hdul = fits.open(counts_file)
data = hdul[0].data
hdul_std = fits.open(std_file)
std_data = hdul_std[0].data
counts, std = extract_central_region(data, extraction_region, std_data, subtract_background = True)
hdul = fits.open(flux_file)
data = hdul[0].data
flux, flux_std = extract_central_region(data, extraction_region, std_data, subtract_background = True)

np.save(f"{output_dir}/counts.npy", counts)
np.save(f"{output_dir}/flux.npy", flux)
np.save(f"{output_dir}/std.npy", std)
np.save(f"{output_dir}/wavelength.npy", counts_wavelength)

### SKYLINES

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

### PLOT

fig, (ax1,ax2) = plt.subplots(2,1, height_ratios = [1,0.3], sharex=True)
ax1.plot(flux_wavelength, flux)
ax2.set_xlabel("$Wavelength\ (\mathring{A})$")
ax1.set_ylabel("$Flux\ (erg \cdot cm^{-2} \cdot s^{-1} \cdot \mathring{A}^{-1})$")
ax1.tick_params(labelbottom=False)
ax2.plot(opacity_data["wavelength"], opacity_data["transmittance"], c='black')
ax2.set_ylabel("Transmission")
ax2.set_xlim(min(counts_wavelength), max(counts_wavelength))
plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/pt_source/output_spectra_raw_flux.png")
plt.show()

fig, (ax1,ax2) = plt.subplots(2,1, height_ratios = [1,0.3], sharex=True)
ax1.plot(counts_wavelength, counts)
ax1.tick_params(labelbottom=False)
ax2.set_xlabel("$Wavelength\ (\mathring{A})$")
ax1.set_ylabel("$Counts$")
ax2.plot(opacity_data["wavelength"], opacity_data["transmittance"], c='black')
ax2.set_ylabel("Transmission")
ax2.set_xlim(min(counts_wavelength), max(counts_wavelength))
plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/pt_source/output_spectra_raw_counts.png")
plt.show()


## Spectrum atmosphere correction

## interpolate transmission
x_new = counts_wavelength
opacity_wave = np.array(opacity_data["wavelength"])
opacity_trans = np.array(opacity_data["transmittance"])
interpolated_wavelengths = interpolate.interp1d(opacity_wave, opacity_trans)
interpolated_trans = np.array(interpolated_wavelengths(x_new))

corrected_counts = counts / interpolated_trans

plt.figure()
plt.plot(counts_wavelength, corrected_counts)
plt.show()



