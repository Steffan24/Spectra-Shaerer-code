from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table, fits, inset_axes, mark_inset,patches, curve_fit, GridSpec, axes3d, CircularAperture,aperture_photometry

import plotting_params
import matplotlib as mpl

#SETUP

output_dir = '/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/resolved_source_redshifted'

counts_file = f'{output_dir}/gauss_source_input_reduced.fits'
flux_file = f'{output_dir}/gauss_source_input_reduced_flux_cal.fits'
std_file = f'{output_dir}/gauss_source_input_std.fits'

output_scale = 7*10**(-3)

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
data = hdul[0].data * (output_scale**2) * (10**(-4))
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
std_data = np.sqrt(np.nansum(data**2, axis=(1,2)))
header = hdul[0].header
crval3 = header.get("CRVAL3", 0)
cdelt3 = header.get("CDELT3", 1) 
crpix3 = header.get("CRPIX3", 1)
n_wave = data.shape[0]
wavelength = crval3 + cdelt3*(np.arange(n_wave) - crpix3)
std_wavelength = wavelength*(10**4)

## APERTURE

extraction_region = 36

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
data = hdul[0].data* (output_scale**2) * (10**(-4))
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

aperture_size = np.linspace(1,20,19,dtype=int)*2
hdul = fits.open(counts_file)
data = hdul[0].data
hdul_std = fits.open(std_file)
std_data = hdul_std[0].data

# total_counts = []
# SNR_array = []

# plt.figure()
# for i in range(len(aperture_size)):
#     print(aperture_size[i])
#     counts_test, std_test = extract_central_region(data, aperture_size[i], std_data, subtract_background = True)
#     counts_test = np.nansum(counts_test)
#     total_counts.append(counts_test)
#     std_test = np.sqrt(np.nansum(std_test**2))
#     SNR = counts_test/std_test
#     SNR_array.append(SNR)

# plt.plot(aperture_size,SNR_array)
# plt.xlabel("Aperture diameter (pixels)")
# plt.ylabel("SNR")
# plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/gaussian_resolved_plots/diameter_SNR_collapsed.png')
# plt.show()

# plt.figure()
# plt.plot(aperture_size, total_counts)
# plt.xlabel("Aperture diameter (pixels)")
# plt.ylabel("Total Counts")
# plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/gaussian_resolved_plots/counts_diameter_collapsed.png')
# plt.show()

# collapsed_image = np.nansum(data, axis=0)
# collapsed_image = np.nan_to_num(collapsed_image)
# ny, nx = collapsed_image.shape
# pixel_scale_mas = 7
# extent = [0, nx * pixel_scale_mas, 0, ny * pixel_scale_mas]
# fig,ax = plt.subplots()
# image = ax.imshow(collapsed_image, origin='lower', cmap='gray', aspect='equal', extent=extent)
# circ = patches.Circle((101,101),36, alpha=0.8, fc='None',ec='red')
# ax.add_patch(circ)
# ax.set_xlabel('$x\ (mas)$')
# ax.set_ylabel('$y\ (mas)$')
# plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/gaussian_resolved_plots/gauss_source_collapsed_image.png')
# plt.show()
    
    

# fig, (ax1,ax2) = plt.subplots(2,1, height_ratios = [1,0.3], sharex=True)
# ax1.plot(flux_wavelength, flux)
# ax2.set_xlabel("$Wavelength\ (\mathring{A})$")
# ax1.set_ylabel("$Flux\ (erg \cdot cm^{-2} \cdot s^{-1} \cdot \mathring{A}^{-1})$")
# ax1.tick_params(labelbottom=False)
# ax2.plot(opacity_data["wavelength"], opacity_data["transmittance"], c='black')
# ax2.set_ylabel("Transmission")
# ax2.set_xlim(min(counts_wavelength), max(counts_wavelength))
# plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/gaussian_resolved_plots/output_spectra_raw_flux.png")
# plt.show()

# fig, (ax1,ax2) = plt.subplots(2,1, height_ratios = [1,0.3], sharex=True)
# ax1.plot(counts_wavelength, counts)
# ax1.tick_params(labelbottom=False)
# ax2.set_xlabel("$Wavelength\ (\mathring{A})$")
# ax1.set_ylabel("$Counts$")
# ax2.plot(opacity_data["wavelength"], opacity_data["transmittance"], c='black')
# ax2.set_ylabel("Transmission")
# ax2.set_xlim(min(counts_wavelength), max(counts_wavelength))
# plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/gaussian_resolved_plots/output_spectra_raw_counts.png")
# plt.show()

### PLOT EMISSION LINES

z_orig = 10.6
z_new = 11.7

spectral_lines_file = ascii.read('/home/steff/hsim/spectra/spectra/GN-z11/lines')
line_name = spectral_lines_file['line']
line_wavelength = spectral_lines_file['l0']

line_wavelength_shift = (line_wavelength/(1+z_orig)) * (1+z_new)

def gaussian(x,A, B, C):
    std = C
    gaussian = (A/(np.sqrt(2*np.pi)*std))*np.exp(-0.5*((x-B)**2)/std**2)
    return gaussian

## SNR of lines ##

def chi_squared(x,y, std, popt):
    chi = ((y - gaussian(x, *popt))**2)/(std**2)
    return np.sum(chi)

def chi_squared_cont(x,y,std, straight_fit):
    chi = ((y - straight_fit)**2) / (std**2)
    #print(chi)
    return np.nansum(chi)

def SNR(wavelength,line_wavelength, counts, std):
    mask_gaussian = (wavelength > line_wavelength - 30) & (wavelength < line_wavelength + 30)
    wavelength_gaussian = wavelength[mask_gaussian]
    counts_gaussian = counts[mask_gaussian]
    std_gaussian = std[mask_gaussian]
    center = line_wavelength
    sigma_guess = 11

    p0 = [np.nanmax(counts_gaussian) - np.nanmedian(counts_gaussian), center, sigma_guess]
    bounds = ([0, center-15, 10], [np.inf, center+15, 13])
    popt, pcov, infodict, mesg, ier = curve_fit(gaussian, wavelength_gaussian,counts_gaussian,sigma = std_gaussian,absolute_sigma = True, p0=p0, bounds=bounds, full_output = True)
    chi_squared_gaussian = chi_squared(wavelength_gaussian, counts_gaussian, std_gaussian, popt) #/ len(wavelength_gaussian)
    #print(f"chi squared gaussian: {chi_squared_gaussian}")

    mask_cont = (wavelength_gaussian > popt[1] + 1.7*popt[2]) | (wavelength_gaussian < popt[1] - 1.7*popt[2])
    wavelength_cont = wavelength_gaussian[mask_cont]
    counts_cont = counts_gaussian[mask_cont]
    std_cont = std_gaussian[mask_cont]
    straight_fit = np.ones(len(wavelength_cont)) * np.nanmedian(np.array(counts_cont))
    chi_squared_con = chi_squared_cont(wavelength_cont, counts_cont, std_cont, straight_fit)# / len(wavelength_cont)
    #print(f"chi squared continuum: {chi_squared_con}")

    # plt.figure()
    # plt.step(wavelength_gaussian, counts_gaussian)
    # plt.plot(wavelength_gaussian, gaussian(wavelength_gaussian, *popt))
    # plt.scatter(wavelength_cont, straight_fit)
    # plt.show()

    SNR = np.sqrt(np.abs(chi_squared_gaussian - chi_squared_con))
    #print(f"CHI SQUARED LINE: {SNR}")
    alt_SNR = max(gaussian(wavelength_gaussian, *popt))/np.nanmedian(np.array(counts_cont))
    print("-------")
    print(f"Alternative snr (corresponds to snr below): {alt_SNR}")
    return SNR

def SNR_lines():
    SNR_values = []
    for i in range(len(line_name)-1):
        SNR_value = SNR(counts_wavelength,line_wavelength_shift[i], counts, std)
        print(f"line name: {line_name[i]}, SNR: {SNR_value}")
        SNR_values.append(SNR_value)
    return SNR_values

n_files = len(line_name) - 1
n_rows = (n_files + 1) // 2
fig, axes = plt.subplots(n_rows, 2, figsize=(10,3 * n_rows))
SNR_values = SNR_lines()
mpl.rcParams['text.usetex'] = False
for i in range(len(line_name)-1):
    print(f"line wavelength under analysis: {line_wavelength_shift[i]}")
    mask = (counts_wavelength < line_wavelength_shift[i] + 70) & (counts_wavelength > line_wavelength_shift[i] - 70)
    wave_line = counts_wavelength[mask]
    counts_line = counts[mask]
    std_line = std[mask]

    row = i % n_rows
    col = i // n_rows

    center = line_wavelength_shift[i]
    sigma_guess = 12

    p0 = [np.nanmax(counts_line) - np.nanmedian(counts_line), center, sigma_guess]
    bounds = ([0, center-15, 10], [np.inf, center+15, 13])
    popt, pcov, infodict, mesg, ier = curve_fit(gaussian, wave_line,counts_line,sigma = std_line,absolute_sigma = True, p0=p0, bounds=bounds, full_output = True)
    y = gaussian(wave_line,*popt)
    axes[row,col].plot(wave_line, y, ls='--', c='red')
    axes[row, col].plot(wave_line, counts_line)
    axes[row, col].text(line_wavelength_shift[i] - 70, max(counts_line) - 0.2*max(counts_line), rf"{line_name[i]}, SNR = {SNR_values[i]:.1f}", bbox=dict(edgecolor='black', fc = 'white'))
    
plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/gaussian_resolved_plots/output_lines.png")
plt.show()


    
    


