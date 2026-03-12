from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table, fits, inset_axes, mark_inset,patches, curve_fit, GridSpec, axes3d, CircularAperture,aperture_photometry

import plotting_params
import matplotlib as mpl

#SETUP

output_dir = '/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/salpeter_run/salpeter_run/run_9'

counts_file = f'{output_dir}/run_9_reduced.fits'
flux_file = f'{output_dir}/run_9_reduced_flux_cal.fits'
std_file = f'{output_dir}/run_9_std.fits'
dir_save = '/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/salpeter_run/salpeter_run/run_9/SNR_values.npy'

extraction_region = 36
centre = 100
d = 14
xcoord = 114
ycoord =114
mass = 6 #solar mass
output_scale = 7*10**(-3)



def import_data(counts_file, flux_file, std_file):
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
    data_std = hdul[0].data
    std_data = np.sqrt(np.nansum(data_std**2, axis=(1,2)))
    header = hdul[0].header
    crval3 = header.get("CRVAL3", 0)
    cdelt3 = header.get("CDELT3", 1) 
    crpix3 = header.get("CRPIX3", 1)
    n_wave = data.shape[0]
    wavelength = crval3 + cdelt3*(np.arange(n_wave) - crpix3)
    std_wavelength = wavelength*(10**4)
    return counts_data, counts_wavelength, flux_data, flux_wavelength, std_data, std_wavelength

## APERTURE


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

def save_galaxy_spectra(counts_file, std_file, output_scale):
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
total_counts = []
#SNR_array = []

def collapsed_image(counts_file):
    hdul = fits.open(counts_file)
    data = hdul[0].data
    hdul_std = fits.open(std_file)
    std_data = hdul_std[0].data

    collapsed_image = np.nansum(data, axis=0)
    collapsed_image = np.nan_to_num(collapsed_image)
    ny, nx = collapsed_image.shape
    pixel_scale_mas = 7
    extent = [0, nx * pixel_scale_mas, 0, ny * pixel_scale_mas]
    # fig,ax = plt.subplots()
    # image = ax.imshow(collapsed_image, origin='lower', cmap='gray', aspect='equal', extent=extent)
    # circ = patches.Circle((101,101),36, alpha=0.8, fc='None',ec='red')
    # ax.add_patch(circ)
    # ax.set_xlabel('$x\ (mas)$')
    # ax.set_ylabel('$y\ (mas)$')
    # #plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/salpeter_run_5_5/M_5_5/mid/full_gauss_source_collapsed_image.png')
    # plt.show()
    

def plot_spectra(counts_file, flux_file, std_file):
    counts_data, counts_wavelength, flux_data, flux_wavelength, std_data, std_wavelength = import_data(counts_file, flux_file, std_file)

    # fig, (ax1,ax2) = plt.subplots(2,1, height_ratios = [1,0.3], sharex=True)
    # ax1.plot(flux_wavelength, flux_data)
    # ax2.set_xlabel("$Wavelength\ (\mathring{A})$")
    # ax1.set_ylabel("$Flux\ (erg \cdot cm^{-2} \cdot s^{-1} \cdot \mathring{A}^{-1})$")
    # ax1.tick_params(labelbottom=False)
    # ax2.plot(opacity_data["wavelength"], opacity_data["transmittance"], c='black')
    # ax2.set_ylabel("Transmission")
    # ax2.set_xlim(min(counts_wavelength), max(counts_wavelength))
    # #plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/salpeter_run_5_5/M_5_5/mid/full_output_spectra_raw_flux.png")
    # plt.show()

    # fig, (ax1,ax2) = plt.subplots(2,1, height_ratios = [1,0.3], sharex=True)
    # ax1.plot(counts_wavelength, counts_data)
    # ax1.tick_params(labelbottom=False)
    # ax2.set_xlabel("$Wavelength\ (\mathring{A})$")
    # ax1.set_ylabel("$Counts$")
    # ax2.plot(opacity_data["wavelength"], opacity_data["transmittance"], c='black')
    # ax2.set_ylabel("Transmission")
    # ax2.set_xlim(min(counts_wavelength), max(counts_wavelength))
    # #plt.savefig("/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/salpeter_run_5_5/M_5_5/mid/full_output_spectra_raw_counts.png")
    # plt.show()

### isolate popIII



def extract_region(data, region_size, centre, d,subtract_background=True):
    nlam, ny, nx = data.shape
    cy, cx = 114/7, 114/7#(centre + d) / 7, (centre + d) / 7
    #cy, cx = ycoord // 7, xcoord // 7
    half = region_size // 2

    y1, y2 = cy - half, cy + half
    x1, x2 = cx - half, cx + half

    #source_region = data[:, y1:y2, x1:x2]

    positions = [(cx,cy)]
    aperture = CircularAperture(positions, r= region_size // 2)

    counts = np.zeros(nlam)
    #error = np.zeros(nlam)
    for i in range(nlam):
        frame = data[i]     
        phot = aperture_photometry(frame, aperture,  method='exact')
        counts[i] = phot['aperture_sum'][0]
        #error[i] = phot['aperture_sum_err'][0]
    return counts#, error

def extract_region_counts(data, region_size,data_std,centre,d, subtract_background=True):
    nlam, ny, nx = data.shape
    cy, cx = 114/7, 114/7#(centre + d) / 7, (centre + d) / 7
    #cy, cx = ycoord // 7, xcoord // 7
    half = region_size // 2

    y1, y2 = cy - half, cy + half
    x1, x2 = cx - half, cx + half

    #source_region = data[:, y1:y2, x1:x2]

    positions = [(cx,cy)]
    aperture = CircularAperture(positions, r= region_size // 2)

    counts = np.zeros(nlam)
    error = np.zeros(nlam)
    for i in range(nlam):
        frame = data[i]     
        phot = aperture_photometry(frame, aperture,error = data_std[i],  method='exact')
        counts[i] = phot['aperture_sum'][0]
        error[i] = phot['aperture_sum_err'][0]
    return counts, error

def extract_region_2(data, region_size,centre, d, subtract_background=True):
    nlam, ny, nx = data.shape
    print(f"lam:{nlam}, ny:{ny}, nx:{nx}")
    cy, cx = (centre - d) / 7, (centre - d) /7
    #cy, cx = (200-ycoord) // 7, (200-xcoord) // 7
    half = region_size // 2

    y1, y2 = cy - half, cy + half
    x1, x2 = cx - half, cx + half

    #source_region = data[:, y1:y2, x1:x2]

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

def extract_region_counts_2(data, region_size,data_std, centre, d, subtract_background=True):
    nlam, ny, nx = data.shape
    print(f"lam:{nlam}, ny:{ny}, nx:{nx}")
    cy, cx = (centre - d) / 7, (centre - d) / 7
    #cy, cx = (200-ycoord) // 7, (200-xcoord) // 7
    half = region_size // 2

    y1, y2 = cy - half, cy + half
    x1, x2 = cx - half, cx + half

    #source_region = data[:, y1:y2, x1:x2]

    positions = [(cx,cy)]
    aperture = CircularAperture(positions, r= region_size // 2)

    counts = np.zeros(nlam)
    error = np.zeros(nlam)
    for i in range(nlam):
        frame = data[i]     
        phot = aperture_photometry(frame, aperture,error = data_std[i], method='exact')
        counts[i] = phot['aperture_sum'][0]
        error[i] = phot['aperture_sum_err'][0]
    return counts, error


def isolate_popIII(flux_file, counts_file, centre, d, output_scale, std_file, region_size):

    hdul = fits.open(flux_file)
    data = hdul[0].data * (output_scale**2) * (10**(-4))
    flux_region = extract_region(data, region_size,centre, d, subtract_background = True)

    flux_nearby = extract_region_2(data, region_size, centre, d,subtract_background = True)

    popIII_flux = flux_region - flux_nearby

    hdul_counts = fits.open(counts_file)
    data_counts = hdul_counts[0].data
    hdul_std = fits.open(std_file)
    data_std = hdul_std[0].data
    counts_region, std_region = extract_region_counts(data_counts, region_size,data_std,centre, d, subtract_background = True)
    counts_nearby, std_nearby = extract_region_counts_2(data_counts, region_size,data_std,centre, d, subtract_background = True)

    popIII_counts = counts_region - counts_nearby
    popIII_std = np.sqrt(std_region**2 + std_nearby**2)



    plt.figure()
    plt.plot(counts_wavelength, np.log10(flux_region), label='PopIII cluster region')
    plt.plot(counts_wavelength, np.log10(flux_nearby), c= 'red', alpha=0.6, label='isolated galaxy region')
    #plt.title("popIII area")
    plt.xlabel("$Wavelength\ (\mathring{A})$")
    plt.ylabel("$log\ Flux\ (erg \cdot cm^{-2} \cdot s^{-1} \cdot \mathring{A}^{-1})$")
    plt.legend()
    #plt.savefig('/home/steff/hsim/zackrisson_pop3_all/code/Report_plots/interim_report/GN_z11_plots/salpeter_run_5_5/M_5_5/mid/popIII_nopopIII_log.png')
    plt.show()

    plt.figure()
    plt.plot(counts_wavelength, flux_nearby)
    plt.title("general area")
    plt.show()

    plt.figure()
    plt.plot(counts_wavelength, popIII_flux)
    #plt.ylim(-0.2*10**(-19), 0.2*10**(-19))
    plt.show()
    return flux_region, flux_nearby, popIII_flux, counts_region, counts_nearby, popIII_counts, std_region, std_nearby, popIII_std


def gaussian(x, A, B, C):
    gaussian = (A/(np.sqrt(2*np.pi)*B))*np.exp(-0.5*((x-C)**2)/B**2)
    return gaussian

def chi_squared_gauss(x,y,std, popt):
    chi = ((y - gaussian(x, *popt))**2)/(std**2)
    return np.sum(chi)

def chi_squared_cont(x,y,std, straight_fit):
    chi = ((y - straight_fit)**2) / (std**2)
    #print(chi)
    return np.nansum(chi)

def analyse_lines(flux_region, flux_nearby, popIII_flux, counts_region, counts_nearby, popIII_counts, mass, dir_save, dir_save_image):
    flux_region = flux_region*10**(20)
    flux_nearby = flux_nearby*10**(20)
    popIII_flux = popIII_flux*10**(20)

    lines_dir = ascii.read('/home/steff/hsim/zackrisson_pop3_all/code/GN-z11_codes/lines.txt')
    line_names = lines_dir['line']
    line_wavelength = (np.array(lines_dir['wavelength']) / (1+10.6)) * (1+11.7)

    spectral_lines_file = ascii.read('/home/steff/hsim/spectra/spectra/GN-z11/lines')
    line_names = spectral_lines_file['line']

    line_wavelength = spectral_lines_file['l0']
    line_wavelength = ((line_wavelength)/(1+10.6))*(1+11.7)

    remove_idx = [2, 4, 8]

    line_names = np.delete(line_names, remove_idx)
    line_wavelength = np.delete(line_wavelength, remove_idx)
    line_wavelength[0] -= 45
    print(line_wavelength)

    n_files = len(line_names)
    n_rows = n_files
    n_cols = 3
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8,6))
    mpl.rcParams['text.usetex'] = False
    axes[0,0].text(0.25,1.15, "Cluster spectrum",transform=axes[0,0].transAxes)
    axes[0,1].text(0.25,1.15, "Galaxy spectrum",transform=axes[0,1].transAxes)
    axes[0,2].text(0.25,1.15, "Subtracted",transform=axes[0,2].transAxes)
    axes[3,0].text(-0.38, -1.8, "$Flux (erg \cdot cm^{-2} \cdot s^{-1} \cdot \mathring{A}^{-1}) \cdot 10^{20}$", transform=axes[3,0].transAxes, rotation='vertical')
    fig.subplots_adjust(hspace=0.6)
    axes[5,1].text(0,-1.1, "$Wavelength (\mathring{A})$", transform=axes[5,1].transAxes)

    SNR_values = []

    for i in range(len(line_names)):
        axes[i,0].plot(counts_wavelength, flux_region)
        axes[i,0].set_xlim(line_wavelength[i]  - 50, line_wavelength[i] + 50)
        axes[i,1].plot(counts_wavelength, flux_nearby)
        axes[i,1].set_xlim(line_wavelength[i]  - 50, line_wavelength[i] + 50)
        axes[i,2].plot(counts_wavelength, popIII_flux)
        axes[i,2].set_xlim(line_wavelength[i]  - 50, line_wavelength[i] + 50)
        axes[i,0].text(-0.42, 0.5,line_names[i],transform=axes[i,0].transAxes,ha="right",va="center")
        mask_lines = (counts_wavelength < line_wavelength[i] + 70) & (counts_wavelength > line_wavelength[i] - 70)
        line_lambda = counts_wavelength[mask_lines]
        line_flux_popIII_coord = flux_region[mask_lines]
        line_counts_popIII_coord = counts_region[mask_lines]
        line_std_popIII_coord = std_region[mask_lines]
        line_flux_nearby = flux_nearby[mask_lines]
        line_counts_nearby = counts_nearby[mask_lines]
        line_std_nearby = std_nearby[mask_lines]
        line_flux_popIII = popIII_flux[mask_lines]
        line_counts_popIII = popIII_counts[mask_lines]
        line_std_popIII = popIII_std[mask_lines]
        center = line_wavelength[i]
        cluster_mass = (10**(mass))*(1.99*10**30)
        G = 6.67*10**(-11)
        r = 4 * (3.086*10**16)
        sigma_gal = np.sqrt(cluster_mass * G / r)
        c = 3*10**8 #m/s
        sigma_line = (sigma_gal/c) * 1640*(1+11.7)
        sigma_guess = sigma_line
        p0 = [np.nanmax(line_flux_popIII) - np.nanmedian(line_flux_popIII), sigma_guess, center]  # [amp, mu, sigma]
        bounds = ([0,sigma_guess - 0.5, center-0.5], [np.inf, sigma_guess + 0.5, center + 0.5])
        popt, pcov, infodict, mesg, ier = curve_fit(gaussian, line_lambda,line_flux_popIII, p0=p0, bounds=bounds, full_output = True)
        axes[i,2].plot(line_lambda, gaussian(line_lambda, *popt), c='red', ls='--')

        p0 = [np.nanmax(line_counts_popIII) - np.nanmedian(line_counts_popIII), sigma_guess, center]  # [amp, mu, sigma]
        bounds = ([0,sigma_guess - 0.5, center-0.5], [np.inf, sigma_guess + 0.5, center + 0.5])
        popt, pcov, infodict, mesg, ier = curve_fit(gaussian, line_lambda,line_counts_popIII, p0=p0, bounds=bounds, full_output = True)

        line_lambda = np.array(line_lambda)
        popt_0 = np.array(popt[2])
        popt_1 = np.array(popt[1])
        peak = popt[0] / (np.sqrt(2*np.pi)*(popt[1]))
        mask_outside = (line_lambda < (popt_0-4*popt_1)) | (line_lambda > (popt_0 + 4*popt_1))
        line_cont = line_lambda[mask_outside]
        flux_cont = line_flux_popIII[mask_outside]
        counts_cont = line_counts_popIII[mask_outside]
        std_cont = line_std_popIII[mask_outside]
        median_cont = np.nanmedian(counts_cont)
        residuals = counts_cont - median_cont
        std_cont = np.std(counts_cont)
        straight_fit = np.ones(len(line_cont))*median_cont
        line = gaussian(line_lambda, *popt)
        # chi_line = chi_squared_gauss(line_lambda, line_counts_popIII,line_std_popIII, popt)
        # print(f"chi line: {chi_line}")
        # chi_cont = chi_squared_cont(line_lambda,line_counts_popIII,line_std_popIII ,straight_fit)
        # print(f"chi cont: {chi_cont}")
        SNR = peak / std_cont #np.sqrt(abs(chi_line - chi_cont))
        print(f"SNR: {SNR}")
        SNR_values.append(SNR)
        # plt.figure()
        # plt.plot(line_lambda, line_counts_popIII)
        # plt.plot(line_lambda, gaussian(line_lambda, *popt), ls='--', c='green')
        # plt.scatter(line_cont, straight_fit)
        # plt.show()
        axes[i, 2].text(1.03, 0.5, f"SNR:{SNR:.1f}", transform=axes[i,2].transAxes)
        if line_flux_popIII_coord.size > 0:
            axes[i,0].set_ylim(
            np.nanmin(line_flux_popIII) - 0.15*np.nanmedian(line_flux_popIII),
                np.nanmax(line_flux_popIII_coord) + 0.15*np.nanmedian(line_flux_popIII)
            )
            axes[i,1].set_ylim(
                np.nanmin(line_flux_popIII) - 0.15*np.nanmedian(line_flux_popIII),
                np.nanmax(line_flux_popIII_coord) + 0.15*np.nanmedian(line_flux_popIII)
            )
            axes[i,2].set_ylim(
                np.nanmin(line_flux_popIII) - 0.15*np.nanmedian(line_flux_popIII),
                np.nanmax(line_flux_popIII_coord)+ 0.15*np.nanmedian(line_flux_popIII)
            )

    for i in range(n_rows):
        for j in range(n_cols):
            ax = axes[i, j]

            # Hide y-axis numbers except first column
            if j != 0:
                ax.tick_params(labelleft=False)
    plt.savefig(dir_save_image)
    plt.show()

    SNR_values = np.array(SNR_values)

    np.save(f"{dir_save}", SNR_values)
    

counts_data, counts_wavelength, flux_data, flux_wavelength, std_data, std_wavelength = import_data(counts_file, flux_file, std_file)

save_galaxy_spectra(counts_file, std_file, output_scale)
collapsed_image(counts_file)
plot_spectra(counts_file, flux_file, std_file)

region_size = 2

flux_region, flux_nearby, popIII_flux, counts_region, counts_nearby, popIII_counts, std_region, std_nearby, popIII_std = isolate_popIII(flux_file, counts_file, centre, d, output_scale, std_file, region_size)

mass = 7
dir_save_image = f'/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/SNR_arrays_move_cen_back/spectra.png'
analyse_lines(flux_region, flux_nearby, popIII_flux, counts_region, counts_nearby, popIII_counts, mass, dir_save, dir_save_image)

