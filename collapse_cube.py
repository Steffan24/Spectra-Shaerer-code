
from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table, fits, inset_axes, mark_inset,patches, curve_fit, GridSpec, axes3d, CircularAperture,aperture_photometry
from constants import T_sun, c_m, T_100M, M_sun_kg, G, kb, c, h, pc, AU, d, R_sun, M_sun
import plotting_params
from variables import output_file, dir_basic, output_mass, corresponding_V, output_exp_time, output_flux, output_std, output_scale, input_file, input_scale


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

    print(counts)
    return counts, error

   
    
    
def collapse_cube_data(output_file, output_flux, output_std, output_scale):
    #COUNTS
    fits_image_filename = output_file
    hdul = fits.open(fits_image_filename)
    data = hdul[0].data
    data_spectrum = np.nansum(data, axis=(1,2))
    header = hdul[0].header
    crval3 = header.get("CRVAL3", 0)
    cdelt3 = header.get("CDELT3", 1) 
    crpix3 = header.get("CRPIX3", 1)
    n_wave = data.shape[0]
    wavelength = crval3 + cdelt3*(np.arange(n_wave) - crpix3)
    wavelength_angstrom = wavelength*(10**4)

    #FLUX
    fits_image_filename = output_flux
    hdul_flux = fits.open(fits_image_filename)
    data_flux = hdul_flux[0].data * (output_scale**2) * (10**(-4))
    flux_data = np.nansum(data_flux, axis=(1,2))
    header = hdul_flux[0].header
    crval3 = header.get("CRVAL3", 0)
    cdelt3 = header.get("CDELT3", 1) 
    crpix3 = header.get("CRPIX3", 1)
    n_wave = data_flux.shape[0]
    wavelength_flux = crval3 + cdelt3*(np.arange(n_wave) - crpix3)
    wavelength_angstrom_flux = wavelength_flux*(10**4)

    #STD
    fits_image_filename = output_std
    hdul_std = fits.open(fits_image_filename)
    data_std = hdul_std[0].data
    std_data = np.sqrt(np.nansum(data_std**2, axis=(1,2)))
    header = hdul_std[0].header
    crval3 = header.get("CRVAL3", 0)
    cdelt3 = header.get("CDELT3", 1) 
    crpix3 = header.get("CRPIX3", 1)
    n_wave = data_std.shape[0]
     #* (output_scale**2)
    wavelength_std = crval3 + cdelt3*(np.arange(n_wave) - crpix3)
    wavelength_angstrom_std = wavelength_flux*(10**4)

    return wavelength_angstrom, data_spectrum,data, wavelength_angstrom_flux,flux_data,data_flux, wavelength_angstrom_std, std_data, data_std

                       
def extracting_over_aperture(wavelength_angstrom, data_spectrum,data, wavelength_angstrom_flux,flux_data,data_flux, wavelength_angstrom_std, std_data, data_std, dir_basic, output_mass, output_exp_time):
    fig = plt.figure(111)#, (ax1,ax3) = plt.subplots(2,1, height_ratios = [2,1])
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 0.6])

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    array = np.linspace(1,11,11, dtype = int)
    SNR_totals = []
    total_counts = []
    std_totals = []
    wavelength_helium_mask = (wavelength_angstrom > 20775) & (wavelength_angstrom < 20825)
    wavelength_helium = wavelength_angstrom[wavelength_helium_mask]
    data_helium = data[wavelength_helium_mask]
    std_helium = data_std[wavelength_helium_mask]
    for i in range(len(array)):
        data_central, central_std = extract_central_region(data_helium, array[i]*2,std_helium, True)
        central_std_collapsed = np.sqrt(np.nansum(central_std**2))#, axis=(1,2)))
        data_central = data_central
        data_central_collapsed = np.nansum(data_central)#, axis=(1,2))
        SNR = data_central_collapsed / central_std_collapsed
        total_counts.append(data_central_collapsed)
        SNR_totals.append(SNR)
        std_totals.append(central_std_collapsed)


    pixel_to_arcsec = 7*10**(-3)

    arcsecs = (array * pixel_to_arcsec)
    
    ax1.plot(array, np.array(total_counts) * (10**-4))
    ax1.set_ylabel("\(Counts (\cdot 10^4)\)")
    ax1.set_xlabel("\(Aperture\ radius\ (Pixels)\)")
    ax1.set_xlim(1,max(array))
    ax1.vlines(x=3,ymin=0, ymax = np.max(total_counts)*(10**-4) , color = 'red', ls = '--')
    ax4 = ax1.twiny()
    ax4.set_xlim(min(arcsecs)*100, max(arcsecs)*100)
    ax4.set_xlabel("\(Aperture\ radius\ \cdot\ 10^{{-2}}\ (arcsec)\)")

    #ax2 = ax1.twinx()
    

    SNR_totals = np.array(SNR_totals) 

    ax2.plot(array, np.array(std_totals) * (10**(-4)))
    ax2.set_ylabel("\(\sigma_{counts} (\cdot 10^4)\)")
    ax2.set_xlabel("\(Aperture\ radius\ (Pixels)\)")
    ax2.set_xlim(1, max(array))
    ax2.vlines(x=3,ymin=0, ymax = np.max(std_totals)*10**(-4) , color = 'red', ls = '--')
    ax3.vlines(x=3,ymin=np.min(SNR_totals), ymax = np.max(SNR_totals) , color = 'red', ls = '--', label = '\(R_{aperture} = 3\ pixels\)')
    ax3.legend()
    ax5 = ax2.twiny()
    ax5.set_xlim(min(arcsecs)*100, max(arcsecs)*100)
    ax5.set_xlabel("\(Aperture\ radius\ \cdot\ 10^{{-2}}\ (arcsec)\)")
    

    ax3.plot(array, SNR_totals)
    ax3.set_xlabel("\(Aperture\ radius\ (Pixels)\)")
    ax3.set_ylabel("\(SNR\)")
    ax3.set_xlim(1, np.max(array))
    ax6 = ax3.twiny()
    ax6.set_xlim(min(arcsecs)*100, max(arcsecs)*100)
    ax6.set_xlabel("\(Aperture\ radius\ \cdot\ 10^{{-2}}\ (arcsec)\)")

    plt.subplots_adjust(hspace=1, wspace=0.6)
    plt.show()

    
    extraction_region = 6
                       
    #COUNTS
    data_central, data_central_std =  extract_central_region(data, extraction_region,data_std, subtract_background=True)
    spectrum = data_central #np.nansum(data_central)#, axis=(1,2))

    #FLUX
    data_central_flux, flux_std_not_used = extract_central_region(data_flux, extraction_region,data_std, subtract_background=True)
    spectrum_flux =data_central_flux # np.nansum(data_central_flux)#, axis=(1,2))

    #STD
    #data_central_std =   extract_central_region(data_std, extraction_region, subtract_background=True)
    spectrum_std = data_central_std #np.sqrt(np.nansum(data_central_std**2))#, axis=(1,2)))

    np.save(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/spectrum_counts", spectrum)
    np.save(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/spectrum_flux", spectrum_flux)
    np.save(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/spectrum_std", spectrum_std)
    np.save(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/wavelength_angstrom", wavelength_angstrom)
    np.save(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/data", data)
    np.save(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/data_flux", data_flux)
    np.save(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/data_std", data_std)
    print("saved file")
    plt.figure()
    plt.plot(wavelength_angstrom, spectrum)
    plt.show()

    return spectrum, spectrum_flux, spectrum_std



wavelength_angstrom, data_spectrum,data, wavelength_angstrom_flux,flux_data,data_flux, wavelength_angstrom_std, std_data, data_std = collapse_cube_data(output_file, output_flux, output_std, output_scale)

spectrum, spectrum_flux, spectrum_std = extracting_over_aperture(wavelength_angstrom, data_spectrum,data, wavelength_angstrom_flux,flux_data,data_flux, wavelength_angstrom_std, std_data, data_std, dir_basic, output_mass, output_exp_time)

#collapse_input_cube(input_file, input_scale)

