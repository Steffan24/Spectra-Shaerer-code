# functions.py
# list of functions used

from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table, fits, inset_axes, mark_inset,patches, curve_fit, GridSpec, axes3d, CircularAperture,aperture_photometry
from constants import T_sun, c_m, T_100M, M_sun_kg, G, kb, c, h, pc, AU, d, R_sun, M_sun
import plotting_params
#from SNR_calculation_chi import calcSNR_chi2


 
def spectral_analysis(spectrum, spectrum_flux, spectrum_std, wavelength_angstrom, data, data_std):

    #EW Calculation
    mask = (wavelength_angstrom > 23000) & (wavelength_angstrom < 23300)
    wavelength_zoomed_flux = wavelength_angstrom[mask]
    counts_zoomed = spectrum[mask]
    std_zoomed = spectrum_std[mask]
    flux_zoomed = spectrum_flux[mask]
    # plt.figure()
    # plt.plot(wavelength_zoomed_flux, counts_zoomed, label='\(Observed\ flux\)')
    # plt.xlabel("\(wavelength (\mathring{A})\)")
    # plt.ylabel("\(Flux (\cdot\ 10^{-10})\)")

    mask_gauss =  (wavelength_zoomed_flux > 23190) & (wavelength_zoomed_flux < 23220)
    wavelength_gauss = wavelength_zoomed_flux#[mask_gauss]
    counts_gauss = counts_zoomed#[mask_gauss]
    std_gauss = std_zoomed#[mask_gauss]
    flux_gauss = flux_zoomed    
    def gaussian(x,A):
        std = 3.244
        B = 23205
        gaussian = (A/(np.sqrt(2*np.pi)*std))*np.exp(-0.5*((x-B)**2)/std**2)
        return gaussian
    popt, pcov = curve_fit(gaussian, wavelength_gauss, counts_gauss, p0 = [1.5*10**5]) #sigma=std_gauss, absolute_sigma=True)
    #print(f"{popt[0]:.2f} exp(-0.5 ((x-{popt[1]:.2f})^2)/{popt[2]:.2f}")
    FWHM_w = 2*np.sqrt(2*np.log(2))*3.244
    plt.plot(wavelength_gauss, gaussian(wavelength_gauss,*popt),ls = '--', c='red', label=f'$Gaussian\ fit: FWHM = {FWHM_w:.2f}\mathring{{A}}$')
    plt.legend(loc='upper right')
    #plt.show()


    velocity_dispersion = c_m * (FWHM_w/2.355)/23204
    print(f"volcity dispersion: {velocity_dispersion}m/s")
    chi_squared_lin = ((counts_gauss - gaussian(wavelength_gauss, *popt))**2 / (std_gauss)**2)
    chi_squared_line = np.sum(chi_squared_lin)
    print(f"difference: {gaussian(wavelength_gauss, *popt) - counts_gauss}")
    print(f"std: {std_gauss}")
    print(f"chi squared gauss: {chi_squared_line}")

    #flux of line
    popt_flux, pcob_flux = curve_fit(gaussian, wavelength_gauss, flux_gauss, p0=[4*10**(-14)])
    peak = popt_flux[0]
    # plt.figure()
    # plt.plot(wavelength_gauss, flux_gauss)
    # plt.plot(wavelength_gauss, gaussian(wavelength_gauss, *popt_flux), label='gauss')
    # plt.ylabel("flux")
    # plt.legend()
    # plt.xlabel
    ("wavelength")
    plt.show()
                                     
    width = 3.244
    total_flux_gaussian = peak

    def straight_fit(x, A):
        y = A
        return y
    
    mask_2 = ((wavelength_zoomed_flux > 23215))# & (wavelength_zoomed_flux < 23300))
    mask_3 = ((wavelength_zoomed_flux < 23190))# & (wavelength_zoomed_flux > 23000))#17707))
    wavelength_zoomed_1 = wavelength_zoomed_flux[mask_2]
    wavelength_zoomed_2 = wavelength_zoomed_flux[mask_3]
    zoomed_flux_1 = counts_zoomed[mask_2]
    zoomed_flux_2 = counts_zoomed[mask_3]
    std_zoom_1 = std_zoomed[mask_2]
    std_zoom_2 = std_zoomed[mask_3]
    wavelength_line_fit = wavelength_zoomed_flux#np.concatenate([wavelength_zoomed_1, wavelength_zoomed_2])
    flux_line_fit = counts_zoomed#np.concatenate([zoomed_flux_1, zoomed_flux_2])
    std_line_fit = std_zoomed#np.concatenate([std_zoom_1, std_zoom_2])
    popt, pcov = curve_fit(straight_fit, wavelength_line_fit, flux_line_fit)
    # plt.figure()
    # plt.plot(wavelength_line_fit, flux_line_fit)
    # plt.hlines(straight_fit(wavelength_zoomed_flux, *popt), xmin=0, xmax=1)
    # plt.title("chi squared continuum range")
    # plt.show()
    print(f"A: {popt[0]}")
    straight_fit_array = straight_fit(wavelength_line_fit, *popt) * np.ones(len(flux_line_fit))
    chi_squared_continuu = ((flux_line_fit - straight_fit_array)**2 / (std_line_fit)**2)
    chi_squared_continuum = np.sum(chi_squared_continuu)
    print(f"diffence cont: {straight_fit_array - flux_line_fit}")
    print(f"std cont: {std_line_fit}")
    print(f"chi cont: {chi_squared_continuum}")
    delta_chi = abs(chi_squared_continuum - chi_squared_line)
    print(delta_chi)
    SNR_He = np.sqrt(delta_chi)
    print(f"SIGNAL TO NOISE: {SNR_He}")

    


    #SUBPLOT 1:
    extraction_region = 3
    # fig = plt.figure(figsize=(8, 6))
    # gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 0.6])

    # ax1 = fig.add_subplot(gs[0, 0])
    # ax2 = fig.add_subplot(gs[0, 1])
    # ax3 = fig.add_subplot(gs[1, :])

    # wavelength_cut = wavelength_angstrom[mask]
    # signal_cut, noise_cut = extract_central_region(data, extraction_region,data_std, True)
    # signal_cut = signal_cut[mask]
    # noise_cut = noise_cut[mask]
    # #signal_cut = np.nansum(signal_cut)[mask]#, axis=(1,2))[mask]
    # ax1.plot(wavelength_cut, signal_cut)
    # ax1.set_xlabel("\(Wavelength\ (\mathring{A})\)")
    # ax1.set_ylabel("\(Counts\)")
    
    # #noise_cut = extract_central_region(data_std, extraction_region, True)[mask]
    # #noise_cut = np.sqrt(np.nansum(noise_cut**2))[mask]#, axis=(1,2)))[mask]
    # ax2.plot(wavelength_cut, noise_cut)
    # ax2.set_xlabel("\(Wavelength\ (\mathring{A})\)")
    # ax2.set_ylabel("\(\sigma_{counts}\)")

    # SNR_cut = signal_cut / noise_cut
    # ax3.plot(wavelength_cut, SNR_cut)
    # ax3.set_xlabel("\(Wavelength\ (\mathring{A})\)")
    # ax3.set_ylabel("\(SNR\)")
    # plt.subplots_adjust(hspace=0.4, wspace=0.3)
    # plt.show()

    return SNR_He, total_flux_gaussian

     

    

def collapse_all(output_array, median_magnitudes,std_array, V_array, SNR_FULL):
    fig, axes = plt.subplots(3, 2, height_ratios=[1,1,1], sharex=True, sharey = True)
    axes = axes.flatten()
    for i in range(len(output_array)):
        fits_image_filename = output_array[i]
        fits_std_filename = std_array[i]
        with fits.open(fits_image_filename) as hdul:
            data = hdul[0].data
            header = hdul[0].header
        with fits.open(fits_std_filename) as hstd:
            std_data = hstd[0].data

        data_central, std_central = extract_central_region(data, 6,std_data, subtract_background=True)
        spectrum = data_central#np.nansum(data_central)#, axis=(1, 2))
        crval3 = header.get("CRVAL3", 0)
        cdelt3 = header.get("CDELT3", 1) 
        crpix3 = header.get("CRPIX3", 1)
        n_wave = data.shape[0]
        wavelength = crval3 + cdelt3*(np.arange(n_wave) - crpix3)
        wavelength_angstrom = wavelength*(10**4)
        axes[i].plot(wavelength_angstrom,(spectrum*10**(-6)))
        #axes[i].set_ylim(32,15)
        axes[i].text(22000,2.1,fr'$M_{{stellar}} = 10^{median_magnitudes[i]:.1f}\ M_{{\odot}}$',bbox=dict(edgecolor='None', fc = 'None'))#, fontsize = 18)
    fig.supylabel('\(Counts (\cdot 10^{6})\)', x=0.05, rotation=90)
    plt.show()
    
    n_files = len(output_array)
    fig, axes = plt.subplots(n_files, 4, figsize=(10, 3 * n_files),
                             gridspec_kw={'width_ratios': [2, 1.5,1.2,1.2]})#, constrained_layout=True)
    if n_files == 1:
        axes = np.array([axes])  # ensure 2D array

    for i in range(len(output_array)):
        fits_image_filename = output_array[i]
        fits_std_filename = std_array[i]
        with fits.open(fits_image_filename) as hdul:
            data = np.squeeze(hdul[0].data)
            header = hdul[0].header
        with fits.open(fits_std_filename) as hstd:
            data_std = np.squeeze(hstd[0].data)
        

        # make sure data has shape (nz, ny, nx)
        if data.ndim == 2:
            data = data[np.newaxis, :, :]

        # --- SPECTRUM EXTRACTION ---
        data_central, std_central = extract_central_region(data, 6, data_std)
        spectrum = data_central#np.nansum(data_central)#, axis=(1, 2))

        crval3 = header.get("CRVAL3", 0)
        cdelt3 = header.get("CDELT3", 1)
        crpix3 = header.get("CRPIX3", 1)
        n_wave = data.shape[0]
        wavelength = crval3 + cdelt3 * (np.arange(n_wave) - (crpix3 - 1))
        wavelength_angstrom = wavelength * 1e4


        # --- COLLAPSED IMAGE ---
        collapsed_image = np.nansum(data, axis=0)  # integrate over all wavelengths (x,y)
        collapsed_image = np.nan_to_num(collapsed_image)
        ny, nx = collapsed_image.shape
        pixel_scale_mas = 7
        extent = [0, nx * pixel_scale_mas, 0, ny * pixel_scale_mas]

        im_ax = axes[i, 1]
        im = im_ax.imshow(collapsed_image, origin='lower', cmap='gray',
                          aspect='equal', extent = extent)  # <- equal keeps pixel aspect ratio
        #im_ax.set_title(f"Collapsed Image {i+1}", fontsize=10)
        # ny, nx = collapsed_image.shape
        # cy, cx = ny // 2, nx // 2
        # half = 6 // 2   # <-- same region_size as in extract_central_region (6)
        # y1, y2 = cy - half, cy + half
        # x1, x2 = cx - half, cx + half

        # # --- Add red rectangle around that region ---
        # rect = patches.Rectangle(
        #     (x1, y1),        # (x, y) of lower-left corner
        #     x2 - x1,         # width
        #     y2 - y1,         # height
        #     linewidth=2,
        #     edgecolor='red',
        #     facecolor='none'
        # )
        # #im_ax.add_patch(rect)
        im_ax.set_xlabel("\(X (mas)\)")
        #im_ax.set_ylabel("Y (mas)")
        im_ax.axis('on')

        # optional colorbar
        #cbar = fig.colorbar(im, ax=im_ax, fraction=0.046, pad=0.04)
        #cbar.ax.tick_params(labelsize=8)

        # --- SPECTRUM PLOT ---
        sp_ax = axes[i, 0]
        sp_ax.plot(wavelength_angstrom, spectrum*10**(-6))
        sp_ax.set_ylim(-5*10**(-3), np.max(spectrum*10**(-6))/2)
        #sp_ax.set_xlabel("Wavelength (Å)")
        #sp_ax.set_ylabel("Flux × 10²¹")
        #sp_ax.text(0.05, 0.9,
         #          fr'$M_{{\mathrm{{theory}}}} = {median_magnitudes[i]:.1f}$',
          #         transform=sp_ax.transAxes,
           #        bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
        if i == n_files - 1:
            sp_ax.set_xlabel(r"$Wavelength\ (\mathrm{\AA})$")
        else:
            sp_ax.set_xticklabels([])
        #sp_ax.set_ylabel(r"$Counts\ (\cdot 10^{6})$")

        zoom = axes[i,2]
        zoom.plot(wavelength_angstrom, spectrum*10**(-6))
        zoom_center = 23207
        zoom_halfwidth = 50 
        zoom.set_xlim(zoom_center - zoom_halfwidth, zoom_center + zoom_halfwidth)
        mask = (wavelength_angstrom > zoom_center - zoom_halfwidth) & (wavelength_angstrom < zoom_center + zoom_halfwidth)
        if np.any(mask):
            counts_region = spectrum[mask]*10**(-6)
            diff = counts_region.max() - counts_region.min()
            zoom.set_ylim(counts_region.min() - (diff/10), counts_region.max() + (diff/10))
        #zoom.yaxis.tick_right()
        #zoom.yaxis.set_label_position("right")
        if i == n_files - 1:
            zoom.set_xlabel(r"$Wavelength\ (\mathrm{\AA})$")
        else:
            zoom.set_xticklabels([])

        collapsed_image = np.nansum(data[mask,:,:], axis=0)  # integrate over all wavelengths (x,y)
        collapsed_image = np.nan_to_num(collapsed_image)

        im2_ax = axes[i, 3]
        ny, nx = collapsed_image.shape
        pixel_scale_mas = 7
        extent = [0, nx * pixel_scale_mas, 0, ny * pixel_scale_mas]
        im2 = im2_ax.imshow(collapsed_image, origin='lower', cmap='gray',
                            aspect='equal', extent=extent)  # <- equal keeps pixel aspect ratio
        #im_ax.set_title(f"Collapsed Image {i+1}", fontsize=10)

        constants = [7.5, 7, 6.5,6, 5.5]
        constants = np.array(constants)
        im2_ax.set_xlabel("\(X (mas)\)")
        im2_ax.axis('on')
        im2_ax.text(1.05, 0.5,f"$M_{{cluster}} = 10^{{{constants[i]:.2f}}} M_{{\\odot}}$",
                    transform=im2_ax.transAxes, ha='left', va='center', fontsize=30)

    #SNR_total_counts = calcSNR_chi2(wavelength_angstrom,spectrum, spectrum_std)
        
    fig.subplots_adjust(
    left=0.1, 
    right=0.85, 
    bottom=0.1,
    top=0.95, 
    wspace=0.25,
    hspace=0.35
    )
    fig.text(0.03, 0.5, r"$\mathrm{Counts\ (\times10^{6})}$", 
             ha='center', va='center', rotation='vertical', fontsize=30)
    plt.show()

    fig, axes = plt.subplots(n_files, 2, figsize=(10, 3 * n_files),
                             gridspec_kw={'width_ratios': [1.5, 1]})#, constrained_layout=True)

    for i in range(len(output_array)):
        fits_image_filename = output_array[i]
        fits_std_filename = std_array[i]
        with fits.open(fits_image_filename) as hdul:
            data = np.squeeze(hdul[0].data)
            header = hdul[0].header
        with fits.open(fits_std_filename) as hstd:
            data_std = np.squeeze(hstd[0].data)
        

        # make sure data has shape (nz, ny, nx)
        if data.ndim == 2:
            data = data[np.newaxis, :, :]

        # --- SPECTRUM EXTRACTION ---
        data_central, std_central = extract_central_region(data, 6, data_std)
        spectrum = data_central#np.nansum(data_central)#, axis=(1, 2))

        crval3 = header.get("CRVAL3", 0)
        cdelt3 = header.get("CDELT3", 1)
        crpix3 = header.get("CRPIX3", 1)
        n_wave = data.shape[0]
        wavelength = crval3 + cdelt3 * (np.arange(n_wave) - (crpix3 - 1))
        wavelength_angstrom = wavelength * 1e4
        
        SNR_val,total_counts = calcSNR_chi2(wavelength_angstrom, spectrum, std_central)

        zoom = axes[i,0]
        def gaussian(x,A):
            std = 3.244
            B = 23207.3
            gaussian = (A/(np.sqrt(2*np.pi)*std))*np.exp(-0.5*((x-B)**2)/std**2)
            return gaussian
        popt, pcov = curve_fit(gaussian, wavelength_angstrom, spectrum, p0 = [1.5*10**5])
        zoom.plot(wavelength_angstrom, spectrum*10**(-6))
        zoom.plot(wavelength_angstrom, gaussian(wavelength_angstrom, *popt)*10**(-6), ls='--', c='red')
        zoom_center = 23207
        zoom_halfwidth = 50 
        zoom.set_xlim(zoom_center - zoom_halfwidth, zoom_center + zoom_halfwidth)
        mask = (wavelength_angstrom > zoom_center - zoom_halfwidth) & (wavelength_angstrom < zoom_center + zoom_halfwidth)
        if np.any(mask):
            counts_region = spectrum[mask]*10**(-6)
            diff = counts_region.max() - counts_region.min()
            zoom.set_ylim(counts_region.min() - (diff/10), counts_region.max() + (diff/10))
        #zoom.yaxis.tick_right()
        #zoom.yaxis.set_label_position("right")
        if i == n_files - 1:
            zoom.set_xlabel(r"$Wavelength\ (\mathrm{\AA})$")
        else:
            zoom.set_xticklabels([])

        collapsed_image = np.nansum(data[mask,:,:], axis=0)  # integrate over all wavelengths (x,y)
        collapsed_image = np.nan_to_num(collapsed_image)

        im2_ax = axes[i, 1]
        ny, nx = collapsed_image.shape
        pixel_scale_mas = 7
        extent = [0, nx * pixel_scale_mas, 0, ny * pixel_scale_mas]
        im2 = im2_ax.imshow(collapsed_image, origin='lower', cmap='gray',
                            aspect='equal', extent=extent)  # <- equal keeps pixel aspect ratio
        #im_ax.set_title(f"Collapsed Image {i+1}", fontsize=10)

        constants = [7.5, 7, 6.5,6, 5.5]
        constants = np.array(constants)
        if i == n_files - 1:
            im2_ax.set_xlabel("\(X (mas)\)")
        else:
            zoom.set_xticklabels([])
        im2_ax.axis('on')
        im2_ax.text(1.2, 0.5,fr"$M_{{cluster}} = 10^{{{constants[i]:.2f}}} M_{{\odot}}$"
"\n"
fr"$V_{{median}} = {V_array[i]}$"
                    "\n"
                    fr"$SNR = {SNR_val:.2f}$",
                    transform=im2_ax.transAxes, ha='left', va='center', fontsize=30)
        
    fig.subplots_adjust(
    left=0.1, 
    right=0.85, 
    bottom=0.1,
    top=0.95, 
    wspace=0.15,
    hspace=0.35
    )
    fig.text(0.03, 0.5, r"$\mathrm{Counts\ (\times10^{6})}$", 
             ha='center', va='center', rotation='vertical', fontsize=30)
    plt.show()

    
        
    
    
    
    
    
    
    
    
    

