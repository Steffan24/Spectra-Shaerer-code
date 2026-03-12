from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table, fits, inset_axes, mark_inset,patches, curve_fit, GridSpec, axes3d, CircularAperture,aperture_photometry
from variables import output_file, output_flux, output_std, output_scale, dir_basic, output_mass, output_exp_time,  mass_output_array, input_scale, dir_basic_logA, dir_basic_logE, dir_basic_high_res
from constants import T_sun, c_m, T_100M, M_sun_kg, G, kb, c, h, pc, AU, d, R_sun, M_sun
import plotting_params


def spatial_plots(spectrum, spectrum_flux,data_flux, spectrum_std, wavelength_angstrom):
    image_collapsed = np.nansum(data_flux, axis=0)
    ny, nx = image_collapsed.shape
    cy, cx = ny // 2, nx // 2

    x_cut = image_collapsed[cy, :] #* 10**(13)
    y_cut = image_collapsed[:, cx] #* 10**(13)

    x_coords = np.arange(nx) - cx
    y_coords = np.arange(ny) - cy

    std_dev = np.std(x_cut)
    def gaussian_fit(x_coords, stand):
        gaussian = (max(x_cut)*np.exp(-0.5*((x_coords - 0)**2)/stand**2))
        return gaussian

    x_values = np.linspace(-15,15,100)
    x_values_as = x_values * (7*10**(-3))
    popt, pcov = curve_fit(gaussian_fit, x_coords, x_cut)
    #print(popt)
    
    
    lambda_he = 20800*10**(-10)
    D = 39
    FWHM = 1.22*(lambda_he/D)
    print(f"FWHM: {FWHM}")
    FWHM_pixels = FWHM*206265 #/ (7*10**(-3))
    gaussian_seeing = np.max(x_cut)*np.exp(-0.5*((x_values_as)**2)/(FWHM_pixels/(2*np.sqrt(2*np.log(2))))**2)

    pixel_to_arcsec = 7*10**(-3)

    arcsecs = (x_values * pixel_to_arcsec)
    popt_as = popt[0] * pixel_to_arcsec
    sigma_dl = (FWHM_pixels/(2*np.sqrt(2*np.log(2))))
    sigma_dl_as = sigma_dl * pixel_to_arcsec

    JWST_noise_FWHM = 0.065 #arcsec
    JWST_FWHM_pixels = JWST_noise_FWHM / (7*10**(-3))
    JWST_gaussian = np.max(x_cut)*np.exp(-0.5*((x_values_as)**2)/(JWST_noise_FWHM/(2*np.sqrt(2*np.log(2))))**2)
    sigma_JW = (JWST_noise_FWHM/(2*np.sqrt(2*np.log(2))))

    x_coords_as = x_coords * (7*10**(-3))
    
        
    fig, ax1 = plt.subplots()
    ax1.step(x_coords_as, x_cut*10**(18), label='$PSF_{observed}$')
    #ax1.plot(x_values, gaussian_fit(x_values, *popt)*10**(18), c='orange', label=f'\(Gaussian\ fit:\ \sigma\ =\)$ {popt[0]:.1f}$')
    ax1.plot(x_values_as, gaussian_seeing*10**(18), ls='--', c='green', label= '$PSF_{DL}$')
    ax1.plot(x_values_as, JWST_gaussian*10**(18), ls='--', c='red', label='$PSF_{NIRSPEC}$')
    ax1.set_ylabel("\(Flux (\cdot 10^{-18}\ erg\cdot s^{-1}\ \cdot cm^{-2}\ \cdot \mathring{{A}}^{-1})\)")
    ax1.set_xlim(np.min(x_values_as), np.max(x_values_as))
    ax2 = ax1.twiny()
    ax2.set_xlim(np.min(x_values), np.max(x_values))
    ax1.set_xlabel("\(Angular\ separation\ (arcsec)\)")
    ax2.set_xlabel("$Pixel\ separation\ (7x7mas\ config)$")
    ax1.legend(loc='upper center',ncols = 3, bbox_to_anchor = [0.5,1.4])
    ax2.text(-13, 4,rf"$\sigma_{{ELT}} = {popt_as:.3f}''$")
    ax2.text(-13, 3.5,rf"$\sigma_{{DL}} = {sigma_dl:.3f}''$")
    ax2.text(-13, 3.0,rf"$\sigma_{{JW}} = {sigma_JW:.3f}''$")
    plt.savefig("/home/steff/hsim/report_plots/spatial_PSF.jpg")
    plt.show()


    ### 2D CONTOUR
    X, Y = np.meshgrid(x_coords, y_coords)
    plt.figure()
    contour = plt.contourf(X, Y, image_collapsed * 1*10**(13), levels=50, cmap='inferno')
    plt.colorbar(contour, label='\(Intensity\)')
    plt.show()

    ### 3D CONTOUR ###

    mask_x = (X[0, :] > -10) & (X[0, :] < 10)
    mask_y = (Y[:, 0] > -10) & (Y[:, 0] < 10)

    X_cut = X[np.ix_(mask_y, mask_x)]
    Y_cut = Y[np.ix_(mask_y, mask_x)]
    image_cut = image_collapsed[np.ix_(mask_y, mask_x)]

    ax = plt.figure().add_subplot(projection='3d')
    ax.plot_surface(X_cut, Y_cut , image_cut * (1*10**(13)), cmap='viridis')
    ax.set_ylim(-10,10)
    ax.set_xlim(-10,10)

    z_min, z_max = 0, np.max(image_cut * (1*10**(13)))
    region_size = 6
    r = region_size / 2     # or any radius you want in plot units

    # vertical limits
    h_bottom = z_min
    h_top    = z_max

    # cylinder resolution
    theta = np.linspace(0, 2*np.pi, 120)
    z_vals = np.linspace(h_bottom, h_top, 60)
    
    Theta, Z = np.meshgrid(theta, z_vals)

    # cylinder coordinates
    X = r * np.cos(Theta)
    Y = r * np.sin(Theta)

    # plot only the side wall
    ax.plot_surface(
        X, Y, Z,
        color='red',
        alpha=0.3,
        zorder=3,
        linewidth=0
    )

    ax.set_ylabel("\(yaxis\ pixels\)", labelpad=15)
    ax.set_xlabel("\(xaxis\ pixels\)", labelpad=15)
    ax.set_zlabel("\(Flux\ (\cdot 10^{-13}\ erg\ \cdot\ s^{-1}\ \cdot\ cm^{-2}\ \cdot\ \mathring{{A}}^{-1})\)", labelpad=15)
    
    
    plt.show()


def gaussian_line(x,A, B, C):
    std = C
    B = B
    gaussian = (A/(np.sqrt(2*np.pi)*std))*np.exp(-0.5*((x-B)**2)/std**2)
    return gaussian


def spectral_psf(spectrum, spectrum_flux,data_flux, spectrum_std, wavelength_angstrom):
    mask = (wavelength_angstrom > (20800 - 50)) & (wavelength_angstrom < (20800 + 50))
    wavelength = wavelength_angstrom[mask]
    flux_data = spectrum_flux[mask]
    print(f"wavelength: {len(wavelength)} & flux: {len(flux_data)}")

    
    centre = 20800
    p0 = [np.nanmax(flux_data)*10**18, centre, 5]
    popt,pcov = curve_fit(gaussian_line, wavelength, flux_data*10**18,p0=p0 )
    plt.figure(figsize=(4,3.1))
    plt.step(wavelength, flux_data*10**(18))
    plt.plot(wavelength, gaussian_line(wavelength, *popt), c='red', ls='--')
    plt.ylabel("\(Flux (\cdot 10^{-18}\ erg\cdot s^{-1}\ \cdot cm^{-2}\ \cdot \mathring{{A}}^{-1})\)")
    plt.xlabel("$Wavelength (\mathring{A})$")
    plt.text(20750, 0.33, rf"$\sigma_{{HeII}} = {popt[2]:.2f}$")
    plt.savefig("/home/steff/hsim/report_plots/spectral_PSF.png", bbox_inches="tight", dpi=300)
    plt.show()
    



## IMPORT DATA

output_mass = 7.0
output_exp_time = 360

data = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/data.npy")
spectrum = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/spectrum_counts.npy")
spectrum_flux = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/spectrum_flux.npy")
spectrum_std = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/spectrum_std.npy")
wavelength_angstrom = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/wavelength_angstrom.npy")
data_flux = np.load(f"{dir_basic}/{output_exp_time}ks_exposures/M_{output_mass}_output/data_flux.npy")


spatial_plots(spectrum, spectrum_flux,data_flux, spectrum_std, wavelength_angstrom)
spectral_psf(spectrum, spectrum_flux,data_flux, spectrum_std, wavelength_angstrom)
