from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table, fits, inset_axes, mark_inset,patches, curve_fit, GridSpec, axes3d, CircularAperture,aperture_photometry
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

    x_values = np.linspace(-20,20,100)
    popt, pcov = curve_fit(gaussian_fit, x_coords, x_cut)
    #print(popt)
    
    
    lambda_he = 23210*10**(-10)
    D = 39
    FWHM = 1.2*(lambda_he/D)
    print(f"FWHM: {FWHM}")
    FWHM_pixels = FWHM*206265 / (7*10**(-3))
    gaussian_seeing = np.max(x_cut)*np.exp(-0.5*((x_values)**2)/(FWHM_pixels/(2*np.sqrt(2*np.log(2))))**2)

    pixel_to_arcsec = 7*10**(-3)

    arcsecs = (x_values * pixel_to_arcsec)
    
        
    fig, ax1 = plt.subplots()
    ax1.plot(x_coords, x_cut*10**(13), label='\(Observed\ flux\)')
    ax1.plot(x_values, gaussian_fit(x_values, *popt)*10**(13), c='orange', label=f'\(Gaussian\ fit:\ \sigma\ =\)$ {popt[0]:.1f}$')
    ax1.plot(x_values, gaussian_seeing*10**(13), ls='--', c='green', label= '\(Diffraction\ limited\ PSF\)')
    ax1.set_ylabel("\(Flux (\cdot 10^{-13}\ erg\cdot s^{-1}\ \cdot cm^{-2}\ \cdot \mathring{{A}}^{-1})\)")
    ax1.set_xlim(np.min(x_values), np.max(x_values))
    ax2 = ax1.twiny()
    ax2.set_xlim(np.min(arcsecs), np.max(arcsecs))
    ax2.set_xlabel("\(Angular\ separation\ (arcsec)\)")
    ax1.set_xlabel("\(Pixels\ from\ source\)")
    ax1.legend(loc='upper left')
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
