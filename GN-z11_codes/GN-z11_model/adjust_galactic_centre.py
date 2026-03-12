
import plotting_params
from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table, fits, inset_axes, mark_inset,patches, curve_fit, GridSpec, axes3d, CircularAperture,aperture_photometry
from collapse_pop_output import import_data, isolate_popIII, analyse_lines

output_dir = '/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/salpeter_run/salpeter_run/run_6'

counts_file = f'{output_dir}/run_6_reduced.fits'
flux_file = f'{output_dir}/run_6_reduced_flux_cal.fits'
std_file = f'{output_dir}/run_6_std.fits'

extraction_region = 36
xcoord = 114
ycoord =114
centre = 100
d = 14
mass = 6 #solar mass
output_scale = 7*10**(-3)
region_size = 2

counts_data, counts_wavelength, flux_data, flux_wavelength, std_data, std_wavelength = import_data(counts_file, flux_file, std_file)

delta_centre = np.linspace(0,15,16)
print(f"delta center: {delta_centre}")

for i in range(len(delta_centre)):
    dir_save = f'/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/SNR_arrays_move_cen/delta_centre_{delta_centre[i]}.npy'
    dir_save_image = f'/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/SNR_arrays_move_cen/spectra_{delta_centre[i]}.png'
    new_centre = centre - 2*delta_centre[i]
    flux_region, flux_nearby, popIII_flux, counts_region, counts_nearby, popIII_counts, std_region, std_nearby, popIII_std = isolate_popIII(flux_file, counts_file, new_centre, d, output_scale, std_file, region_size)

    analyse_lines(flux_region, flux_nearby, popIII_flux, counts_region, counts_nearby, popIII_counts, mass, dir_save, dir_save_image)
    #plot collapsed image with circ
    hdul = fits.open(counts_file)
    data = hdul[0].data
    collapsed_image = np.nansum(data, axis=0)
    collapsed_image = np.nan_to_num(collapsed_image)
    ny, nx = collapsed_image.shape
    pixel_scale_mas = 7
    extent = [0, nx*pixel_scale_mas, 0, ny*pixel_scale_mas]
    fig, ax = plt.subplots()
    image = ax.imshow(collapsed_image, origin='lower', cmap='gray_r', aspect='equal', extent=extent)
    circ = patches.Circle((114,114),region_size, alpha=0.8, fc='None',ec='red', label = '$Cluster$')
    circ_ref = patches.Circle((new_centre - d,new_centre -d), region_size, alpha=0.8, fc='None', ec='blue', label='$Ref$')
    ax.add_patch(circ)
    ax.add_patch(circ_ref)
    ax.set_xlabel('$x\ (mas)$')
    ax.set_ylabel('$y\ (mas)$')
    ax.legend(loc='upper center',ncols = 2, bbox_to_anchor=(0.5, 1.2), frameon=False)
    plt.savefig(f"/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/SNR_arrays_move_cen/aperture_{new_centre}.png")
    plt.show()



SNR_HeII = []
SNR_NIV = []
SNR_CIV = []
SNR_OIII = []
SNR_NIII = []
SNR_Ly_a = []

for j in range(len(delta_centre)):
    dir_save = f'/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/SNR_arrays_move_cen/delta_centre_{delta_centre[j]}.npy'
    SNR_values = np.load(f"{dir_save}")
    SNR_HeII.append(SNR_values[3])
    SNR_NIV.append(SNR_values[1])
    SNR_CIV.append(SNR_values[2])
    SNR_OIII.append(SNR_values[4])
    SNR_NIII.append(SNR_values[5])


delta_centre_arcsec = delta_centre*(output_scale/7)
delta_centre_parsec = (delta_centre_arcsec / 0.016) * 64
print(f"delta centre parsec: {delta_centre_parsec}")
plt.figure()
plt.plot(delta_centre_arcsec, SNR_HeII, label = '$He_{II}$')
plt.plot(delta_centre_arcsec, SNR_NIV, label = '$NIV$')
plt.plot(delta_centre_arcsec, SNR_CIV, label = '$CIV$')
plt.plot(delta_centre_arcsec, SNR_OIII, label = '$OIII$')
plt.plot(delta_centre_arcsec, SNR_NIII, label = '$NIII$')
plt.hlines(y=5, xmin=0, xmax=max(delta_centre_arcsec) + 0.0002, colors = 'red',ls='--')
#plt.vlines(x = 0.0063, ymin=0, ymax = 15, colors='red', ls='--')
plt.text(0.0005,5.2,"No detection")
plt.xlim(0, max(delta_centre_arcsec) + 0.0002)
plt.ylim(0, max(SNR_HeII) + 1)
plt.xlabel("$\sigma_{gal, cen}\ (arcsec)$")
plt.ylabel("$SNR_{popIII}$")
plt.legend(loc='upper center',ncols = 3, bbox_to_anchor=(0.5, 1.35), frameon=False)
#plt.text(0.0075, 23, "$IMF: Salpeter$\n$M: 10^{7}M_{\odot}$", bbox=dict(boxstyle='round', facecolor='white',alpha=0.8))
plt.savefig('/home/steff/hsim/HSIM/hsim/output_cubes/GN-z11_outputs/SNR_arrays_move_cen/figure.png')
plt.show()
    
    
    
    

