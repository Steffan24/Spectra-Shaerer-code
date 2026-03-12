from modules import np,plt, ascii, latex, os, cosmo, interpolate
import plotting_params
from scipy import datasets
from scipy.ndimage import gaussian_filter

## variables
ttt = 'ge0' 
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'
n_single = 31

z_gal = 10.6
input_scale = 1*10**(-3)

## import data

file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_single}"
data = ascii. read(file_loc, guess = True, data_start = 0)
wavelength = np.array(data['col1'])
flux = np.array(data['col3'])
cont_data = {"wavelength": np.array(wavelength),
             "flux": np.array(flux)}

# interpolate to harmoni

GN_z11_resolution = 0.0779 # at z=0
wavelength_h = np.arange(min(wavelength), max(wavelength), GN_z11_resolution)
interpolated = interpolate.interp1d(cont_data["wavelength"], cont_data["flux"])
flux_h = interpolated(wavelength_h)

## import lines

lines_file = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.22"
data = ascii.read(lines_file,guess = True, data_start = 0)
print(data.colnames)
age_log = data['col1']
H_beta = data['col5']
H_lya = data['col7'] * H_beta
HeII_1640 = data['col15'] * H_beta

# insert as gaussian

def insert_one(sigma_gal, z_gal, input_scale):
    lambda_lines = np.array([1215, 1640])
    c = 3*10**8 #m/s
    sigma_line = (sigma_gal/c) * lambda_lines #m/s
    H_lya_peak = (H_lya[0])
    HeII_1640_peak = (HeII_1640[0])

    flux_H_lya = (H_lya_peak/(np.sqrt(2*np.pi)*sigma_line[0])) * np.exp((wavelength_h - lambda_lines[0])**2 / (-2*sigma_line[0]**2))
    flux_HeII_1640 = (HeII_1640_peak/(np.sqrt(2*np.pi)*sigma_line[1])) * np.exp((wavelength_h - lambda_lines[1])**2 / (-2*sigma_line[1]**2))

    total_flux = []
    for i in range(len(wavelength_h)):
        flux_add = flux_h[i] + flux_H_lya[i] + flux_HeII_1640[i]
        total_flux.append(flux_add)

    total_flux = np.array(total_flux)

    plt.figure()
    plt.plot(wavelength_h, total_flux)
    plt.xlabel("$wavelength (\mathring{A})$")
    plt.ylabel("$Flux (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot M_{\odot}^{-1})$")
    plt.show()

    comoving_d = cosmo.comoving_distance(z_gal)
    print(f"comoving dist: {comoving_d}Mpc")
    d_lumo_cm = (1+z_gal) * (comoving_d * (3.086*10**24))

    wavelength = wavelength_h * (1+z_gal)
    flux = total_flux * (1/(4*np.pi*d_lumo_cm**2)) * (1/(1+z_gal))

    flux_pas = flux / (input_scale**2)

    wavelength = (wavelength/(1+z_gal)) * (1+11.7)

    mask_wavelength = (wavelength > 15000) & (wavelength < 22500)
    wavelength = np.array(wavelength[mask_wavelength])
    flux_pas = np.array(flux_pas[mask_wavelength])

    return wavelength, flux_pas

# import GN-z11 spectra

spectrum_wavelength = np.load('/home/steff/hsim/GNz11/model_spectra/wavelength.npy')
spectrum_flux = np.load('/home/steff/hsim/GNz11/model_spectra/flux.npy')
spectrum_flux_pa = spectrum_flux / (input_scale**2)
spectrum_wavelengt = (spectrum_wavelength / (1+z_gal)) * (1+11.7)
mask_wavelength = (spectrum_wavelengt > 15000) & (spectrum_wavelengt < 22500)
spectrum_wavelength = spectrum_wavelengt[mask_wavelength]
spectrum_flux_pas = spectrum_flux_pa[mask_wavelength]

### create data cube

cube_length = 0.2

n_spaxels = int(round(cube_length / input_scale))
n_wavelength = len(spectrum_wavelength)
min_wavelength = min(spectrum_wavelength)
central_spaxel = n_spaxels // 2
spectral_res = spectrum_wavelength[26] - spectrum_wavelength[25]
print(f"Number of wavelength points: {n_wavelength}")
print(f"Number of spaxels: {n_spaxels}")
print(f"Minimum wavelength: {min_wavelength}")
print(f"central spaxel no: {central_spaxel}")
print(f"spectral resolution: {spectral_res}")

bytes_cube = n_wavelength * n_spaxels * n_spaxels * 4  # float32
print(f"cube size    = {bytes_cube/1e9:.2f} GB")
print("=======================")

print("Creating input cube")
cube = np.zeros((n_wavelength, n_spaxels, n_spaxels), dtype=np.float32)
print("created")
# create the resolved GN-z11 gal
print("Adding pt source")
cube[:, central_spaxel, central_spaxel] = spectrum_flux_pas.astype(np.float32)
print("Added")
half_light_radius = 0.016 #arcsec
resolved_fwhm = half_light_radius*2 #arcsec
sigma = resolved_fwhm / (2*np.sqrt(2*np.log(2))) # arcsec
#convert to pixels using spatial resolution (angstrom/pixel)
spatial_res = 0.001 #arcsec/pixel
sigma_pixel = sigma/spatial_res
print(f"gaussian sigma: {sigma_pixel}")

print("gaussian smooth")
new_cube = gaussian_filter(cube, sigma=(0, sigma_pixel, sigma_pixel)).astype(np.float32)
print("smoothed")
# add pop III clusters
print("Adding popIII clusters")
pop_file = '/home/steff/hsim/zackrisson_pop3_all/code/GN-z11_codes/GN-z11_model/popII_coords.txt'
data = ascii.read(pop_file,guess = True, data_start = 1)
mass = data["mass"]
coordx = data["coordx"]
coordy = data["coordy"]

for i in range(len(mass)):
    print(f"{mass[i]}")
    mass = np.array(mass[i])
    cluster_mass = (10**(mass))*(1.99*10**30)
    G = 6.67*10**(-11)
    r = 4 * (3.086*10**16)
    sigma_gal = np.sqrt(cluster_mass * G / r)
    print("Inserting")
    wavelength, flux_pas_popIII = insert_one(sigma_gal, z_gal, input_scale)
    current_flux = new_cube[:, coordx[i], coordy[i]]
    flux_pas_combine = current_flux + flux_pas_popIII

    new_cube[:, coordx[i], coordy[i]] = flux_pas_combine.astype(np.float32)
    print("Done")

print("Writing file")
hdu = fits.PrimaryHDU(data=new_cube) # create cube
#header adjustments
header = hdu.header
header.clear()
header["SIMPLE"] = bool(SIMPLE)
header["BITPIX"] = bool(BITPIX)
header["NAXIS"] = NAXIS
header["NAXIS1"] = cube.shape[2]
header["NAXIS2"] = cube.shape[1]
header["NAXIS3"] = cube.shape[0]
header["EXTEND"] = bool(EXTEND)
header["CTYPE1"] = CTYPE1
header["CTYPE2"] = CTYPE2
header["CTYPE3"] = CTYPE3
header["CUNIT1"] = CUNIT1
header["CUNIT2"] = CUNIT2
header["CUNIT3"] = CUNIT3
header["CDELT1"] = CDELT1
header["CDELT2"] = CDELT2
header["CRVAL3"] = min_wavelength
header["CDELT3"] = spectral_res
header["CRPIX3"] = CRPIX3
header["BUNIT"] = BUNIT
header["SPECRES"] = 3*spectral_res
hdul = fits.HDUList([hdu])
hdul.writeto(f'/home/steff/hsim/GNz11/model_spectra/gauss_source_input.fits', overwrite = True)
print("Input cube saved")
    
    
    







    

    




    
    

                                
