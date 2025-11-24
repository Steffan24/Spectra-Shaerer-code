from functions import collapse_cube_data, extracting_over_aperture
from variables import output_file, output_flux, output_std, output_scale, dir_basic, output_mass
from modules import np, plt, curve_fit

spectrum = np.load(f"{dir_basic}/M_{output_mass}_output/spectrum_counts.npy")
spectrum_flux = np.load(f"{dir_basic}/M_{output_mass}_output/spectrum_flux.npy")
spectrum_std = np.load(f"{dir_basic}/M_{output_mass}_output/spectrum_std.npy")
wavelength_angstrom = np.load(f"{dir_basic}/M_{output_mass}_output/wavelength_angstrom.npy")

def gaussian(x,A):
    std = 3.244
    B = 23205
    gaussian = (A/(np.sqrt(2*np.pi)*std))*np.exp(-0.5*((x-B)**2)/std**2)
    return gaussian

def chi_squared_gauss(x,y, std, popt):
    chi = ((y - gaussian(x, *popt))**2)/(std**2)
    return np.sum(chi)

def chi_squared_cont(x,y,std, straight_fit):
    chi = ((y - straight_fit)**2) / (std**2)
    #print(chi)
    return np.nansum(chi)

def calcSNR_chi2(wavelength_angstrom,spectrum, spectrum_std):
    mask_gauss = (wavelength_angstrom > 23000) & (wavelength_angstrom < 23300)
    wavelength_gauss = wavelength_angstrom[mask_gauss]
    counts_gauss = spectrum[mask_gauss]
    std_gauss = np.array(spectrum_std[mask_gauss])
    print(spectrum_std.shape)
    popt, pcov, infodict, mesg, ier = curve_fit(gaussian, wavelength_gauss,counts_gauss,sigma = std_gauss,absolute_sigma = True, p0=[1], full_output = True)
    print(mesg)
    print(ier)
    chi_squared_gaussian = chi_squared_gauss(wavelength_gauss,counts_gauss, std_gauss, popt)
    print(f"chi squared gauss: {chi_squared_gaussian}")
    straight_fit = np.ones(len(wavelength_gauss)) * np.nanmedian(np.array(counts_gauss))
    #print(f"straight fit: {straight_fit}")
    chi_squared_cont_value = chi_squared_cont(wavelength_gauss, counts_gauss, std_gauss, straight_fit)
    print(f"chi squared cont: {chi_squared_cont_value}")

    SNR = np.sqrt(np.abs(chi_squared_gaussian - chi_squared_cont_value))
    print(f"SNR_chi_method: {SNR}")

    plt.plot(wavelength_gauss, counts_gauss)
    plt.show()
    plt.figure()
    plt.plot(wavelength_gauss, counts_gauss)
    plt.plot(wavelength_gauss, straight_fit)
    plt.show()
    plt.figure()
    plt.plot(wavelength_angstrom, spectrum)
    plt.show()
    plt.figure()
    plt.plot(wavelength_gauss, gaussian(wavelength_gauss, *popt), c = 'red')
    plt.plot(wavelength_gauss, counts_gauss)
    plt.show()


def calcSNR_peak(wavelength_angstrom, spectrum, spectrum_std):
    mask_gauss = (wavelength_angstrom > 23190) & (wavelength_angstrom < 23220)
    wavelength_gauss = wavelength_angstrom[mask_gauss]
    counts_gauss = spectrum[mask_gauss]
    std_gauss = spectrum_std[mask_gauss]
    popt, pcov = curve_fit(gaussian, wavelength_gauss,counts_gauss, p0=[1])
    A = popt[0]
    peak = A / (np.sqrt(2*np.pi)*3.24)
    std = np.nanstd(spectrum[wavelength_angstrom > 23300])
    return peak / std*np.sqrt(2.355*(3.244/3.3))

print(f"SNR method 2: {calcSNR_peak(wavelength_angstrom, spectrum, spectrum_std)}")
    

calcSNR_chi2(wavelength_angstrom,spectrum, spectrum_std)

    
    
