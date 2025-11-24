from modules import np, plt, Normalize
from variables import n, save, ttt, imf, mup, low, sfh, n_single, n_array
from functions import import_lines, import_data
from constants import c_m

age_log, H_beta, H_lya, H_alpha, H_beta_, HeI_4471, HeII_1640, HeII_4686, HeII_3203, HeII_4541 = import_lines(ttt,imf,mup,low,sfh)

HeII_1640_EW = HeII_1640 / H_beta

SED_data = import_data(n,save,ttt,imf,mup,low,sfh,n_single,n_array)
wavelength = SED_data["wavelengths"]


plt.figure()
plt.plot(age_log, HeII_1640_EW)
plt.ylabel(f"\(EW\)")
plt.xlabel("\(\log_{10}{t} (yr)\)")
plt.show()

def line_fit(HeII_1640, wavelength, c_m, age_log):
    lambda_peak = 1640
    sigma_gal = 30000
    sigma_line = (sigma_gal/c_m) * lambda_peak
    plt.figure()
    for j in range(len(HeII_1640)):
        flux_HeII_1640 = []
        for i in range(len(wavelength)):
            exp = np.exp(((wavelength[i] - lambda_peak)**2)/(-2*(sigma_line)**2))
            flux = (HeII_1640 / (np.sqrt(2*np.pi)*sigma_line))*exp
            flux_HeII_1640.append(flux)
        plt.plot(wavelength, flux_HeII_1640, label = f"age = {age_log[j]}")
    plt.show()

line_fit(HeII_1640, wavelength, c_m, age_log)
    
