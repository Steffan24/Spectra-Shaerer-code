from modules import np, plt, Normalize, ascii
from variables import n, save, ttt, imf, mup, low, sfh, n_single, n_array
from import_spectrum import import_data
from constants import c_m
import os
import matplotlib.cm as cm

imf_array = ['sal', 'sca', 'logA', 'logB', 'logE', 'l05']


def import_EW(ttt,imf,mup,low,sfh):
    file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.22"
    print(file_loc)
    if os.path.exists(file_loc):
        data = ascii.read(file_loc,guess = True, data_start = 0)
        age_log = data['col1']
        HeII_1640 = data['col14']
        return age_log, HeII_1640
    else:
        return None
    
    

def IMF_EW(imf_array):
    colors = cm.rainbow(np.linspace(0, 1, len(imf_array)))
    plt.figure()
    for i in range(len(imf_array)):
        result  = import_EW(ttt, imf_array[i], mup, low, sfh)
        if result is None:
            print(f"skipping IMF: {imf_array[i]}")
            continue
        age_log, HeII_1640 = result
        plt.plot(age_log, HeII_1640,c = colors[i], label = f'\(IMF = {imf_array[i]}\)')
    plt.ylabel('\(EW (\mathring{A})\)')
    plt.xlabel('\(\log_{10}{t} (yr)\)')
    plt.legend()
    plt.show()
        
IMF_EW(imf_array) 
    

age_log, HeII_1640_EW = import_EW(ttt,imf,mup,low,sfh)

SED_data = import_data(n,save)
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
        plt.scatter(wavelength, flux_HeII_1640, label = f"age = {age_log[j]}")
    plt.show()

line_fit(HeII_1640, wavelength, c_m, age_log)
    
