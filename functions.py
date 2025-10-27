# functions.py
# list of functions used

from modules import np, plt, ScalarMappable, Normalize, ascii, latex, os, cosmo, interpolate, mticker, Table
from constants import T_sun, c_m, T_100M, M_sun_kg, G, kb, c, h, pc, AU, d, R_sun, M_sun
import plotting_params

def import_data(n, save, ttt, imf, mup, low, sfh, n_single, n_array):
    #read data and determine if single or multiple files
    if n == 'single':
        file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_single}"
        if os.path.exists(file_loc):
            age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
            age_data = np.array(age['6'][0])
            data = ascii.read(file_loc,guess = True, data_start = 0)
            wavelength = np.array(data['col1'])
            flux_raw = np.array(data['col3'])
            SED_flux = flux_raw #/ (4*np.pi*(d**2))
            age_repeat = np.repeat(age_data, len(wavelength))
            SED_data = {"ages": np.array(age_repeat),
                        "wavelengths": np.array(wavelength),
                        "SED_flux": np.array(SED_flux)}
    elif n == 'multi':
        all_ages = []
        all_wavelengths = []
        all_fluxes = []
        for i in range(len(n_array)):
            file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_array[i]}"
            if os.path.exists(file_loc):
                age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
                age_data = np.array(age['6'][0])
                data = ascii.read(file_loc,guess = True, data_start = 2)
                wavelength = np.array(data['col1'])
                all_wavelengths.append(wavelength)
                flux_raw = np.array(data['col3'])
                SED_flux = flux_raw #/ (4*np.pi*(d**2))
                age_repeat = np.repeat(age_data, len(wavelength))
                all_fluxes.append(SED_flux)
                all_ages.append(age_repeat)
        SED_data = {"ages": np.array(all_ages),
                    "wavelengths": all_wavelengths,
                    "SED_flux": all_fluxes}
    if save == True:
        np.save(f'{n}_saved_data_model_{ttt}_{imf}_{mup}_{low}_{sfh}', SED_data)
    return SED_data

def import_opacity():
    file_loc = "/home/steff/hsim/zackrisson_pop3_all/chtrans_nir_18504ft_lon30d_pwv0900um_zenith_smoothed2e-3um-1.dat"
    data = ascii.read(file_loc, guess = True, data_start = 2)
    wavelength_opacity = np.array(data['Wavelength'])
    wavelength_opacity_angstroms = wavelength_opacity * (10**4)
    transmittance_opacity = np.array(data['Transmittance'])
    opacity_data = {"wavelength": wavelength_opacity_angstroms,
                    "transmittance": transmittance_opacity}
    return opacity_data

def import_OH():
    file_loc = "/home/steff/hsim/zackrisson_pop3_all/skylineOH.txt"
    data = ascii.read(file_loc, guess = True, data_start = 2)
    wavelength_1 = data["wave1"]
    wavelength_2 = data["wave2"]
    wavelength_full = np.concatenate((wavelength_1, wavelength_2))
    flux_1 = data["flux"]
    flux = np.concatenate((flux_1, flux_1))
    sort_indices = np.argsort(wavelength_full)
    wavelength_full = wavelength_full[sort_indices]
    flux = flux[sort_indices]
    skyline_data = {"wavelength": np.array(wavelength_full),
                    "flux": np.array(flux)}
    return skyline_data
    

def plot(n, SED_data, lambda_sun, B, log_B, lambda_blackbody, blackbody, log_blackbody):
    if n == 'single':
        log_SED = np.log10(SED_data["SED_flux"])
        log_wavelength = np.log10(SED_data["wavelengths"])
        plt.figure()
        plt.plot(np.log10(lambda_sun), log_B - 30, label = "\(Blackbody\ 1M_{\odot}\)", linestyle = '--', color = 'red', zorder = 3)
        plt.plot(log_wavelength, log_SED - 30,c='darkblue',label = f'\({SED_data["ages"][0]} Myr\ since\ ZAMS\)', zorder = 1)
        #plt.xlim(0,7500)
        plt.ylim(0,8)
        plt.xlim(2,5.3)
        plt.xlabel("\(\log{\lambda}\ (\mathring{A})\)")
        plt.ylabel("\(logF_{\lambda}\ 1e+10\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
        plt.legend(bbox_to_anchor = [0.85,1], ncols = 2)
        plt.show()
        
    elif n == 'multi':
        fig,(ax1,ax2) = plt.subplots(1,2,width_ratios=[0.95,0.05])
        ax1.plot(np.log10(lambda_sun), log_B - 30, label = '\(Blackbody\ 1M_{\odot}\)', linestyle = '--', color='red', zorder = 3)
        ax1.plot(np.log10(lambda_blackbody), log_blackbody - 30, label = '\(Blackbody\ 100M_{\odot}\)', linestyle = '--', color = 'orange', zorder = 2)
        t_100 = np.log10((10**10)*(100)**(-2.5))
        t_50 = np.log10((10**10)*(50)**(-2.5))
        t_10 = np.log10((10**10)*(10)**(-2.5))
        t_5 = np.log10((10**10)*(5)**(-2.5))
        cmap = 'viridis'
        ages = np.array([a[0] for a in SED_data["ages"]])
        #norm = Normalize(np.min(ages), np.max(ages))
        norm = Normalize(vmin=ages.min(), vmax=ages.max())
        cmap = plt.cm.viridis
        colours = cmap(norm(ages))
        bar = ScalarMappable(cmap=cmap, norm=norm)
        for i in range(len(SED_data["ages"])):
            log_wavelength = np.log10(SED_data["wavelengths"][i])
            log_SED = np.log10(SED_data["SED_flux"][i])
            age = SED_data["ages"][i][0]
            ax1.plot(log_wavelength, log_SED - 30, c=colours[i])
        ax1.set_xlabel("\(\log{\lambda}\ (\mathring{A})\)")
        ax1.set_ylabel("\(logF_{\lambda}\ 1e+10\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot M_{\odot}^{-1})\)")
        ax1.set_ylim(-1.8,7)
        ax1.set_xlim(2,5.2)
        ax1.legend(loc='upper right')
        ax1.text(4.776,-0.38 + 2.6, '\(t_{100M_{\odot}}\)', bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy=(4.383,0.7), xycoords='data', xytext=(4.776,-0.38 + 2.6), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(4.4,0.2 + 2.6, '\(t_{50M_{\odot}}\)', bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy=(4.024,1.14), xycoords='data', xytext=(4.387,0.32 + 2.6), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(3.994,2.68 + 2.55, '\(t_{10M_{\odot}}\)', bbox=dict(edgecolor='black', fc='None'))
        ax1.annotate("",xy = (3.363,2.26), xycoords='data', xytext=(3.994,2.68 + 2.55), textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(4.186, 1.28 + 2.55, '\(t_{5M_{\odot}}\)', bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy = (3.252,1.65), xycoords='data', xytext=(4.186,1.28 + 2.55), textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle='arc3'))
        plt.colorbar(bar,cax=ax2,location = 'right', orientation = 'vertical')
        labels = [f'\({min(ages)}\)','\(t_{100M_{\odot}}\)','\(t_{50M_{\odot}}\)','\(t_{10M_{\odot}}\)','\(t_{5M_{\odot}}\)', f'\({max(ages)}\)']
        ax2.set_yticks([min(ages),t_100, t_50, t_10, t_5,max(ages)], labels=labels)
        ax2.set_ylabel('\(\log{t}\ since\ ZAMS\ (yr)\)')
        plt.show()
                     


def sun_type_star():
    lambda_sun = np.linspace(0, 300000 , 1*10**7) #angstrom
    lambda_sun_cm = lambda_sun * 1*10**(-8)
    B = ((2*h*c**2)/(lambda_sun_cm**5))*(1/(np.exp(h*c/(lambda_sun_cm*kb*T_sun)) - 1))
    B = B * (1*10**(-8))
    B = 4*np.pi * B * (R_sun)**2# /((d**2))
    log_B = np.log10(B)
    return lambda_sun, B, log_B



def blackbody(T):
    lambda_blackbody = np.linspace(0, 300000 , 1*10**5) #angstrom
    lambda_blackbody_m = lambda_blackbody * 1*10**(-8)
    B = ((2*h*c**2)/(lambda_blackbody_m**5))*(1/(np.exp(h*c/(lambda_blackbody_m*kb*T)) - 1))
    blackbody = B * (1*10**(-8))
    blackbody = 4*np.pi * blackbody* (R_sun*(13.8))**2#/ ((d**2)*100)
    log_blackbody = np.log10(blackbody)
    return lambda_blackbody, blackbody, log_blackbody

def interpolate_SED(SED_data, n, z, R):
    if n == 'single':
        mask = (SED_data["wavelengths"] < 8000)
        y = SED_data["SED_flux"][mask]
        x = SED_data["wavelengths"][mask]
        resolution = 1500/R
        print(f"Resolution: {resolution}")
        sampling = resolution/(5)
        print(f"Sampling pre redshift: {sampling}")
        x_new = np.arange(min(x), max(x), sampling)
        interpolated_wavelengths = interpolate.interp1d(x,y)
        interpolated_fluxes = interpolated_wavelengths(x_new)
        #print(f"INTERPOLATED: {interpolated_fluxes}")
        plt.figure()
        plt.plot(x, y, c='green')
        plt.plot(x_new, interpolated_fluxes, c='blue')
        plt.show()
        SED_data = {"ages": SED_data["ages"],
                    "wavelengths": x_new,
                    "SED_flux": np.array(interpolated_fluxes)}
        return SED_data
    
                    
    

def import_lines(ttt, imf, mup, low, sfh):
    #read data
    file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.22"
    if os.path.exists(file_loc):
        data = ascii.read(file_loc,guess = True, data_start = 0)
       # print(data['col15'])
        age_log = data['col1']
        H_beta = data['col5']
        H_lya = data['col7'] * H_beta
        H_alpha = data['col9'] * H_beta
        H_beta_ = data['col11'] * H_beta
        HeI_4471 = data['col13'] * H_beta
        HeII_1640 = data['col15'] * H_beta
        HeII_4686 = data['col17'] * H_beta
        HeII_3203 = data['col19'] * H_beta
        HeII_4541 = data['col21'] * H_beta
    return age_log, H_beta, H_lya, H_alpha, H_beta_, HeI_4471, HeII_1640, HeII_4686, HeII_3203, HeII_4541

def gaussian_profile(M, R, wavelength, age_log, H_beta, H_lya, H_alpha, H_beta_, HeI_4471, HeII_1640, HeII_4686, HeII_3203, HeII_4541):
    print(f"INPUT PARAMS: M = {M/1000} kg, R = {R/100} m")
    print(f"LEN WAVELENGTH: {len(wavelength)}")
    sigma_gal = 30000 #np.sqrt((M/1000)*G/(R/100)) # m/s (or 30000)
    print(f"SIGMA GAL: {sigma_gal} m/s")
    #array in peak [H_beta, H_lya, H_alpha, HEI_4471, HeII1640, HeII_4686, HeII_3203, HeII_4541]]
    lambda_peak_array = [4861, 1215, 6563, 4471, 1640, 4686, 3203, 4541]
    lambda_peak_array = np.array(lambda_peak_array)
    sigma_line = (sigma_gal/c_m) * np.array(lambda_peak_array)
    print(f"SIGMA LINE: {sigma_line} Angstrom")
    H_beta_peak_Intensity = (H_beta[0])#/(4*np.pi*(d**2))
    #print(f"H_beta_peak_Intensity: {np.log10(H_beta_peak_Intensity) + 30}")
    H_lya_peak_intensity = (H_lya[0])#/(4*np.pi*(d**2))
    #print(f"H_lya_peak_Intensity: {np.log10(H_lya_peak_intensity) + 30}")
    H_alpha_peak_intensity = (H_alpha[0])#/(4*np.pi*(d**2))
    #print(f"H_alpha_peak_Intensity: {np.log10(H_alpha_peak_intensity) + 30}")
    HeI_4471_peak_intensity = (HeI_4471[0])#/(4*np.pi*(d**2))
    #print(f"H_4471_peak_Intensity: {np.log10(HeI_4471_peak_intensity) + 30}")
    HeII_1640_peak_intensity = (HeII_1640[0])#/(4*np.pi*(d**2))
    print(f"H_1640_peak_Intensity: {np.log10(HeII_1640_peak_intensity) - 30}")
    HeII_3203_peak_intensity = (HeII_3203[0])#/(4*np.pi*(d**2))
    #print(f"H_3203_peak_Intensity: {np.log10(HeII_3203_peak_intensity) + 30}")
    HeII_4541_peak_intensity = (HeII_4541[0])#/(4*np.pi*(d**2))
    #print(f"H_4541_peak_Intensity: {np.log10(HeII_4541_peak_intensity) + 30}")
    HeII_4686_peak_intensity = (HeII_4686[0])#/(4*np.pi*(d**2))
    #print(f"H_4686_peak_Intensity: {np.log10(HeII_4686_peak_intensity) + 30}")
    flux_H_beta = []
    flux_H_lya = []
    flux_H_alpha = []
    flux_H_4471 = []
    flux_HII_1640 = []
    flux_HII_4686 = []
    flux_HII_3203 = []
    flux_HII_4541 = []
    for i in range(len(wavelength)):
        exponential_1 = np.exp(((wavelength[i] - lambda_peak_array[0])**2)/(-2*(sigma_line[0])**2))
        flux_H_beta_1 = ((H_beta_peak_Intensity/(np.sqrt(2*np.pi)*sigma_line[0])) * exponential_1)
        flux_H_beta.append(flux_H_beta_1)
        exponential_2 = np.exp(((wavelength[i] - lambda_peak_array[1])**2)/(-2*(sigma_line[1])**2))
        flux_H_lya_1 = ((H_lya_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[1]))*exponential_2)
        flux_H_lya.append(flux_H_lya_1)
        exponential_3 = np.exp(((wavelength[i] - lambda_peak_array[2])**2)/(-2*(sigma_line[2])**2))
        flux_H_alpha_1 = ((H_alpha_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[2]))*exponential_3)
        flux_H_alpha.append(flux_H_alpha_1)
        exponential_4 = np.exp(((wavelength[i] - lambda_peak_array[3])**2)/(-2*(sigma_line[3])**2))
        flux_H_4471_1 = ((HeI_4471_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[3]))*exponential_4)
        flux_H_4471.append(flux_H_4471_1)
        exponential_5 = np.exp(((wavelength[i] - lambda_peak_array[4])**2)/(-2*(sigma_line[4])**2))
        flux_H_1640_1 = ((HeII_1640_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[4])) * exponential_5)
        flux_HII_1640.append(flux_H_1640_1)
        exponential_6 = np.exp(((wavelength[i] - lambda_peak_array[5])**2)/(-2*(sigma_line[5])**2))
        flux_H_4686_1 = ((HeII_4686_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[5]))*exponential_6)
        flux_HII_4686.append(flux_H_4686_1)
        exponential_7 = np.exp(((wavelength[i] - lambda_peak_array[6])**2)/(-2*(sigma_line[6])**2))
        flux_H_3203_1 = ((HeII_3203_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[6]))*exponential_7)
        flux_HII_3203.append(flux_H_3203_1)
        exponential_8 = np.exp(((wavelength[i] - lambda_peak_array[7])**2)/(-2*(sigma_line[7])**2))
        flux_H_4541_1 = ((HeII_4541_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[7]))*exponential_8)
        flux_HII_4541.append(flux_H_4541_1)
    flux_H_beta = np.array(flux_H_beta)
    flux_H_lya = np.array(flux_H_lya)
    flux_H_alpha = np.array(flux_H_alpha)
    flux_H_4471 = np.array(flux_H_4471)
    flux_HII_1640 = np.array(flux_HII_1640)
    flux_HII_4686 = np.array(flux_HII_4686)
    flux_HII_3203 = np.array(flux_HII_3203)
    flux_HII_4541 = np.array(flux_HII_4541)

    plt.figure()
    plt.plot(wavelength, np.log10(flux_HII_1640) - 30)
    plt.plot(wavelength, np.log10(flux_H_lya) - 30)
    plt.ylim(0,7)
    plt.show()


    sigma_gal = [30000, 150000, 300000]
    data_ = Table()
    for i in range(len(sigma_gal)):
        sigma_line = (sigma_gal[i]/c_m) * np.array(lambda_peak_array)
        print(f"SIGMA LINE: {sigma_line} Angstrom")
        H_beta_peak_Intensity = (H_beta[0])#/(4*np.pi*(d**2))
        HeII_1640_peak_intensity = (HeII_1640[0])#/(4*np.pi*(d**2))
        print(f"H_1640_peak_Intensity: {HeII_1640_peak_intensity}")
        flux_HII_1640 = []
        for j in range(len(wavelength)):
            exponential_5 = np.exp(((wavelength[j] - lambda_peak_array[4])**2)/(-2*(sigma_line[4])**2))
            flux_H_1640_1 = ((HeII_1640_peak_intensity/(np.sqrt(2*np.pi)*sigma_line[4])) * exponential_5)
            flux_HII_1640.append(flux_H_1640_1)
        data_[f'sigma_{sigma_gal[i]}']= np.array(flux_HII_1640)
    data_['wavelength'] = wavelength
    ascii.write(data_, f'sigma_check_F_{HeII_1640_peak_intensity}.dat', format='basic', overwrite=True)
    print("saved file")
    data_ = Table.read('sigma_check.dat', format='ascii.basic')
    plt.figure(figsize=(8,5))
    plt.plot(wavelength, data_['sigma_30000'])
    plt.plot(wavelength, data_['sigma_150000'])
    plt.plot(wavelength, data_['sigma_300000'])
    plt.xlabel('Wavelength [Å]', fontsize=12)
    plt.ylabel('Flux [arbitrary units]', fontsize=12)
    plt.title('Gaussian Line Profiles for Different Velocity Dispersions', fontsize=13)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    

    


    
    return flux_H_beta, flux_H_lya, flux_H_alpha, flux_H_4471, flux_HII_1640, flux_HII_3203, flux_HII_4541, flux_HII_4686


def full_spectra(n, wavelength, flux_H_beta, flux_H_lya, flux_H_alpha, flux_H_4471, flux_HII_1640, flux_HII_3203, flux_HII_4541, flux_HII_4686, SED_data):
    if n == 'single':
        total_flux_lines = []
        for i in range(len(wavelength)):
            flux = (SED_data["SED_flux"][i]) + flux_H_beta[i] + flux_H_lya[i] + flux_H_alpha[i] + flux_H_4471[i] + flux_HII_1640[i] + flux_HII_3203[i] + flux_HII_4541[i] + flux_HII_4686[i]
            total_flux_lines.append(flux)
        total_flux_lines = np.array(total_flux_lines)
    if n == 'multi':
        total_flux_lines = []
        for i in range(len(SED_data["ages"])):
            wavelength = (SED_data["wavelengths"][i])
            SED_fluxes = (SED_data["SED_flux"][i])
            age = SED_data["ages"][i][0]
            flux = SED_fluxes + flux_H_beta[i] + flux_H_lya[i] + flux_H_alpha[i] + flux_H_4471[i] + flux_HII_1640[i] + flux_HII_3203[i] + flux_HII_4541[i] + flux_HII_4686[i]
            total_flux_lines.append(flux)
        total_flux_lines = np.array(total_flux_lines)
    return total_flux_lines
    

def plot_full_spectra(n, total_flux_lines, SED_data, wavelength_z, z):
    if n == 'single':
        wavelength = SED_data["wavelengths"]
        fig, ax1 = plt.subplots()
        ax1.plot(np.log10(wavelength), np.log10(total_flux_lines) - 30, lw = 2)
        ax1.set_ylabel("logs")
        ax1.set_xlim(2.5, 4.5)
        ax1.set_ylim(2,9)
        ax1.set_xlabel("\(\log{\lambda}\ (\mathring{A})\)")
        ax1.set_ylabel("\(logF_{\lambda}\ 1e+10\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
        ax1.text(3.696,0.542 + 2.55, r'\(H\)$\beta$', bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy=(3.689,1.337 + 2.55), xycoords='data', xytext=(3.729,0.852 + 2.55), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(2.85,5.23 + 2.55,r'\(Ly-\)$\alpha$',bbox=dict(edgecolor='black', fc = 'None'))
        ax1.text(3.856,2.795 + 2.55,r'\(H\)$\alpha$',bbox=dict(edgecolor='black', fc = 'None'))
        ax1.text(3.4,1 + 2.55,'\(HII_{4471}\)',bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy=(3.645,1.395 + 2.55), xycoords='data', xytext=(3.575,1.014 + 2.55), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(3.136,4.35 + 2.6,'\(HeII_{1640}\)',bbox=dict(edgecolor='black', fc = 'None'))
        #ax1.annotate("",xy=(3.210, 3.176 + 2.7), xycoords='data', xytext=(3.216,4.227 + 2.55), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(3.3,5.792,'\(HII_{3203}\)',bbox=dict(edgecolor='black', fc = 'None'))
        ax1.text(3.58,3.38 + 2.55,'\(HII_{4541}\)',bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy=(3.668,2.36 + 2.55), xycoords='data', xytext=(3.659,3.275 + 2.55), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(3.488,-0.06 + 2.55,'\(HII_{4686}\)',bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy=(3.657,1.37 + 2.55), xycoords='data', xytext=(3.5,-0.06 + 2.55), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
        ax2 = ax1.twiny()
        ax2.plot(np.log10(wavelength_z), total_flux_lines, alpha = 0)
        ax2.set_xlim(3.54, 5.54)
        ax2.set_xlabel(
    rf"\textbf{{$\log \lambda\, (\mathrm{{\AA}})\ redshift\ z = {z}$}}"
)
        NIR = np.log10(np.array([0.8, 2.5])*(10**4))
        ax2.plot(NIR, [8.661, 8.661], c = 'red', lw = 3)
        ax2.text(4.443,8.5, '\(:\ HARMONI\ range\)', bbox=dict(edgecolor='white', fc = 'None'))
        plt.show()
    elif n == 'multi':
        fig,(ax1,ax2) = plt.subplots(1,2,width_ratios=[0.95,0.05])
        all_ages = np.array([a[0] for a in SED_data["ages"]])
        norm = Normalize(vmin=all_ages.min(), vmax=all_ages.max())
        cmap = plt.cm.viridis
        colours = cmap(norm(all_ages))
        bar = ScalarMappable(cmap=cmap, norm=norm)
        for i in range(len(all_ages)):
            wavelength = SED_data["wavelengths"][i]
            flux = total_flux_lines[i]
            ax1.plot(np.log10(wavelength), np.log10(flux) + 10, color=colours[i], lw = 2)

        ax1.set_xlim(2, 5.5)
        ax1.set_ylim(-2, 7)
        ax1.set_xlabel(r"$\log{\lambda}\ (\mathring{A})$")
        ax1.set_ylabel(r"$\log{F_{\lambda}}\ 1e+10\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})$")

        cbar = plt.colorbar(bar, cax=ax2, orientation='vertical')
        cbar.set_label(r"$\log{t}\ \mathrm{since\ ZAMS\ (yr)}$")

        plt.tight_layout()
        plt.show()
            
            
def redshifting(n, total_flux_lines, SED_data, z):
     if n == 'single':
         comoving_d = cosmo.comoving_distance(z)
         print(f"comoving_d: {comoving_d} Mpc")
         d_lumo_cm = (1 + z)* (comoving_d * (3.086*10**24))

         wavelength = SED_data["wavelengths"]

         flux_z = total_flux_lines * (1/(d_lumo_cm**2)) * (1/(1 + z)) 
         wavelength_z = wavelength * (1 + z)

         flux_z = flux_z.value

         return flux_z, wavelength_z



        
def AB_magnitude_conversion(flux_z, wavelength_z):
    flux_zab = -2.5*np.log10(flux_z) - 2.402 - 5.0*np.log10(wavelength_z)
    return flux_zab

def magnitudes_to_flux(mag, wavelength_z):
    return 10 ** (-(mag + 2.402 + 5.0*np.log10(wavelength))/2.5)

def import_harmoni_res(LR_IZJ_min, LR_IZJ_max,LR_HK_min,LR_HK_max,MR_IZ_min,MR_IZ_max,MR_J_min,MR_J_max,MR_H_min,MR_H_max,MR_K_min,MR_K_max):
    LR_IZJ = np.linspace(LR_IZJ_min, LR_IZJ_max,3)* 10**4
    LR_HK = np.linspace(LR_HK_min, LR_HK_max,3)* 10**4
    MR_IZ = np.linspace(MR_IZ_min, MR_IZ_max,3)* 10**4
    MR_J = np.linspace(MR_J_min, MR_J_max,3)* 10**4
    MR_H = np.linspace(MR_H_min, MR_H_max,3)* 10**4
    MR_K = np.linspace(MR_K_min, MR_K_max,3)* 10**4
    return LR_IZJ, LR_HK, MR_IZ, MR_J, MR_H, MR_K

def plot_spectra_redshifted(flux_z, wavelength_z, flux_zab, LR_IZJ, LR_HK, MR_IZ, MR_J, MR_H, MR_K, opacity_data, skyline_data, z, R):
    fig, (ax1, ax3, ax4) = plt.subplots(3,1, height_ratios=[3,1,1], sharex = True)
    ax4.set_xlabel("\(\lambda (\mathring{A})\)")
    ax1.set_ylabel("\(\log{f_{\lambda}} (erg/s/cm^2/\mathring{A}/M_{\odot})\)", labelpad = 12)
    ax1.plot((wavelength_z), np.log10(flux_z))


    def AB_magnitude_conversion_single(log_flux, wavelength_ref):
        flux = 10**log_flux
        return -2.5*np.log10(flux) - 2.402 - 5.0*np.log10(wavelength_ref)

    def magnitudes_to_flux_single(mag, wavelength_ref):
        flux = 10 ** (-(mag + 2.402 + 5.0*np.log10(wavelength_ref))/2.5)
        return np.log10(flux)

    mean_lambda = 17000
    print(f"mean wavelength: {mean_lambda}")

    def flux_to_mag(y):
        return AB_magnitude_conversion_single((y), mean_lambda)

    def mag_to_flux(y):
        return magnitudes_to_flux_single((y), mean_lambda)

    ax1.set_ylim(-29.5, -24)

    ax2 = ax1.secondary_yaxis('right', functions=(flux_to_mag, mag_to_flux))

    ax2.set_ylabel(r"$AB\ Magnitude$")

    print(LR_IZJ)

    ax1.plot(LR_IZJ, np.array([-28.5, -28.5, -28.5]), linewidth = 3, c='blue')
    ax1.text(10680,-28.377,'\(IZJ\)',bbox=dict(edgecolor='None', fc = 'None'), fontsize = 18)
    ax1.plot(LR_HK, np.array([-28.5, -28.5, -28.5]), linewidth = 3, c= 'red')
    ax1.text(19040,-28.377,'\(HK\)',bbox=dict(edgecolor='None', fc = 'None'), fontsize = 18)
    ax1.plot(MR_IZ, np.array([-28.75, -28.75, -28.75]), linewidth = 3, c ='green')
    ax1.text(9150,-29,'\(IZ\)',bbox=dict(edgecolor='None', fc = 'None'), fontsize = 18)
    ax1.plot(MR_J, np.array([-28.75, -28.75, -28.75]), linewidth = 3, c='purple')
    ax1.text(11750,-29,'\(J\)',bbox=dict(edgecolor='None', fc = 'None'), fontsize = 18)
    ax1.plot(MR_H, np.array([-28.75, -28.75, -28.75]), linewidth = 3, c='blue')
    ax1.text(16110,-29,'\(H\)',bbox=dict(edgecolor='None', fc = 'None'), fontsize = 18)
    ax1.plot(MR_K, np.array([-28.75, -28.75, -28.75]), linewidth = 3, c = 'red')
    ax1.text(22000,-29,'\(K\)',bbox=dict(edgecolor='None', fc = 'None'), fontsize = 18)
    ax1.text(13500, -25.790, r'\(Ly-\)$\alpha$', bbox=dict(edgecolor='black', fc='None'))
    ax1.text(17580, -26, r'\(HeII_{1640}\)', bbox=dict(edgecolor='black', fc='None'))
   # angstrom_Ly = (angstrom_redshifted/18010) * 13350
    #ax1.text(21200, -26, rf"$\begin{{array}}{{c}} " \
    #    rf"Z = {z} \\ R = {R} \\ \sigma_{{He_{{II1640}}}} = {angstrom_redshifted:.2f}\,\mathring{{A}} \\ \sigma_{{Ly-\alpha}} = {angstrom_Ly:.2f}\,\mathring{{A}}" \
     #   rf"\end{{array}}$", bbox=dict(edgecolor='white', fc='None'))
    NIR = (np.array([0.8, 2.5])*(10**4))
    ax1.set_xlim(np.min(NIR),np.max(NIR))
    flux = skyline_data["flux"]#*(10**-7)
    norm_flux = (flux) #*(10**(-10))
    for i in range(len(skyline_data['wavelength'])):
                   f = flux[i]
                   w = skyline_data['wavelength'][i]
                   ax3.vlines(w, 0, f, color='blue', alpha=0.6, linewidth=0.8)
    ax4.plot((opacity_data["wavelength"]), opacity_data["transmittance"])
    #ax3.plot((skyline_data["wavelength"]),norm_flux)
    #ax5 = ax3.secondary_yaxis('right', functions=(flux_to_mag, mag_to_flux))
    ax4.set_ylim(0.25,1.1)
    ax4.set_ylabel("\(T (\%)\)", labelpad = 6)
    ax3.set_ylabel("\(f\)", labelpad = -5)
    pos1 = ax1.get_position()
    pos3 = ax3.get_position()
    pos4 = ax4.get_position()

    #ax1.spines['bottom'].set_visible(False)
    #ax3.spines['top'].set_visible(False)
    
    y_solid = (pos1.y0 + pos3.y1) / 2

    fig.add_artist(plt.Line2D([0.125, 0.9], [y_solid, y_solid],
                              transform=fig.transFigure, color='black', lw=1.5, linestyle = '--'))

    plt.show()



    fig, (ax1, ax3, ax4) = plt.subplots(3,1, height_ratios=[3,1,1], sharex = True)
    ax4.set_xlabel("\(\lambda (\mathring{A})\)")
    ax1.set_ylabel("\(\log{f_{\lambda}} (erg/s/cm^2/\mathring{A})\)")
    ax1.plot((wavelength_z), np.log10((flux_z*(10**8))))
    
    ax1.set_ylim(-23, -19)

    ax2 = ax1.secondary_yaxis('right', functions=(flux_to_mag, mag_to_flux))

    ax2.set_ylabel(r"$AB\ Magnitude$")

    print(LR_IZJ)

    NIR = (np.array([0.8, 2.5])*(10**4))
    ax1.set_xlim(np.min(NIR),np.max(NIR))
    flux = skyline_data["flux"]#*(10**3)
    norm_flux = (flux) #*(10**(-10))
    for i in range(len(skyline_data['wavelength'])):
                   f = flux[i]
                   w = skyline_data['wavelength'][i]
                   ax3.vlines(w, 0, f, color='blue', alpha=0.6, linewidth=0.8)
    ax4.plot((opacity_data["wavelength"]), opacity_data["transmittance"])
    #ax3.plot((skyline_data["wavelength"]), norm_flux)
    #ax5 = ax3.secondary_yaxis('right', functions=(flux_to_mag, mag_to_flux))
    ax4.set_ylim(0.25,1.05)
    ax4.set_ylabel("\(T (\%)\)")
    ax3.set_ylabel("\(f\)")
    pos1 = ax1.get_position()
    pos3 = ax3.get_position()
    pos4 = ax4.get_position()

    ax1.spines['bottom'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    
    y_solid = (pos1.y0 + pos3.y1) / 2

    ax1.axhline(mag_to_flux(27.2), 0,1, c='red', ls = '--', label = '\(10mas\ limiting\ magnitude\)')
    ax1.axhline(mag_to_flux(26.3), 0,1, c='orange', ls = '--', label = '\(4mas\ limiting\ magnitude\)')
    ax1.legend()

    fig.add_artist(plt.Line2D([0.125, 0.9], [y_solid, y_solid],
                              transform=fig.transFigure, color='black', lw=1.5, linestyle = '--'))

    plt.show()
    

