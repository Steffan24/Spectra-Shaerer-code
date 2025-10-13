### MODULES ###
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from astropy.io import ascii
import latex
import os

### PLOTTING PARAMS ###
plt.rcParams["font.size"] = 30
plt.rcParams["legend.fontsize"] = 30
plt.rcParams["figure.frameon"] = False
plt.rcParams["figure.titlesize"] = 16
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.subplot.wspace"] = 0.03
plt.rcParams["legend.frameon"] = False
plt.rcParams["axes.linewidth"] = 1.75

### CONSTANTS USED ###
T_sun = 5772 #k
kb = 1.38*10**(-16)
c = 3*10**10
c_m = 3*10**8
h = 6.63*10**(-27)
T_100M = 10**(4.7)
pc = 3.08*10**18
AU = 1.496*10**13
d = 1 * AU
R_sun = 6.957*10**10
M_sun = 1.989*10**33
M_sun_kg = 1.989*10**30
G = 6.67*10**(-11)


### VARIABLES (unchanged as of 07/10/25)  ###
ttt = 'ge0' #Options - ge0: no mass loss, mdt: with mass loss
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'


n_single = 31 #SED of population (single plotting)

n_array = np.linspace(31,1030, 999, dtype='int') # SED (multiple plotting)


### FUNCTION FOR SINGLE LINE PLOT ####

def single_plotting(n_single):
    #read data
    
    file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_single}"
    if os.path.exists(file_loc):
        age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
        print(age)
        age_data = age['6'][0]
        print(f"AGE: {age_data}")
        age_data = (10**(age_data))*10**(-6)
        data = ascii.read(file_loc,guess = True, data_start = 2)
        wavelength = data['col1']
        total_flux = data['col3']
        total_flux = total_flux / (4*np.pi*(d**2)*M_sun)
        log_flux = np.log10(total_flux)
        log_flux = log_flux

        #PLOT 1: logf vs lambda
        plt.figure()
        plt.plot(lambda_sun, log_B+30, label = "\(Blackbody\ 1M_{\odot}\)", linestyle = '--', color = 'red', zorder = 2)
        plt.plot(wavelength, log_flux+30,c='darkblue',label = f'\({age_data} Myr\ since\ ZAMS\)', zorder = 1)
        #plt.xlim(0,7500)
        #plt.ylim(0,8)
        plt.xlim(0,6500)
        plt.xlabel("\(\lambda (\mathring{A})\)")
        plt.ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
        plt.legend(bbox_to_anchor = [0.85,1], ncols = 2)
        plt.show()

        #PLOT 2: logf vs loglambda
        plt.figure()
        plt.plot(np.log10(lambda_sun), log_B+ 30, label = "\(Blackbody\ 1M_{\odot}\)", linestyle = '--', color = 'red', zorder = 2)
        plt.plot(np.log10(wavelength), log_flux + 30,c='darkblue',label = f'\({age_data} Myr\ since\ ZAMS\)', zorder = 1)
        #plt.xlim(0,7500)
        plt.ylim(0,8)
        plt.xlim(2,5.3)
        plt.xlabel("\(\log{\lambda}\ (\mathring{A})\)")
        plt.ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
        plt.legend(bbox_to_anchor = [0.85,1], ncols = 2)
        plt.show()

### FUNCTION FOR MANY PLOTS ###

def multiple_plotting(n_array):
    #PLOT 1: logf vs lambda
    fig,(ax1,ax2) = plt.subplots(1,2,width_ratios=[0.95,0.05])
    ax1.plot(lambda_sun, log_B + 30, label = '\(Blackbody\ 1M_{\odot}\)', linestyle = '--', color='red', zorder = 2)
    colour = plt.cm.viridis(np.linspace(0,1,len(n_array)))
    cmap = 'viridis'
    norm = Normalize(0,max(n_array))
    bar = ScalarMappable(cmap=cmap, norm=norm)
    all_ages = []
    for i in range(len(n_array)):
        #read dara
        file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_array[i]}"
        if os.path.exists(file_loc):
            age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
            age_data = age['6'][0]
            all_ages.append(age_data)
            data = ascii.read(file_loc,guess = True, data_start = 2)
            wavelength = data['col1']
            total_flux = data['col3']
            total_flux = total_flux / (4*np.pi*(d**2)*M_sun)
            log_flux = np.log10(total_flux)
            ax1.plot(wavelength, log_flux + 30,c=colour[i], zorder = 1)
    all_ages_raw = []
    for i in range(len(all_ages)):
        all_ages_raw.append((10**(-6))*10**(all_ages[i]))
    ax1.set_xlim(0,7500)
    ax1.set_xlabel("\(\lambda (\mathring{A})\)")
    ax1.set_ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
    ax1.set_ylim(0,8)
    ax1.legend(bbox_to_anchor = [1.15,1.11])
    plt.colorbar(bar,cax=ax2,location = 'right', orientation = 'vertical')
    labels = [f'\({min(all_ages_raw)}\)', f'\({max(all_ages_raw)}\)']
    ax2.set_yticks([0,1030], labels=labels)
    ax2.set_ylabel('\(Time\ since\ ZAMS\ (Myr)\)')
    plt.show()

    #PLOT 2: f vs lambda
    fig,(ax1,ax2) = plt.subplots(1,2,width_ratios=[0.95,0.05])
    ax1.plot(np.log10(lambda_sun), log_B + 30, label = '\(Blackbody\ 1M_{\odot}\)', linestyle = '--', color='red', zorder = 2)
    ax1.plot(np.log10(lambda_blackbody), log_blackbody + 30, label = '\(Blackbody\ 100M_{\odot}\)', linestyle = '--', color = 'orange', zorder = 2)
    
    all_ages = []
    for i in range(len(n_array)):
        #read dara
        file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_array[i]}"
        if os.path.exists(file_loc):
            age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
            age_data = age['6'][0]
            all_ages.append(age_data)
    colour = plt.cm.viridis(np.linspace(0,1,len(all_ages)))
    cmap = 'viridis'
    norm = Normalize(min(all_ages),max(all_ages))
    bar = ScalarMappable(cmap=cmap, norm=norm)
    for i in range(len(all_ages)):
        #read dara
        file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_array[i]}"
        if os.path.exists(file_loc):
            data = ascii.read(file_loc,guess = True, data_start = 2)
            wavelength = data['col1']
            total_flux = data['col3']
            total_flux = total_flux / (4*np.pi*(d**2)*M_sun)
            log_flux = np.log10(total_flux)
            ax1.plot(np.log10(wavelength), log_flux + 30 ,c=colour[i], zorder = 1)
    all_ages_raw = []
    for i in range(len(all_ages)):
        all_ages_raw.append((10**(-6))*10**(all_ages[i]))
    #ax1.set_xlim(0,7500)
    ax1.set_xlabel("\(\log{\lambda}\ (\mathring{A})\)")
    ax1.set_ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
    ax1.set_ylim(-2,6)
    ax1.set_xlim(2,5.2)
    ax1.legend(bbox_to_anchor = [1,1])
    ax1.text(4.4,0.2, '\(0.01 Myr\)', bbox=dict(edgecolor='black', fc = 'None'))
    ax1.annotate("",xy=(4.286,0.18), xycoords='data', xytext=(4.387,0.32), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
    plt.colorbar(bar,cax=ax2,location = 'right', orientation = 'vertical')
    #lifetimes
    t_100 = np.log10((10**10)*(100)**(-2.5))
    t_50 = np.log10((10**10)*(50)**(-2.5))
    t_10 = np.log10((10**10)*(10)**(-2.5))
    t_5 = np.log10((10**10)*(5)**(-2.5))
    labels = [f'\({min(all_ages)}\)','\(t_{100M_{\odot}}\)','\(t_{50M_{\odot}}\)','\(t_{10M_{\odot}}\)','\(t_{5M_{\odot}}\)', f'\({max(all_ages)}\)']
    ax2.set_yticks([min(all_ages),t_100, t_50, t_10, t_5,max(all_ages)], labels=labels)
    ax2.set_ylabel('\(\log{t}\ since\ ZAMS\ (yr)\)')
    plt.show()
    return all_ages
        
        
### FUNCTION FOR SUN TYPE STAR ###
def sun_type_star():
    lambda_sun = np.linspace(0, 300000 , 1*10**5) #angstrom
    lambda_sun_cm = lambda_sun * 1*10**(-8)
    B = ((2*h*c**2)/(lambda_sun_cm**5))*(1/(np.exp(h*c/(lambda_sun_cm*kb*T_sun)) - 1))
    B = B * (1*10**(-8))
    B = 4*np.pi * B * (R_sun)**2 /((d**2)*M_sun)
    log_B = np.log10(B)
    return lambda_sun, B, log_B

### FUNCTION FOR A MODEL BLACKBODY ###
def blackbody(T):
    lambda_blackbody = np.linspace(0, 300000 , 1*10**5) #angstrom
    lambda_blackbody_m = lambda_blackbody * 1*10**(-8)
    B = ((2*h*c**2)/(lambda_blackbody_m**5))*(1/(np.exp(h*c/(lambda_blackbody_m*kb*T)) - 1))
    blackbody = B * (1*10**(-8))
    blackbody = 4*np.pi * blackbody* (R_sun*(13.8))**2/ ((d**2)*100*M_sun)
    log_blackbody = np.log10(blackbody)
    return lambda_blackbody, blackbody, log_blackbody
    
#### IMPORT RECOMBINATION LINES ##

def recombination():
    file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.22"
    if os.path.exists(file_loc):
        data = ascii.read(file_loc,guess = True, data_start = 2)
        print(data)
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

def gaussian_profile(M, R, age_1):
    sigma_gal =np.sqrt(M*G/R)
    print(sigma_gal)
    #Hbeta
    lambda_data = np.linspace(0, 6000 , 1*10**6)
    lambda_peak_H_beta = 4861 #angstrom
    sigma_line = (sigma_gal/c_m) * lambda_peak_H_beta
    print(sigma_line)
    H_beta_peak_Intensity = (H_beta[0])/(4*np.pi*(d**2)*M_sun)
    print(H_beta_peak_Intensity)
    flux_H_beta = (H_beta_peak_Intensity/(np.sqrt(2*np.pi)*sigma_line)*np.exp((-0.5)*((np.subtract(lambda_data, lambda_peak_H_beta))**2)/(sigma_line**2)))
    print(flux_H_beta)
    plt.figure()
    plt.plot(lambda_data, flux_H_beta)
    plt.show()
            
    
    

def single_plot_with_recomb():
    log_wavelength, log_flux = single_plotting(n_single)
    

#lambda_sun, B, log_B = sun_type_star()
#lambda_blackbody, blackbody, log_blackbody = blackbody(T_100M)
        
#all_ages = multiple_plotting(n_array)        
age_log, H_beta, H_lya, H_alpha, H_beta_, HeI_4471, HeII_1640, HeII_4686, HeII_3203, HeII_4541 = recombination()        
age_1 = age_log[0]
gaussian_profile(M_sun_kg*(10**6) , 100*pc/100, age_1)

#single_plotting(n_single)



